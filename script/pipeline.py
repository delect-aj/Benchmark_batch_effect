"""
MicroCoRE Framework -- Main Pipeline
Chains all four stages and exposes a unified end-to-end API.

corr_mode can be a single string ('binary', 'abundance', 'fusion'),
a list of strings, or the special value 'all'.  When multiple modes are
requested, each mode runs an independent pipeline (separate consensus matrix
and graph autoencoder), enabling fair comparison.
"""

import numpy as np

from correlation_builder import (build_batch_correlation_matrices,
                                  build_all_batch_correlation_matrices)
from consensus_fusion import build_consensus_matrix
from representation_learner import train_graph_autoencoder
from sample_aggregator import compute_sample_embeddings

_VALID_MODES = {"binary", "abundance", "fusion"}


def _normalise_modes(corr_mode):
    """Return a list of mode strings from whatever the user passed."""
    if corr_mode == "all":
        return ["binary", "abundance", "fusion"]
    if isinstance(corr_mode, str):
        return [corr_mode]
    modes = list(corr_mode)
    unknown = set(modes) - _VALID_MODES
    if unknown:
        raise ValueError(f"Unknown corr_mode(s): {unknown}. "
                         f"Choose from {_VALID_MODES}.")
    return modes


class MicroCoREModel:
    """
    End-to-end MicroCoRE pipeline with optional multi-mode comparison.

    Single-mode usage (original API, unchanged)
    -------------------------------------------
    model = MicroCoREModel(corr_mode="fusion", ...)
    model.fit(batch_data_list)
    sample_embeddings, batch_labels = model.transform(batch_data_list)
    # model.taxon_embeddings  -> np.ndarray (n_taxa, embed_dim)
    # model.consensus_matrix  -> np.ndarray (n_taxa, n_taxa)

    Multi-mode usage (new)
    ----------------------
    model = MicroCoREModel(corr_mode="all", ...)   # or a list
    model.fit(batch_data_list)
    embeddings_dict, batch_labels = model.transform(batch_data_list)
    # embeddings_dict = {"binary": ..., "abundance": ..., "fusion": ...}
    # model.taxon_embeddings_per_mode  -> dict of np.ndarray
    # model.consensus_matrix_per_mode  -> dict of np.ndarray
    # model.final_train_loss_per_mode  -> dict of float
    """

    def __init__(self,
                 embed_dim=64,
                 corr_mode="fusion",
                 binary_threshold=5,
                 fusion_weights=(0.4, 0.6),
                 cv_threshold=0.5,
                 min_abs_corr=0.1,
                 sparsity_threshold=0.1,
                 n_epochs=2000,
                 lr=0.01,
                 smoothness_weight=0.1,
                 weight_mode="rank",
                 device="cpu",
                 random_seed=42):
        self.embed_dim = embed_dim
        self.corr_mode = corr_mode          # original value, kept for inspection
        self._modes = _normalise_modes(corr_mode)
        self.binary_threshold = binary_threshold
        self.fusion_weights = fusion_weights
        self.cv_threshold = cv_threshold
        self.min_abs_corr = min_abs_corr
        self.sparsity_threshold = sparsity_threshold
        self.n_epochs = n_epochs
        self.lr = lr
        self.smoothness_weight = smoothness_weight
        self.weight_mode = weight_mode
        self.device = device
        self.random_seed = random_seed

        # Populated after fit(); always dicts keyed by mode string
        self.taxon_embeddings_per_mode = {}
        self.consensus_matrix_per_mode = {}
        self.final_train_loss_per_mode = {}

    # ------------------------------------------------------------------
    # Convenience properties for single-mode backward compatibility
    # ------------------------------------------------------------------
    @property
    def taxon_embeddings(self):
        if len(self._modes) == 1:
            return self.taxon_embeddings_per_mode.get(self._modes[0])
        return self.taxon_embeddings_per_mode

    @property
    def consensus_matrix(self):
        if len(self._modes) == 1:
            return self.consensus_matrix_per_mode.get(self._modes[0])
        return self.consensus_matrix_per_mode

    @property
    def final_train_loss(self):
        if len(self._modes) == 1:
            return self.final_train_loss_per_mode.get(self._modes[0])
        return self.final_train_loss_per_mode

    # ------------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------------
    def fit(self, batch_data_list, phylo_matrix=None):
        """
        Learn taxon embeddings for every requested corr_mode.

        Parameters
        ----------
        batch_data_list : list[np.ndarray]
            One count matrix per batch; all must share the same n_taxa.
        phylo_matrix : np.ndarray or None
            Optional phylogenetic proximity matrix, shape (n_taxa, n_taxa).
        """
        batch_sizes = [d.shape[0] for d in batch_data_list]

        # ------------------------------------------------------------------
        # Stage 1: per-batch correlation matrices
        # ------------------------------------------------------------------
        print("=" * 60)
        print("Stage 1: Within-batch correlation matrix construction")
        print("=" * 60)

        # When more than one mode is needed, compute phi and CLR once and
        # reuse them rather than repeating heavy per-batch computation.
        needs_multi = len(self._modes) > 1 or set(self._modes) == {"fusion"}
        if len(self._modes) > 1:
            all_corr = build_all_batch_correlation_matrices(
                batch_data_list,
                binary_threshold=self.binary_threshold,
                fusion_weights=self.fusion_weights)
        else:
            single_list = build_batch_correlation_matrices(
                batch_data_list,
                mode=self._modes[0],
                binary_threshold=self.binary_threshold,
                fusion_weights=self.fusion_weights)
            all_corr = {self._modes[0]: single_list}

        # ------------------------------------------------------------------
        # Stages 2 & 3: consensus + GAE, independently for each mode
        # ------------------------------------------------------------------
        for mode in self._modes:
            print("=" * 60)
            print(f"Stage 2 [{mode}]: Cross-batch consensus correlation fusion")
            print("=" * 60)
            consensus = build_consensus_matrix(
                all_corr[mode],
                batch_sizes,
                cv_threshold=self.cv_threshold,
                min_abs_corr=self.min_abs_corr,
                sparsity_threshold=self.sparsity_threshold)
            self.consensus_matrix_per_mode[mode] = consensus

            print("=" * 60)
            print(f"Stage 3 [{mode}]: Taxon graph embedding learning")
            print("=" * 60)
            embeddings, loss = train_graph_autoencoder(
                consensus,
                embed_dim=self.embed_dim,
                n_epochs=self.n_epochs,
                lr=self.lr,
                phylo_matrix=phylo_matrix,
                smoothness_weight=self.smoothness_weight,
                device=self.device,
                random_seed=self.random_seed)
            self.taxon_embeddings_per_mode[mode] = embeddings
            self.final_train_loss_per_mode[mode] = loss
            print(f"  [{mode}] taxon embedding shape: {embeddings.shape}, "
                  f"final loss: {loss:.6f}")

        return self

    def transform(self, batch_data_list):
        """
        Map raw count data to sample embeddings.

        Returns
        -------
        Single-mode  : (sample_embeddings np.ndarray, batch_labels np.ndarray)
        Multi-mode   : (dict {mode: sample_embeddings}, batch_labels np.ndarray)
        """
        if not self.taxon_embeddings_per_mode:
            raise RuntimeError("Call .fit() before .transform().")

        print("=" * 60)
        print("Stage 4: Sample-level representation aggregation")
        print("=" * 60)

        results = {}
        batch_labels = None
        for mode in self._modes:
            emb_list, batch_labels = compute_sample_embeddings(
                batch_data_list,
                self.taxon_embeddings_per_mode[mode],
                weight_mode=self.weight_mode)
            results[mode] = np.vstack(emb_list)

        if len(self._modes) == 1:
            return results[self._modes[0]], batch_labels
        return results, batch_labels

    def fit_transform(self, batch_data_list, phylo_matrix=None):
        """Convenience method: fit and transform in one call."""
        self.fit(batch_data_list, phylo_matrix=phylo_matrix)
        return self.transform(batch_data_list)
