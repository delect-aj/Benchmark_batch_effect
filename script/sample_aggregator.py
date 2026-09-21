"""
MicroCoRE Framework -- Stage 4
Sample-level representation aggregation

Core idea:
    Represent each sample as a (weighted) combination of taxon embeddings.
    The weighting function g(·) should be insensitive to multiplicative
    batch biases.

    Two aggregation strategies are provided:

    First-order aggregation ('rank', 'clr', 'binary', 'raw'):
        s_i = sum_k  g(O_{i,k}) * e_k
        Captures differences in abundance *levels*.
        'rank' is invariant to multiplicative biases within a batch.

    Second-order co-occurrence aggregation ('cooccurrence'):
        S_i = E^T diag(b_i) E   (flattened to a vector)
        Captures which embedding directions co-occur in sample i.
        This is the recommended mode when the phenotype signal resides
        primarily in co-variation structure rather than abundance levels.
"""

import numpy as np


def within_batch_rank_transform(abundance_matrix, normalize=True):
    """
    Replace each column (taxon) with within-batch sample ranks.

    Multiplicative bias does not change within-batch rank order, so this
    transform is invariant to multiplicative batch effects.

    Parameters
    ----------
    abundance_matrix : np.ndarray, shape (n_samples, n_taxa)
    normalize : bool
        If True, rescale ranks to [0, 1].

    Returns
    -------
    rank_matrix : np.ndarray
    """
    n_samples, n_taxa = abundance_matrix.shape
    rank_matrix = np.zeros_like(abundance_matrix, dtype=np.float64)
    for col in range(n_taxa):
        rank_matrix[:, col] = np.argsort(np.argsort(abundance_matrix[:, col]))
    if normalize and n_samples > 1:
        rank_matrix = rank_matrix / (n_samples - 1)
    return rank_matrix


def clr_transform(abundance_matrix, pseudocount=0.5):
    """
    Centered log-ratio (CLR) transformation.
    Provides partial robustness to multiplicative biases.
    """
    counts = abundance_matrix.astype(np.float64) + pseudocount
    row_sums = counts.sum(axis=1, keepdims=True)
    log_rel = np.log(counts / row_sums)
    return log_rel - log_rel.mean(axis=1, keepdims=True)


def binary_transform(abundance_matrix, threshold=5):
    """
    Binarise counts at a detection threshold.
    Maximally robust to multiplicative biases.
    """
    return (abundance_matrix >= threshold).astype(np.float64)


def compute_sample_embeddings(batch_data_list, taxon_embeddings,
                               weight_mode="rank"):
    """
    Aggregate taxon embeddings into one vector per sample.

    Parameters
    ----------
    batch_data_list : list[np.ndarray]
        One array per batch, each of shape (n_samples_in_batch, n_taxa).
    taxon_embeddings : np.ndarray, shape (n_taxa, embed_dim)
        Output of Stage 3.
    weight_mode : str
        'rank'        -- within-batch rank transform (invariant to
                         multiplicative bias; captures abundance-level signal)
        'clr'         -- CLR transform
        'binary'      -- binary presence/absence
        'raw'         -- relative abundance (no transformation)
        'cooccurrence' -- second-order: E^T diag(b_i) E, flattened.
                          Captures co-variation structure differences.
                          Recommended when phenotype signal is in co-variation
                          patterns rather than abundance levels.

    Returns
    -------
    sample_embedding_list : list[np.ndarray]
        One array per batch, shape (n_samples_in_batch, feature_dim).
    batch_labels : np.ndarray
        Integer array of length n_total_samples indicating batch membership.
    """
    valid_modes = {"rank", "clr", "binary", "raw", "cooccurrence"}
    if weight_mode not in valid_modes:
        raise ValueError(f"Unknown weight_mode: '{weight_mode}'. "
                         f"Choose from {valid_modes}.")

    sample_embedding_list = []
    batch_label_list = []

    for batch_idx, batch_data in enumerate(batch_data_list):
        n_samples, n_taxa = batch_data.shape
        embed_dim = taxon_embeddings.shape[1]

        if weight_mode == "cooccurrence":
            # Second-order representation:
            #   For sample i with binary presence vector b_i:
            #   S_i = E^T @ diag(b_i) @ E  (shape: embed_dim x embed_dim)
            #   Flattened to a vector of length embed_dim^2.
            #
            # This captures 'which embedding directions co-occur in sample i',
            # reflecting sample-specific co-variation patterns.
            binary_mat = (batch_data >= 1).astype(np.float64)
            feat_dim = embed_dim * embed_dim
            embeddings = np.zeros((n_samples, feat_dim), dtype=np.float64)

            for s in range(n_samples):
                # Weighted embedding: diag(b_i) @ E
                weighted_emb = binary_mat[s:s + 1].T * taxon_embeddings
                # Outer-product summary: E^T @ (diag(b_i) @ E)
                outer = taxon_embeddings.T @ weighted_emb
                embeddings[s] = outer.flatten()

            # Within-batch z-score normalisation to remove batch-level
            # intensity differences
            col_mean = embeddings.mean(axis=0, keepdims=True)
            col_std = embeddings.std(axis=0, keepdims=True)
            col_std[col_std < 1e-10] = 1e-10
            embeddings = (embeddings - col_mean) / col_std

        else:
            # First-order aggregation: s_i = weights_i @ E
            if weight_mode == "rank":
                weights = within_batch_rank_transform(batch_data, normalize=True)
            elif weight_mode == "clr":
                weights = clr_transform(batch_data)
            elif weight_mode == "binary":
                weights = binary_transform(batch_data)
            else:  # raw
                row_sums = batch_data.sum(axis=1, keepdims=True)
                row_sums[row_sums == 0] = 1.0
                weights = batch_data / row_sums

            embeddings = weights @ taxon_embeddings   # (n_samples, embed_dim)

            # L2 normalisation
            norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
            norms[norms < 1e-10] = 1e-10
            embeddings = embeddings / norms

        sample_embedding_list.append(embeddings)
        batch_label_list.extend([batch_idx] * n_samples)

    batch_labels = np.array(batch_label_list)
    total_samples = sum(d.shape[0] for d in batch_data_list)
    feat_dim = sample_embedding_list[0].shape[1]
    print(f"  Sample embeddings built -- mode: {weight_mode}, "
          f"total samples: {total_samples}, feature dim: {feat_dim}")

    return sample_embedding_list, batch_labels
