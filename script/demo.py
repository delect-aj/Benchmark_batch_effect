"""
MicroCoRE Framework -- Demo Script
Semisynthetic data generation + end-to-end pipeline + comparison with CLR baseline

Synthetic data design
---------------------
The demo uses a semisynthetic 'pure co-variation signal' scenario:
  - Background samples are bootstrapped from the Wirbel et al. 2019 CRC
    meta-analysis dataset (849 samples, 8 studies, ~300 species) to inherit
    realistic sparsity, taxon prevalence distributions, and taxon-taxon
    correlation structure.
  - A phenotype signal is injected as a mean-zero shared CLR-space factor
    on the signal taxa for positive samples.  This increases pairwise
    correlations without shifting marginal means.
  - Batch effects are applied as known multiplicative biases, per-sample
    log-space noise, and batch-specific detection-limit truncation.
  - Counts are drawn with Dirichlet-Multinomial sampling (overdispersed).

Methods compared
----------------
  CLR baseline   -- CLR transform + logistic regression (abundance-level)
  MicroCoRE/binary   -- phi coefficient correlation -> GAE
  MicroCoRE/abundance -- CLR-Pearson correlation -> GAE
  MicroCoRE/fusion   -- weighted phi+CLR-Pearson -> GAE   (recommended)
"""

import sys
import os
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from pipeline import MicroCoREModel

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


# -------------------------------------------------------------------------
# Reference data loading (Wirbel et al. 2019 CRC meta-analysis)
# -------------------------------------------------------------------------

def _fetch_reference_data(cache_dir):
    """
    Download and cache the Wirbel et al. 2019 CRC meta-analysis feature
    matrix and metadata from the zellerlab/crc_meta GitHub repository.
    """
    base_url = ("https://raw.githubusercontent.com/zellerlab/crc_meta"
                "/master/data/")
    files = {"feat": "feat_matrix.tsv", "meta": "meta_data.tsv"}
    paths = {}

    for key, fname in files.items():
        local = os.path.join(cache_dir, f"wirbel_{fname}")
        if not os.path.exists(local):
            url = base_url + fname
            print(f"  Downloading {fname} from zellerlab/crc_meta ...")
            try:
                urllib.request.urlretrieve(url, local)
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to download {url}.\n"
                    f"Please download manually and save to {local}.\n"
                    f"(Error: {exc})") from exc
        paths[key] = local

    return paths["feat"], paths["meta"]


def load_reference_data(cache_dir=None, n_taxa=120):
    """
    Load the Wirbel et al. 2019 CRC meta-analysis dataset as a simulation
    foundation.  Downloads and caches the files on the first call.

    The feature matrix is taxa × samples (relative abundances); the top
    n_taxa taxa by prevalence are selected.  Rows are renormalized to sum
    to one after subsetting.

    Parameters
    ----------
    n_taxa : int
        Number of taxa to retain (selected by cross-sample prevalence).

    Returns
    -------
    rel_abund : np.ndarray, shape (n_samples, n_taxa)
        Per-sample relative abundances (each row sums to 1).
    study_labels : np.ndarray, shape (n_samples,)
        Integer study IDs (0-indexed).
    """
    if cache_dir is None:
        cache_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "_ref_data")
    os.makedirs(cache_dir, exist_ok=True)

    feat_path, meta_path = _fetch_reference_data(cache_dir)

    feat_df = pd.read_csv(feat_path, sep="\t", index_col=0)  # taxa × samples
    meta_df = pd.read_csv(meta_path, sep="\t", index_col=0)  # samples × meta

    # Align common samples
    common = feat_df.columns.intersection(meta_df.index)
    feat_df = feat_df[common]
    meta_df = meta_df.loc[common]

    # Select top n_taxa by prevalence (fraction of samples where taxon > 0)
    prevalence = (feat_df > 0).mean(axis=1)
    top_taxa = prevalence.nlargest(n_taxa).index
    feat_df = feat_df.loc[top_taxa]

    # (n_samples, n_taxa) relative abundance matrix, renormalized
    rel_abund = feat_df.T.values.astype(np.float64)
    row_sums = rel_abund.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    rel_abund /= row_sums

    # Study labels — look for a 'Study' column in metadata (case-insensitive)
    study_col = next(
        (c for c in meta_df.columns if c.lower() == "study"), None)
    if study_col is None:
        study_labels = np.zeros(len(common), dtype=int)
    else:
        study_labels = np.array(pd.Categorical(meta_df[study_col]).codes)

    n_studies = len(np.unique(study_labels))
    print(f"  Reference data loaded: {rel_abund.shape[0]} samples, "
          f"{rel_abund.shape[1]} taxa, {n_studies} studies")
    print(f"  Zero rate in reference data: {(rel_abund == 0).mean():.3f}")

    return rel_abund, study_labels


# -------------------------------------------------------------------------
# Synthetic data generator
# -------------------------------------------------------------------------

def generate_synthetic_data(n_batches=4, n_samples_per_batch=80,
                             n_signal_taxa=15, signal_corr_increment=0.22,
                             bias_strength=2.0,
                             dm_concentration=0.1,
                             sample_noise_std=0.3,
                             sequencing_depth=10000,
                             n_taxa=120,
                             reference_data=None,
                             cache_dir=None,
                             random_seed=2024):
    """
    Semisynthetic microbiome count data generation.

    Background samples are bootstrapped from real microbiome data (Wirbel
    et al. 2019 CRC meta-analysis by default) to inherit realistic sparsity,
    taxon prevalence distributions, and taxon-taxon correlation structure.
    A phenotype signal is injected in CLR space; batch effects are applied
    as known multiplicative biases.

    Parameters
    ----------
    n_signal_taxa : int
        Number of taxa (selected by prevalence rank) that carry the
        phenotype signal.
    signal_corr_increment : float
        Magnitude of the shared CLR-space factor added to signal taxa for
        positive samples.  A zero-mean factor is used so that marginal CLR
        means are preserved; the variance increase per taxon is
        signal_corr_increment**2 (typically < 5 % of real CLR variance).
    dm_concentration : float
        Dirichlet-Multinomial concentration scale.  Total Dirichlet mass =
        dm_concentration * n_taxa.  Smaller → more overdispersed / sparse.
    sample_noise_std : float
        Std of per-(sample, taxon) log-space noise added on top of the
        shared batch bias (models individual technical variation).
    n_taxa : int
        Number of taxa to retain when loading the reference dataset.
        Ignored if reference_data is provided explicitly.
    reference_data : np.ndarray or None
        Real relative abundance matrix, shape (n_ref_samples, n_taxa).
        Each row must sum to 1.  If None, Wirbel et al. 2019 is loaded
        automatically via load_reference_data().
    cache_dir : str or None
        Directory for caching the downloaded reference dataset.

    Returns
    -------
    batch_data_list : list[np.ndarray]   observed count matrices per batch
    label_list      : list[np.ndarray]   binary phenotype labels per batch
    true_bias_list  : list[np.ndarray]   true multiplicative bias vectors
    """
    rng = np.random.RandomState(random_seed)

    # ------------------------------------------------------------------
    # Load or accept the real-data reference
    # ------------------------------------------------------------------
    if reference_data is None:
        reference_data, _ = load_reference_data(
            cache_dir=cache_dir, n_taxa=n_taxa)

    n_ref, n_taxa = reference_data.shape

    # Ensure rows are valid relative abundances
    row_sums = reference_data.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    ref_rel = reference_data / row_sums

    # ------------------------------------------------------------------
    # CLR helpers
    # ------------------------------------------------------------------
    def _clr(rel, pseudo=5e-4):
        x = rel + pseudo
        log_x = np.log(x)
        return log_x - log_x.mean(axis=1, keepdims=True)

    def _inv_clr(clr_mat):
        e = np.exp(clr_mat)
        return e / e.sum(axis=1, keepdims=True)

    # ------------------------------------------------------------------
    # Simulate batches
    # ------------------------------------------------------------------
    batch_data_list, label_list, true_bias_list = [], [], []

    for batch_idx in range(n_batches):
        n_pos = n_samples_per_batch // 2
        n_neg = n_samples_per_batch - n_pos
        labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
        rng.shuffle(labels)
        pos_idx = np.where(labels == 1)[0]

        # --- bootstrap real samples as background ---
        boot_idx = rng.choice(n_ref, n_samples_per_batch, replace=True)
        clr_bg = _clr(ref_rel[boot_idx])        # (n_samples, n_taxa)

        # --- inject phenotype signal in CLR space ---
        # Add a mean-zero shared factor z to the n_signal_taxa most
        # prevalent taxa for positive samples.  Centering z keeps CLR
        # means invariant; residual variance increase ≈ increment² ≪ 1.
        z = rng.normal(0, 1, n_pos)
        z -= z.mean()
        clr_bg[pos_idx, :n_signal_taxa] += (
            signal_corr_increment * z[:, np.newaxis])

        true_rel = _inv_clr(clr_bg)             # (n_samples, n_taxa)

        # --- multiplicative batch bias (per-taxon, shared within batch) ---
        log_bias = rng.normal(0, bias_strength * 0.5, size=n_taxa)
        bias = np.exp(log_bias)
        true_bias_list.append(bias)
        biased = true_rel * bias[np.newaxis, :]
        biased /= biased.sum(axis=1, keepdims=True)

        # --- per-sample log-space noise ---
        noise = rng.normal(0, sample_noise_std, (n_samples_per_batch, n_taxa))
        biased = biased * np.exp(noise)
        biased /= biased.sum(axis=1, keepdims=True)

        # --- batch-specific detection-limit truncation ---
        batch_rng = np.random.RandomState(random_seed + batch_idx * 7919)
        hard_taxa = batch_rng.choice(
            n_taxa, size=int(n_taxa * 0.4), replace=False)
        for t in hard_taxa:
            cutoff_pct = batch_rng.uniform(20, 50)
            threshold = np.percentile(biased[:, t], cutoff_pct)
            biased[biased[:, t] < threshold, t] = 0.0
        biased_sums = biased.sum(axis=1, keepdims=True)
        biased_sums[biased_sums == 0] = 1.0
        biased /= biased_sums

        # --- Dirichlet-Multinomial count sampling ---
        depth = int(sequencing_depth * batch_rng.uniform(0.5, 1.5))
        dm_alpha_total = dm_concentration * n_taxa
        counts = np.zeros((n_samples_per_batch, n_taxa), dtype=np.int64)
        for s in range(n_samples_per_batch):
            alpha = biased[s] * dm_alpha_total
            sampled_pi = batch_rng.dirichlet(np.maximum(alpha, 1e-6))
            counts[s] = batch_rng.multinomial(depth, sampled_pi)

        batch_data_list.append(counts.astype(np.float64))
        label_list.append(labels)

    return batch_data_list, label_list, true_bias_list


# -------------------------------------------------------------------------
# Evaluation utility
# -------------------------------------------------------------------------

def evaluate_cross_batch(features, labels, batch_labels, method_name="",
                          auto_reduce=True):
    """
    Leave-one-batch-out cross-batch prediction evaluation.

    Parameters
    ----------
    features : np.ndarray, shape (n_samples, n_features)
    labels : np.ndarray, shape (n_samples,)
    batch_labels : np.ndarray, shape (n_samples,)
    method_name : str
    auto_reduce : bool
        If True and n_features > 30, apply PCA to 30 components before
        classification (fitted on training folds only to prevent data leakage).

    Returns
    -------
    mean_auroc : float
    pooled_auroc : float
    """
    logo = LeaveOneGroupOut()
    fold_aurocs = []
    all_probs = np.zeros_like(labels, dtype=np.float64)

    for train_idx, test_idx in logo.split(features, labels, batch_labels):
        steps = [("scaler", StandardScaler())]
        if auto_reduce and features.shape[1] > 30:
            n_components = min(30, len(train_idx) - 1, features.shape[1])
            steps.append(("pca", PCA(n_components=n_components)))
        steps.append(("clf", LogisticRegression(max_iter=2000, C=1.0)))
        pipe = Pipeline(steps)

        pipe.fit(features[train_idx], labels[train_idx])
        probs = pipe.predict_proba(features[test_idx])[:, 1]
        all_probs[test_idx] = probs

        if len(np.unique(labels[test_idx])) > 1:
            auc = roc_auc_score(labels[test_idx], probs)
            fold_aurocs.append(auc)
            held_out = batch_labels[test_idx[0]]
            print(f"    Held-out batch {held_out}  auROC = {auc:.4f}")

    mean_auroc = float(np.mean(fold_aurocs)) if fold_aurocs else float("nan")
    pooled_auroc = (roc_auc_score(labels, all_probs)
                    if len(np.unique(labels)) > 1 else float("nan"))
    print(f"  [{method_name}]  mean auROC = {mean_auroc:.4f}  "
          f"pooled auROC = {pooled_auroc:.4f}")
    return mean_auroc, pooled_auroc


# -------------------------------------------------------------------------
# Main demo
# -------------------------------------------------------------------------

def run_demo():
    print("\n" + "#" * 60)
    print("#  MicroCoRE Framework -- End-to-End Comparison Demo")
    print("#" * 60)
    print("""
Scenario: semisynthetic 'pure co-variation signal'
  - Background: bootstrapped from Wirbel et al. 2019 CRC meta-analysis
    (849 samples, 8 studies); inherits real sparsity and correlation structure.
  - Signal: mean-zero CLR-space shared factor on signal taxa for positive
    samples (increases pairwise correlation, marginal means unchanged).
  - Batch effects: multiplicative bias + per-sample noise + detection truncation.
  - Count model: Dirichlet-Multinomial (overdispersed).

Methods compared
  1. CLR baseline         (abundance-level, expected ~0.50)
  2. MicroCoRE/binary     (phi-coefficient correlation -> GAE)
  3. MicroCoRE/abundance  (CLR-Pearson correlation -> GAE)
  4. MicroCoRE/fusion     (weighted phi + CLR-Pearson -> GAE)
""")

    # ------------------------------------------------------------------
    # Step 1: load reference data + generate semisynthetic batches
    # ------------------------------------------------------------------
    print("[ Step 1 ] Loading reference data and generating semisynthetic batches")
    print("  Reference: Wirbel et al. 2019 CRC meta-analysis (zellerlab/crc_meta)")
    ref_rel, ref_study_labels = load_reference_data(n_taxa=120)

    batch_data_list, label_list, true_bias_list = generate_synthetic_data(
        n_batches=4, n_samples_per_batch=120,
        n_signal_taxa=15, signal_corr_increment=0.22,
        bias_strength=2.0, dm_concentration=0.1, sample_noise_std=0.3,
        sequencing_depth=10000, reference_data=ref_rel, random_seed=2024)

    n_total = sum(d.shape[0] for d in batch_data_list)
    n_taxa  = batch_data_list[0].shape[1]
    print(f"  {len(batch_data_list)} batches, {n_total} total samples, "
          f"{n_taxa} taxa")

    bias_mat = np.array(true_bias_list)
    print(f"  True bias factor cross-batch std: "
          f"{bias_mat.std(axis=0).mean():.3f} ± {bias_mat.std(axis=0).std():.3f}")

    all_labels = np.concatenate(label_list)
    all_batches = np.concatenate([
        np.full(d.shape[0], i) for i, d in enumerate(batch_data_list)])

    # ------------------------------------------------------------------
    # Step 2: CLR baseline
    # ------------------------------------------------------------------
    print("\n[ Step 2 ] Baseline: CLR transform + logistic regression")
    all_counts = np.vstack(batch_data_list) + 0.5
    rel = all_counts / all_counts.sum(axis=1, keepdims=True)
    clr_features = np.log(rel) - np.log(rel).mean(axis=1, keepdims=True)

    clr_mean, clr_pooled = evaluate_cross_batch(
        clr_features, all_labels, all_batches, method_name="CLR baseline")

    # ------------------------------------------------------------------
    # Step 3: all three MicroCoRE variants in one model
    # ------------------------------------------------------------------
    print("\n[ Step 3 ] MicroCoRE pipeline -- all three correlation modes")
    print(f"  Device: {DEVICE}")
    model = MicroCoREModel(
        embed_dim=16,
        corr_mode="all",          # trains binary, abundance, and fusion
        binary_threshold=5,
        fusion_weights=(0.4, 0.6),
        cv_threshold=0.5,
        min_abs_corr=0.05,
        sparsity_threshold=0.1,
        n_epochs=1500,
        lr=0.01,
        smoothness_weight=0.0,
        weight_mode="cooccurrence",
        device=DEVICE,
        random_seed=42)

    embeddings_dict, batch_labels = model.fit_transform(batch_data_list)

    # ------------------------------------------------------------------
    # Step 4: evaluate each MicroCoRE variant
    # ------------------------------------------------------------------
    print("\n[ Step 4 ] Cross-batch prediction -- per-mode evaluation")
    results = {}
    for mode in ["binary", "abundance", "fusion"]:
        print(f"\n  -- MicroCoRE/{mode} --")
        mean_auc, pooled_auc = evaluate_cross_batch(
            embeddings_dict[mode], all_labels, batch_labels,
            method_name=f"MicroCoRE/{mode}")
        results[mode] = (mean_auc, pooled_auc)

    # ------------------------------------------------------------------
    # Step 5: summary table
    # ------------------------------------------------------------------
    print("\n" + "=" * 66)
    print("           Cross-batch prediction performance summary")
    print("=" * 66)
    print(f"  {'Method':<28}  {'mean auROC':>10}  {'pooled auROC':>12}  {'vs CLR':>8}")
    print(f"  {'-'*28}  {'-'*10}  {'-'*12}  {'-'*8}")
    print(f"  {'CLR baseline':<28}  {clr_mean:>10.4f}  {clr_pooled:>12.4f}  {'--':>8}")

    mode_labels = {
        "binary":    "MicroCoRE/binary  (phi)",
        "abundance": "MicroCoRE/abundance (CLR-P)",
        "fusion":    "MicroCoRE/fusion",
    }
    best_mode = max(results, key=lambda m: results[m][0])
    for mode in ["binary", "abundance", "fusion"]:
        m_auc, p_auc = results[mode]
        delta = m_auc - clr_mean
        marker = "  <-- best" if mode == best_mode else ""
        print(f"  {mode_labels[mode]:<28}  {m_auc:>10.4f}  {p_auc:>12.4f}  "
              f"{delta:>+8.4f}{marker}")
    print("=" * 66)

    # ------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------
    best_mean = results[best_mode][0]
    delta_best = best_mean - clr_mean

    print()
    if delta_best > 0.02:
        print(f"  MicroCoRE/{best_mode} gives the largest gain ({delta_best:+.4f}) "
              f"over the CLR baseline.")

        # Compare fusion vs others
        fusion_mean = results["fusion"][0]
        binary_mean = results["binary"][0]
        abund_mean  = results["abundance"][0]

        if best_mode == "fusion":
            print("  Fusion outperforms individual modes, consistent with the")
            print("  hypothesis that combining binary co-occurrence and continuous")
            print("  co-variation captures complementary correlation structure.")
        elif best_mode == "binary":
            print("  phi (binary) mode outperforms others, suggesting that")
            print("  presence/absence information dominates the phenotype signal")
            print("  in this dataset.")
        else:
            print("  CLR-Pearson mode outperforms others, suggesting that")
            print("  continuous abundance co-variation carries more signal than")
            print("  binary co-occurrence in this dataset.")

        if abs(fusion_mean - max(binary_mean, abund_mean)) < 0.01:
            print("  Note: fusion and the best individual mode are within 0.01 --")
            print("  the gain from fusion is marginal here.")
    else:
        print("  No method shows a clear advantage (delta < 0.02).")
        print("  Consider increasing bias_strength or n_samples_per_batch.")


if __name__ == "__main__":
    run_demo()
