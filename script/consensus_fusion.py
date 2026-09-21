"""
MicroCoRE Framework -- Stage 2
Cross-batch consensus correlation fusion

Core idea:
    Merge per-batch correlation matrices into a single consensus matrix.
    Sample-size-weighted aggregation combined with cross-batch stability
    filtering retains only edges that are reproducible across batches,
    removing batch-specific spurious correlations.
"""

import numpy as np


def weighted_aggregate(corr_matrix_list, batch_sample_counts):
    """
    Compute a sample-size-weighted average of per-batch correlation matrices.
    Batches with more samples receive higher weight (more reliable estimates).

    Parameters
    ----------
    corr_matrix_list : list[np.ndarray]
    batch_sample_counts : list[int]

    Returns
    -------
    weighted_mean : np.ndarray
    """
    weights = np.array(batch_sample_counts, dtype=np.float64)
    weights = weights / weights.sum()

    stacked = np.stack(corr_matrix_list, axis=0)   # (n_batches, n_taxa, n_taxa)
    weighted_mean = np.tensordot(weights, stacked, axes=(0, 0))
    return weighted_mean


def cross_batch_stability_filter(corr_matrix_list,
                                 cv_threshold=0.5,
                                 min_abs_corr=0.1):
    """
    Compute the cross-batch coefficient of variation (CV) for every edge
    and mask unstable edges.

    For edge (k1, k2) the CV = std / |mean| across batches measures
    reproducibility. High-CV edges are likely noise or batch-specific
    artefacts and should be discarded.

    Parameters
    ----------
    corr_matrix_list : list[np.ndarray]
    cv_threshold : float
        Edges whose CV exceeds this value are zeroed out.
    min_abs_corr : float
        Edges with |mean| below this value are treated as negligible and
        zeroed out (avoids inflated CV due to near-zero denominators).

    Returns
    -------
    stable_mask : np.ndarray, dtype bool, shape (n_taxa, n_taxa)
        True means the edge is stable and should be kept.
    """
    stacked = np.stack(corr_matrix_list, axis=0)
    edge_mean = stacked.mean(axis=0)
    edge_std = stacked.std(axis=0)

    # Build CV matrix; set very weak edges to inf so they are removed
    denom = np.abs(edge_mean)
    weak_mask = denom < min_abs_corr
    denom[weak_mask] = 1.0                  # placeholder to avoid divide-by-zero
    cv_matrix = edge_std / denom
    cv_matrix[weak_mask] = np.inf           # weak edges are unconditionally removed

    stable_mask = cv_matrix <= cv_threshold
    np.fill_diagonal(stable_mask, True)     # always keep the diagonal
    return stable_mask


def build_consensus_matrix(corr_matrix_list, batch_sample_counts,
                            cv_threshold=0.5, min_abs_corr=0.1,
                            sparsity_threshold=0.0):
    """
    Full consensus matrix construction pipeline.

    Parameters
    ----------
    corr_matrix_list : list[np.ndarray]
    batch_sample_counts : list[int]
    cv_threshold : float
        Cross-batch stability filter threshold.
    min_abs_corr : float
        Weak-edge filter threshold.
    sparsity_threshold : float
        Final entries with |value| below this are zeroed (extra sparsification).

    Returns
    -------
    consensus_matrix : np.ndarray
        Sample-size-weighted, stability-filtered, and sparsified correlation
        network.
    """
    # Step 1: weighted aggregation
    weighted_mean = weighted_aggregate(corr_matrix_list, batch_sample_counts)

    # Step 2: stability filtering
    stable_mask = cross_batch_stability_filter(corr_matrix_list,
                                               cv_threshold=cv_threshold,
                                               min_abs_corr=min_abs_corr)
    consensus = weighted_mean * stable_mask.astype(np.float64)

    # Step 3: optional extra sparsification
    if sparsity_threshold > 0:
        consensus = consensus * (np.abs(consensus) >= sparsity_threshold)

    # Enforce symmetry
    consensus = (consensus + consensus.T) / 2
    np.fill_diagonal(consensus, 1.0)

    n_kept = (np.abs(consensus) > 0).sum() // 2
    n_total = consensus.shape[0] * (consensus.shape[0] - 1) // 2
    print(f"  Consensus matrix built -- edges kept: {n_kept} / {n_total} "
          f"({100 * n_kept / n_total:.1f}%)")
    return consensus
