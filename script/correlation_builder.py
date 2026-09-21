"""
MicroCoRE Framework -- Stage 1
Within-batch correlation matrix construction

Core idea:
    Compute taxon-taxon correlations independently within each batch.
    Because multiplicative batch biases are constant within a batch,
    they cancel out in Pearson correlation, binary judgments, and rank
    orderings. Within-batch correlation structure is therefore a natural
    invariant of multiplicative batch effects.
"""

import numpy as np


def compute_phi_coefficient(abundance_matrix, detection_threshold=5):
    """
    Compute the phi coefficient (binary Pearson correlation) from
    presence/absence data.

    Parameters
    ----------
    abundance_matrix : np.ndarray, shape (n_samples, n_taxa)
        Raw count matrix for one batch.
    detection_threshold : int
        Minimum read count to call a taxon 'present'.

    Returns
    -------
    phi_matrix : np.ndarray, shape (n_taxa, n_taxa)
        Symmetric matrix with ones on the diagonal.
    """
    binary_matrix = (abundance_matrix >= detection_threshold).astype(np.float64)
    n_samples, n_taxa = binary_matrix.shape

    # Column-wise mean and standard deviation (vectorised)
    col_mean = binary_matrix.mean(axis=0)
    centered = binary_matrix - col_mean
    std = np.sqrt((centered ** 2).mean(axis=0))
    # Prevent division by zero for invariant taxa
    std[std < 1e-10] = 1e-10

    cov_matrix = (centered.T @ centered) / n_samples
    phi_matrix = cov_matrix / np.outer(std, std)

    # Numerical clean-up
    phi_matrix = np.nan_to_num(phi_matrix, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(phi_matrix, 1.0)
    phi_matrix = np.clip(phi_matrix, -1.0, 1.0)
    return phi_matrix


def compute_clr_pearson(abundance_matrix, pseudocount=0.5):
    """
    Simplified SparCC-style correlation: apply centered log-ratio (CLR)
    transformation to relative abundances, then compute Pearson correlation.

    Full SparCC iteratively estimates true correlations under a sparse
    background; CLR-Pearson is a cheaper approximation that inherits the
    core advantage of resistance to compositional bias.

    Parameters
    ----------
    abundance_matrix : np.ndarray, shape (n_samples, n_taxa)
    pseudocount : float
        Added to all counts to avoid log(0).

    Returns
    -------
    corr_matrix : np.ndarray, shape (n_taxa, n_taxa)
    """
    # Add pseudocount and row-normalise to relative abundance
    counts = abundance_matrix.astype(np.float64) + pseudocount
    row_sums = counts.sum(axis=1, keepdims=True)
    rel_abund = counts / row_sums

    # CLR: subtract row-wise log mean
    log_rel = np.log(rel_abund)
    clr = log_rel - log_rel.mean(axis=1, keepdims=True)

    # Pearson correlation
    n_samples = clr.shape[0]
    col_mean = clr.mean(axis=0)
    demeaned = clr - col_mean
    std = np.sqrt((demeaned ** 2).mean(axis=0))
    std[std < 1e-10] = 1e-10

    corr_matrix = (demeaned.T @ demeaned) / n_samples / np.outer(std, std)
    corr_matrix = np.nan_to_num(corr_matrix, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr_matrix, 1.0)
    corr_matrix = np.clip(corr_matrix, -1.0, 1.0)
    return corr_matrix


def build_batch_correlation_matrices(batch_data_list, mode="fusion",
                                     binary_threshold=5,
                                     fusion_weights=(0.4, 0.6)):
    """
    Compute a within-batch correlation matrix for every batch independently.

    Parameters
    ----------
    batch_data_list : list[np.ndarray]
        One array per batch, each of shape (n_samples_in_batch, n_taxa).
    mode : str
        'binary'  -- phi coefficient only
        'abundance' -- CLR-Pearson only
        'fusion'  -- weighted combination (recommended)
    binary_threshold : int
        Detection threshold used for the phi coefficient.
    fusion_weights : tuple
        (phi weight, CLR-Pearson weight); only used in 'fusion' mode.

    Returns
    -------
    corr_matrix_list : list[np.ndarray]
        One correlation matrix per batch.
    """
    if mode not in {"binary", "abundance", "fusion"}:
        raise ValueError(f"Unknown mode: {mode}")

    corr_matrix_list = []
    for batch_idx, batch_data in enumerate(batch_data_list):
        if mode == "binary":
            mat = compute_phi_coefficient(batch_data,
                                          detection_threshold=binary_threshold)
        elif mode == "abundance":
            mat = compute_clr_pearson(batch_data)
        else:  # fusion
            phi_mat = compute_phi_coefficient(batch_data,
                                              detection_threshold=binary_threshold)
            clr_mat = compute_clr_pearson(batch_data)
            mat = fusion_weights[0] * phi_mat + fusion_weights[1] * clr_mat

        corr_matrix_list.append(mat)
        n_nonzero = (np.abs(mat) > 0.1).sum() // 2
        print(f"  Batch {batch_idx + 1}/{len(batch_data_list)} -- "
              f"shape {mat.shape}, non-zero edges: {n_nonzero}")

    return corr_matrix_list


def build_all_batch_correlation_matrices(batch_data_list,
                                         binary_threshold=5,
                                         fusion_weights=(0.4, 0.6)):
    """
    Compute phi, CLR-Pearson, and fusion correlation matrices for every batch
    in a single pass.  phi and CLR-Pearson are computed once per batch and
    reused to construct the fusion matrix, avoiding redundant computation.

    Parameters
    ----------
    batch_data_list : list[np.ndarray]
        One array per batch, each of shape (n_samples_in_batch, n_taxa).
    binary_threshold : int
        Detection threshold used for the phi coefficient.
    fusion_weights : tuple
        (phi weight, CLR-Pearson weight) for the fusion matrix.

    Returns
    -------
    dict with keys 'binary', 'abundance', 'fusion', each mapping to a
    list[np.ndarray] of per-batch correlation matrices.
    """
    binary_list, abundance_list, fusion_list = [], [], []

    for batch_idx, batch_data in enumerate(batch_data_list):
        phi_mat = compute_phi_coefficient(batch_data,
                                          detection_threshold=binary_threshold)
        clr_mat = compute_clr_pearson(batch_data)
        fusion_mat = fusion_weights[0] * phi_mat + fusion_weights[1] * clr_mat

        binary_list.append(phi_mat)
        abundance_list.append(clr_mat)
        fusion_list.append(fusion_mat)

        n_b = (np.abs(phi_mat) > 0.1).sum() // 2
        n_a = (np.abs(clr_mat) > 0.1).sum() // 2
        n_f = (np.abs(fusion_mat) > 0.1).sum() // 2
        print(f"  Batch {batch_idx + 1}/{len(batch_data_list)} -- "
              f"shape {phi_mat.shape} | "
              f"binary edges: {n_b}, abundance edges: {n_a}, fusion edges: {n_f}")

    return {"binary": binary_list, "abundance": abundance_list, "fusion": fusion_list}
