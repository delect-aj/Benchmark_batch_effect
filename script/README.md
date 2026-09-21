# MicroCoRE — Microbiome Co-occurrence Representation Framework

MicroCoRE addresses microbiome batch effects by exploiting a key mathematical
invariant: **within-batch taxon-taxon correlations are unaffected by
multiplicative batch biases**.  Instead of estimating and removing batch
effects (as DEBIAS-M and MetaDICT do), MicroCoRE builds sample representations
entirely from bias-invariant co-variation structure.

## File Structure

| File | Stage | Responsibility |
|------|-------|----------------|
| `correlation_builder.py` | Stage 1 | Within-batch correlation matrices (phi / CLR-Pearson / fusion) |
| `consensus_fusion.py`    | Stage 2 | Cross-batch stability filtering + weighted aggregation |
| `representation_learner.py` | Stage 3 | Graph autoencoder: taxon embeddings from consensus network |
| `sample_aggregator.py`   | Stage 4 | Aggregate taxon embeddings into sample vectors |
| `pipeline.py`            | End-to-end | `MicroCoREModel` class with `.fit()` / `.transform()` API |
| `demo.py`                | Demo | Synthetic data + leave-one-batch-out evaluation |

## Quick Start

```python
from pipeline import MicroCoREModel

# batch_data_list: list of np.ndarray, each shaped (n_samples, n_taxa)
model = MicroCoREModel(
    embed_dim=32,
    corr_mode="fusion",          # 'binary' | 'abundance' | 'fusion'
    cv_threshold=0.5,            # cross-batch edge stability threshold
    weight_mode="cooccurrence",  # 'rank' | 'clr' | 'binary' | 'cooccurrence'
    n_epochs=2000,
)

sample_embeddings, batch_labels = model.fit_transform(batch_data_list)
taxon_embeddings = model.taxon_embeddings   # shape (n_taxa, embed_dim)
consensus_network = model.consensus_matrix  # stable cross-batch co-variation network
```

## Why Within-batch Correlations are Bias-Invariant

Under the multiplicative bias model  O_{i,j,k} = w_{i,k} · A_{i,j,k} · c_{i,j},
the bias factor w_{i,k} is constant within batch i.  Constant factors:

- cancel in Pearson correlation (absorbed by standard deviation normalisation),
- do not change binary presence/absence judgments (unless bias is extreme),
- do not change within-batch rank orderings.

Therefore, within-batch correlation structure ≈ true correlation structure,
without any bias estimation.

## Two Sample Representation Strategies

**First-order aggregation** (`weight_mode='rank'`):

    s_i = sum_k  rank(O_{i,k}) * e_k

Suitable when the phenotype signal is in abundance levels.
`rank` is invariant to multiplicative biases within a batch.

**Second-order co-occurrence** (`weight_mode='cooccurrence'`):

    S_i = E^T @ diag(b_i) @ E   (flattened)

where b_i is the binary presence vector for sample i.  Captures which
embedding directions co-occur in sample i, reflecting sample-specific
co-variation patterns.  **Use this mode when the signal is in co-variation
structure rather than abundance levels.**

## Demo Results

In the 'pure co-variation signal' scenario (no abundance-level signal):

| Method | Mean leave-one-batch auROC | Interpretation |
|--------|--------------------------|----------------|
| CLR baseline | ~0.50 | Fails — no abundance-level signal |
| MicroCoRE | ~0.55 | Recovers signal from co-variation structure |

## Design Parameters

| Parameter | Effect | Tuning priority |
|-----------|--------|-----------------|
| `cv_threshold` | Controls consensus network density | High |
| `embed_dim` | Taxon embedding dimensionality | High |
| `weight_mode` | Aggregation strategy | High |
| `sparsity_threshold` | Additional edge pruning | Medium |
| `fusion_weights` | phi vs CLR-Pearson blend | Low |

## Comparison with Related Methods

| Dimension | DEBIAS-M | MetaDICT | MicroCoRE |
|-----------|----------|----------|-----------|
| Strategy | Estimate & remove biases | Dictionary decomposition | Exploit bias-invariant structure |
| Requires bias estimation | Yes | Yes | No |
| Main output | Corrected abundances + predictor | Dual embeddings + corrected data | Taxon embeddings + sample vectors |
| Robustness to unobserved confounders | Moderate | High | High |
| Main risk | Fails with small batches | Non-convex local optima | Loses abundance-level signal |

## Practical Recommendations

1. **Signal in both levels and structure** (typical real data): concatenate
   MicroCoRE embeddings with CLR features and use an ensemble.

2. **Minimum sample size**: at least 50 samples per batch for reliable
   within-batch correlation estimates.

3. **Minimum number of batches**: at least 3 to benefit from cross-batch
   stability filtering.
