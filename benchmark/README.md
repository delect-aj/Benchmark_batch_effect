# Microbiome batch-correction benchmark

Full design: [PLAN.md](PLAN.md).

## Method wrapper contract

Every method is one script in `methods/`, called the same way:

```
Rscript methods/<name>.R  counts.tsv meta.tsv out.tsv
python  methods/<name>.py counts.tsv meta.tsv out.tsv
```

| File | Format |
|---|---|
| `counts.tsv` | rows = samples, first column `sample_id`, remaining columns = taxa, integer read counts |
| `meta.tsv` | first column `sample_id` (same order as counts), required `batch`, `phenotype` (binary; `0` = control) |
| `out.tsv` | rows = same samples/order, first column `sample_id`; columns = taxa (tracks A) or `dim1..k` (track B) |
| `out.tsv.json` | `{"kind": ...}` — one of `counts, relabund, clr, log, percentile, embedding`; tells metrics which scale the output is on |

Rules applied to every wrapper (neutrality):
- Use the authors' default parameters; any deviation is commented in the wrapper.
- If a method accepts a protected biological covariate, pass `phenotype` (`bio_design()` in `_common.R`). If this makes it fail (e.g. fully confounded design), let it fail. A failure is a result.
- A wrapper exits non-zero on failure. Runtime/memory come from the workflow (Snakemake `benchmark:`), not from the wrapper.
- Wrappers are **transductive**: they are fitted on the whole input. For the leave-one-study-out prediction track, fit on the training studies only, then transform the held-out study (added later).

Optional side inputs (`taxonomy.tsv`, `tree.nwk`) sit next to `counts.tsv`; wrappers read them with `side_file()`.

Shared helpers: `methods/_common.R` (`read_input`, `write_output`, `side_file`, `clr`, `bio_design`) and `methods/_common.py`.

## Wrapper status

| Wrapper | Track | Output kind | Notes |
|---|---|---|---|
| raw | baseline | counts | |
| limma | A | clr | `removeBatchEffect` on CLR |
| combat | A | clr | `sva::ComBat` on CLR |
| combatseq | A | counts | `sva::ComBat_seq` |
| mmuphin | A | relabund | `MMUPHin::adjust_batch` |
| conqur | A | counts | reference batch = first level (author default is user-chosen) |
| plsdabatch / wplsdabatch | A | clr | `balance = TRUE / FALSE` |
| percentile | A | percentile | Gibbons 2018; controls = `phenotype == 0` |
| harmony | B | embedding | on the top 20 CLR principal components |
| fastmnn | B | embedding | `batchelor::fastMNN`, d = 20 |
| debiasm | A | relabund | Python, `DebiasMClassifier.transform` on raw counts |
| metadict | A | counts | uses `tree.nwk` or `taxonomy.tsv` next to counts.tsv; with neither, falls back to flat taxon distances |
| cqr | A | counts | **partial reimplementation** of Park 2025 (released code does not run): robust-CV reference batch + ConQuR composite quantile regression, NB step omitted |
| ruviiinb | A | counts | ZINB; replicate sets = `meta$replicate` or phenotype groups; controls = 20% of taxa least associated with phenotype |
| scanvi | B | embedding | Python, scvi-tools SCVI → SCANVI with defaults |

Not written yet: metacal (only runs on mock-community data). Track C methods (MaAsLin2, ANCOM-BC2, BDMMA, SVA) output differential-abundance results rather than tables, so they will get their own contract. ANCOM-BC2 is not installed yet: version 2.8 needs CVXR < 1.0, which conda-forge does not provide; install it from the CRAN archive when Track C starts.

## Pilot

```
bash pilot/run_pilot.sh            # local: uses Rscript / python3 on PATH
RSCRIPT=~/software/miniconda3/envs/bench-r/bin/Rscript PYTHON=~/software/miniconda3/envs/bench-py/bin/python bash pilot/run_pilot.sh
```
This generates a toy dataset (`simulate/toy_sim.R`), runs every wrapper, and records `ok` or `FAIL` for each. `pilot/check.R` then validates the outputs and prints batch/phenotype PERMANOVA R². Outputs go to `results/` (gitignored).

Environments: `bash envs/build.sh <dir with GitHub clones>` builds `bench-r` and `bench-py` (one env per language). The GitHub-only packages are cloned on the login node, because compute nodes cannot reach GitHub.
