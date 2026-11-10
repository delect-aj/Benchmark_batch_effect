# Introduction and Methods

> Draft manuscript sections. The dataset provenance — how each of the thirteen datasets was obtained,
> reconstructed and preprocessed — is written up separately in [`methods_data.md`](methods_data.md) and
> belongs under *Methods → Datasets*; the summary here is deliberately short so the two do not duplicate
> each other. References are numbered across both sections and listed at the end.

---

## Introduction

Microbiome sequencing measurements carry a technical signature of the laboratory that produced them.
DNA extraction chemistry, primer choice, library preparation, sequencing platform and bioinformatic
pipeline each leave a reproducible imprint on the observed community, and controlled multi-laboratory
experiments have shown that this imprint can rival or exceed the biological contrast a study was designed
to detect [1,2]. The consequence is familiar from other high-throughput assays [3]: when samples are
pooled across studies, the dominant axis of variation is often the study itself, and any downstream
inference that ignores this is confounded by construction.

The pressure to pool is nonetheless considerable. Individual microbiome studies are small, their effect
sizes modest, and their findings frequently fail to replicate; meta-analysis across cohorts has become the
standard route to signatures that generalise [4–6]. This has produced a steady supply of batch-correction
methods. Some are imported from bulk transcriptomics, where the location-scale model of ComBat and the
linear-model adjustment of limma were developed [7–9]. Others come from single-cell genomics, where
integration methods such as Harmony, mutual-nearest-neighbour correction and scANVI return a shared
embedding rather than a corrected feature table [10–12]. A third group was built for microbiome data
specifically, and addresses properties the imported methods do not model: compositionality, excess zeros,
overdispersion and a sequencing depth that is an artefact of the run rather than a property of the
community [13]. MMUPHin, ConQuR, PLSDA-batch, percentile normalisation, DEBIAS-M and MetaDICT belong to
this group [14–19].

Choosing among them is not straightforward, and three difficulties recur.

**There is no gold standard, and the proxy criteria need not agree.** Correction is almost always
evaluated by proxy: the variance attributable to batch should fall, the variance attributable to the
biological variable should not, and batches should mix in an ordination. None of these is the quantity of
interest. They are also not one criterion but several, and they are computed in a space the analyst
chooses — Bray–Curtis, Aitchison or Jaccard distance, a silhouette width, a local-neighbourhood mixing
score — with no guarantee that the resulting rankings coincide. A method can be ranked first under one
distance and in the lower half under another, and the published comparisons rarely report whether their
conclusion survives that choice.

**The literature disagrees about whether correction helps downstream.** Removing the batch term reliably
improves the appearance of a pooled ordination, but the evidence that it improves inference is mixed.
Cross-study prediction studies have reported that aggressive correction can remove generalisable signal
along with the batch effect, while others report the opposite [5,6,16]. Differential-abundance results are
more unstable still: methods disagree with each other across datasets even before batch correction is
applied [20]. Without a ground truth, these disagreements cannot be adjudicated, only restated.

**Batch–phenotype confounding is ubiquitous but is rarely treated as an experimental factor.** In practice
the batch is usually the study, and studies differ in the populations they recruit; batch and phenotype are
therefore partially confounded in almost every real meta-analysis. Most correction methods accept the
biological variable as a protected covariate, which is exactly the right thing to do when the design is
balanced and increasingly dangerous as it is not: conditioning on a phenotype that batch nearly determines
can move batch structure onto the phenotype axis rather than removing it. Existing benchmarks typically
evaluate at whatever confounding their test datasets happen to have, so the breaking point of each method
is unknown.

These three problems are not independent. The first determines how any result should be read, the second
cannot be settled without the ground truth the first lacks, and the third requires criteria established by
the other two. We therefore designed a benchmark that addresses them in that order, and asks four
questions: **(Q1)** do the evaluation criteria agree, and does the ranking depend on which is chosen;
**(Q2)** with a known truth in hand, what did the methods actually remove and retain, and are the proxy
criteria faithful to it; **(Q3)** how does batch–phenotype confounding erode correction, and where does each
method break; and **(Q4)** does removing the batch term translate into better differential-abundance and
cross-study prediction.

Answering these requires more than one kind of data, because no single kind can answer them all. We
assembled three. Ten real case-control and multi-cohort datasets (7,115 samples) show whether a biological
contrast survives correction, but carry no ground truth. Designed protocol experiments, in which defined
mock communities or repeated aliquots of one specimen are processed under many laboratory protocols,
isolate a named technical factor and supply the true composition. Simulations provide complete ground truth
over a controlled grid, including a confounding axis that real data cannot be made to vary on demand. Every
method is driven through a single wrapper interface, so that the comparison differs only in the method, and
every method is run at its authors' defaults; where a method fails, the failure is recorded as a result
rather than excluded. The work is deliberately neutral: we propose no method of our own here, and no
method in the comparison is ours.

---

## Methods

### Benchmark design

The benchmark is organised around one interface and four families of criteria. Every method is wrapped as a
single script taking a sample-by-taxon count matrix and a sample table carrying a `batch` and a `phenotype`
column, and returning either a corrected table on the same taxa or a low-dimensional embedding on the same
samples. A sidecar declares which scale the output is on (`counts`, `relabund`, `clr`, `log`, `percentile`
or `embedding`), so that each criterion is computed in a space appropriate to the output rather than by
forcing every method onto one scale. Methods returning a table are referred to as **track A**, methods
returning an embedding as **track B**; criteria that have no meaning for an embedding, such as per-taxon
fidelity or differential abundance, are recorded as missing for track B rather than as zero.

Two baselines bound the comparison. `raw` is the uncorrected input and is the reference against which every
change is measured. `oracle` is defined only for simulated data, where it is the batch-free table the same
biology would have produced with no batch effect at all; it is the upper bound of what correction could
achieve, not a competitor.

### Datasets

Thirteen datasets comprising 7,687 biological samples were assembled, of which ten (7,115 samples) enter the
method comparison and three are reserved for calibrating the simulation. Seven are 16S rRNA amplicon
surveys and six are shotgun metagenomes. They fall into three regimes: curated multi-cohort case-control
collections with no ground truth (`crc_mgx`, `ibd_mgx`, `ici_mgx`, `crc_16s`, `ibd_16s`, `hiv_16s`,
`hivrc_16s`, `seed_16s`); designed protocol experiments in which the same material is processed under many
protocols (`mbqc`, `jmbc_extraction`, `jmbc_library`); and defined mock communities with a pipetted
composition (`brooks`, `salim_pig`, `jmbc_library`). Two datasets (`mbqc`, and the two JMBC experiments)
were reconstructed from raw sequencing reads because no usable processed feature table exists for them.

All datasets were brought to a common format: integer counts, amplicon data aggregated to genus and shotgun
data to species except for the defined mocks, which are reported at strain level; taxa present in fewer than
10% of samples removed; no rarefaction, transformation or compositional normalisation, those being choices
that belong to the methods under evaluation. The complete provenance, processing parameters and known
limitations of every dataset are given in *Datasets* (see [`methods_data.md`](methods_data.md)) and
Supplementary Table S1.

Two properties of this collection are used directly. `mbqc` is fully crossed — all 22 specimens processed
by all 13 laboratories — and therefore has Cramér's V = 0, providing the unconfounded anchor of the
confounding axis. `jmbc_library` consists of 28 library preparations of a single DNA mock and contains no
biological variation at all, so any structure a method leaves behind is batch and any difference it creates
is a false signal.

### Methods under evaluation

Fourteen correction methods were evaluated in fifteen configurations, alongside the two baselines and one
negative control (Table 1). Track A comprises limma's `removeBatchEffect` applied to centred log-ratio
coordinates [9], ComBat on the same coordinates and ComBat-seq on counts [7,8,21], MMUPHin's `adjust_batch`
[14], ConQuR in both its tuned and default configurations [15], a composite quantile-regression method
[22], PLSDA-batch in its balanced and weighted variants [16], percentile normalisation [17], MetaDICT [19]
and DEBIAS-M [18]. Track B comprises Harmony applied to the leading 20 principal components of the CLR
table [10], `fastMNN` with d = 20 [11,23] and scANVI [12,24].

Three rules were applied uniformly to keep the comparison neutral. First, every method runs at its authors'
default parameters; any deviation is documented in the wrapper. Second, where a method accepts a protected
biological covariate, the phenotype is passed; if this causes the method to fail — as it must in a fully
confounded design — the failure is recorded rather than worked around. Third, a wrapper that exits
non-zero is recorded as a failure for that dataset and contributes missing values, so that fragility is
measured rather than silently removed from the comparison. Runtime and peak memory are recorded by the
workflow, not self-reported by the wrapper.

Two configurations of ConQuR are included because they answer different questions. `conqur` runs the
authors' own tuning grid, which selects the quantile-regression specification and the reference batch; this
is the method as its authors intend it, and it accounts for more than half of the benchmark's total compute.
`conqur_default` runs the untuned method with the first batch level as reference, which is the setting most
users would apply in practice. A third variant, `conqur_permcov`, is a negative control rather than a
method: it is the untuned method given a randomly permuted phenotype as its protected covariate, and it
isolates how much of the observed differential-abundance inflation is attributable to conditioning on the
phenotype at all. It is shown only where that comparison is the point, and is excluded from rankings and
from the practicality comparison.

Two wrappers exist but do not contribute results. RUV-III-NB [25] could not be run: both documented entry
points fail inside their own dependencies under the current Bioconductor release, and the wrapper is kept
with the failure documented rather than deleted. The composite quantile-regression method [22] is a partial
reimplementation — the released code does not execute, and the negative-binomial step is omitted — and is
labelled as such wherever it appears.

### Evaluation criteria

Four families of criteria are computed, plus practicality.

**Batch removal.** The proportion of variance attributable to batch, estimated by PERMANOVA [26,27] under
three distances — Aitchison [28], Bray–Curtis [29] and Jaccard — together with the average silhouette width
of the batch labels [30] and the integration local inverse Simpson index (iLISI) [10]. The three distances
are reported separately rather than averaged, because whether they agree is itself one of the questions.
Each distance is computed only where it is defined for the output kind: compositions are mapped to Aitchison
coordinates, log- and CLR-scale outputs are mapped back to the simplex by softmax for Bray–Curtis, and
Jaccard is computed only for outputs that retain presence/absence structure, so that eight of the fifteen
configurations are necessarily unscored under it.

**Signal retention.** The same five quantities computed on the phenotype labels: PERMANOVA R² under the
three distances, the phenotype silhouette width, and cLISI, oriented so that larger is better.

**Fidelity to the truth.** For simulated data, the Mantel correlation between the sample-distance matrix of
the corrected table and that of the batch-free oracle table, and the mean per-taxon Spearman correlation
between corrected and oracle abundances. For the designed mock datasets, the Aitchison and Bray–Curtis
distance between each corrected sample and its pipetted composition; the per-taxon log-ratio of observed to
expected abundance, which the multiplicative bias model predicts should be constant within a batch [31]; the
residual variance attributable to batch among samples of identical material, which is false structure by
construction; and the number of taxa called different between batches of identical material, by
Kruskal–Wallis test.

**Downstream consequences.** Differential abundance is tested by Wilcoxon rank-sum test on the corrected
table with Benjamini–Hochberg control at 5% [32], and scored against two definitions of the truth. The first
is the set of taxa the simulation planted. The second is the set the identical test calls on the oracle
table — the answer the data would have given with no batch effect. The second is reported alongside the
first because the first is not a clean criterion: planting fold changes changes the community's total mass
by 1.5–2×, the simulator renormalises, and taxa that were never planted genuinely do shift, so counting
them as errors charges the method for a property of the simulation. False discovery rate, power and average
precision are reported against both. Cross-study prediction is evaluated by leave-one-study-out auROC with
L1-penalised logistic regression [33] and random forests [34,35], restricted to the seven datasets whose
batch is the study and whose phenotype is a disease contrast.

**Practicality.** Wall-clock time, peak resident memory and failure rate, recorded by the workflow for every
method × dataset combination.

The leave-one-study-out figures are an explicit **upper bound**, not an estimate. Correction is run once on
the whole dataset by the ordinary workflow, and this script only re-splits the already-corrected table, so
the correction saw the held-out study's labels. Two things make the bound informative nonetheless: `raw`
involves no correction and therefore no leakage, so a corrected table that does not beat `raw` even with the
leakage in its favour did not help; and the auROC is computed on exactly the table the other criteria are
computed on, making the comparison between "removes batch" and "predicts across studies" like for like.
Folds in which the penalised model selects no features are reported as degenerate rather than scored as
chance, which would otherwise erase the comparison.

### Aggregation and ranking

Within each dataset, each criterion is oriented so that larger is better and min–max scaled across methods,
following the convention established for single-cell integration benchmarks [36]; a family score is the mean
of its scaled criteria. Datasets on which every method ties, or on which a criterion does not apply,
contribute nothing rather than contributing a tie. Where the comparison concerns order rather than distance,
competition ranks are used instead. Batch–phenotype association is quantified throughout by Cramér's V
[37], normalised by `min(dim) − 1` so that it is bounded at 1 for phenotypes with more than two levels.

### Simulation

Simulated datasets were generated with MIDASim in parametric mode [38], fitted separately to a 16S and a
shotgun template drawn from the real collection, so that sparsity and taxon–taxon correlation reflect real
data of each assay type. Differentially abundant taxa were drawn from those with prevalence ≥ 20% and their
mean relative abundance shifted by a specified log2 fold change in cases, with direction randomised.

Batch effects are superimposed on the true relative abundances, so that the batch-free table is known
exactly, under five models applied per batch: per-taxon multiplicative bias with log-normal magnitude (the
McLaren model [31]); per-taxon variance scaling about each taxon's batch mean; batch-specific excess
dropout of present taxa; restriction of the bias to a random subset of taxa, giving non-systematic effects;
and a between-batch ratio of sequencing depth. Counts are drawn multinomially from the resulting relative
abundances, and the oracle table is drawn from the uncontaminated ones at the same library sizes. The same
prevalence filter applied to the real data is then applied to both.

Batch–phenotype confounding is a designed parameter rather than an incidental property: the phenotype is
drawn with a batch-dependent probability, with alternating batches leaning towards case and control, so
that the realised Cramér's V spans 0 to 1 on demand.

The grid comprises **213 scenarios** of three kinds. Ten *sweep* scenarios vary confounding alone across
V ∈ {0, 0.3, 0.6, 0.9, 1} on both templates with every other factor at its default, and carry the gradient
analysis. Three *null* scenarios contain no biological signal at all and quantify how much phenotype
structure each method invents. Two hundred *LHS* scenarios form a discrete Latin hypercube over twelve
factors — template, number of batches, samples per batch, batch size imbalance, confounding, proportion of
DA taxa, DA effect size, bias magnitude, fraction of taxa affected, variance shift, extra dropout and depth
ratio — in which every level of every factor appears equally often with levels permuted independently
across factors [39], so that interactions are covered without a full factorial. Sweep and null scenarios are
replicated 20 times and LHS scenarios 5 times, with fixed seeds throughout; the scenario table is frozen and
is not regenerated once results exist.

### Workflow and reproducibility

The benchmark is a Snakemake workflow [40] in which every method × dataset combination is one job and every
score is one file, so that a failure is localised and the whole is restartable. The complete run comprises
33,250 jobs over 1,270 dataset instances (10 real and 1,260 simulated), yielding 20,840 method × dataset
score records of which 19,465 completed and 1,375 failed, plus 1,409 leave-one-study-out folds and 52
known-composition evaluations. A method exceeding a six-hour limit is recorded as failed rather than
dropped. BLAS thread counts are pinned per rule so that a job requesting one core does not take the node.

Two conda environments pin the software: one for R and one for Python, built by a script in the repository.
Analyses used R 4.4 with vegan, cluster, glmnet, randomForest and pROC, and Python 3.10 with scvi-tools and
the DEBIAS-M package; MBQC denoising additionally required R 4.1.3 for DADA2. Random seeds are fixed for
simulation, for scenario design and for every stochastic figure element. All code, the assembled datasets,
the scenario table and the scoring scripts are available at
`https://github.com/delect-aj/Benchmark_batch_effect`.

### Statistical reporting

Variance components are PERMANOVA R² with no permutation testing, since the quantity rather than its
significance is of interest. Differential abundance uses two-sided Wilcoxon rank-sum tests with
Benjamini–Hochberg control at 5%. Gradient curves report the mean over replicates with one standard
deviation; box plots show the median, interquartile range and 1.5× IQR whiskers. Error rates pooled across
scenarios are weighted by the number of calls, because most scenarios call only a handful of taxa and a
per-scenario ratio with a denominator of one or two is 0 or 1 by arithmetic rather than by behaviour. No
multiplicity correction is applied across methods or criteria; the benchmark is descriptive and reports
effect sizes rather than testing hypotheses about method superiority.

---

## References

1. Sinha R, Abu-Ali G, Vogtmann E, *et al.* Assessment of variation in microbial community amplicon sequencing by the Microbiome Quality Control (MBQC) project consortium. *Nature Biotechnology* 2017;**35**:1077–1086.
2. Costea PI, Zeller G, Sunagawa S, *et al.* Towards standards for human fecal sample processing in metagenomic studies. *Nature Biotechnology* 2017;**35**:1069–1076.
3. Leek JT, Scharpf RB, Bravo HC, *et al.* Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics* 2010;**11**:733–739.
4. Duvallet C, Gibbons SM, Gurry T, Irizarry RA, Alm EJ. Meta-analysis of gut microbiome studies identifies disease-specific and shared responses. *Nature Communications* 2017;**8**:1784.
5. Thomas AM, Manghi P, Asnicar F, *et al.* Metagenomic analysis of colorectal cancer datasets identifies cross-cohort microbial diagnostic signatures and a link with choline degradation. *Nature Medicine* 2019;**25**:667–678.
6. Wirbel J, Zych K, Essex M, *et al.* Microbiome meta-analysis and cross-disease comparison enabled by the SIAMCAT machine learning toolbox. *Genome Biology* 2021;**22**:93.
7. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 2007;**8**:118–127.
8. Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics* 2020;**2**:lqaa078.
9. Ritchie ME, Phipson B, Wu D, *et al.* limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 2015;**43**:e47.
10. Korsunsky I, Millard N, Fan J, *et al.* Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* 2019;**16**:1289–1296.
11. Haghverdi L, Lun ATL, Morgan MD, Marioni JC. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nature Biotechnology* 2018;**36**:421–427.
12. Xu C, Lopez R, Mehlman E, Regier J, Jordan MI, Yosef N. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. *Molecular Systems Biology* 2021;**17**:e9620.
13. Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology* 2017;**8**:2224.
14. Ma S, Shungin D, Mallick H, *et al.* Population structure discovery in meta-analyzed microbial communities and inflammatory bowel disease. *Genome Biology* 2022;**23**:208.
15. Ling W, Lu J, Zhao N, *et al.* Batch effects removal for microbiome data via conditional quantile regression. *Nature Communications* 2022;**13**:5418.
16. Wang Y, Lê Cao K-A. PLSDA-batch: a multivariate framework to correct for batch effects in microbiome data. *Briefings in Bioinformatics* 2023;**24**:bbac622.
17. Gibbons SM, Duvallet C, Alm EJ. Correcting for batch effects in case-control microbiome studies. *PLoS Computational Biology* 2018;**14**:e1006102.
18. Austin GI, Brown Kav A, Park H, *et al.* Processing-bias correction with DEBIAS-M improves cross-study generalization of microbiome-based prediction models. *Nature Microbiology* 2025. ⚠️
19. MetaDICT — dictionary-learning batch correction and integration for microbiome data. ⚠️
20. Nearing JT, Douglas GM, Hayes MG, *et al.* Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications* 2022;**13**:342.
21. Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. The sva package for removing batch effects and other unwanted variation in high-throughput experiments. *Bioinformatics* 2012;**28**:882–883.
22. Park *et al.* Composite quantile regression for microbiome batch correction, 2025. ⚠️
23. Haghverdi L, Lun ATL, Morgan MD, Marioni JC; batchelor: single-cell batch correction methods. Bioconductor package.
24. Gayoso A, Lopez R, Xing G, *et al.* A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology* 2022;**40**:163–166.
25. Salim A, Molania R, Wang J, *et al.* RUV-III-NB: normalization of single cell RNA-seq data. *Nucleic Acids Research* 2022;**50**:e96.
26. Anderson MJ. A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 2001;**26**:32–46.
27. Oksanen J, Simpson GL, Blanchet FG, *et al.* vegan: Community Ecology Package. R package.
28. Aitchison J. *The Statistical Analysis of Compositional Data.* London: Chapman & Hall, 1986.
29. Bray JR, Curtis JT. An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs* 1957;**27**:325–349.
30. Rousseeuw PJ. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics* 1987;**20**:53–65.
31. McLaren MR, Willis AD, Callahan BJ. Consistent and correctable bias in metagenomic sequencing experiments. *eLife* 2019;**8**:e46923.
32. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* 1995;**57**:289–300.
33. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software* 2010;**33**:1–22.
34. Breiman L. Random forests. *Machine Learning* 2001;**45**:5–32.
35. Robin X, Turck N, Hainard A, *et al.* pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics* 2011;**12**:77.
36. Luecken MD, Büttner M, Chaichoompu K, *et al.* Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* 2022;**19**:41–50.
37. Cramér H. *Mathematical Methods of Statistics.* Princeton: Princeton University Press, 1946.
38. He Y, Zhao N, Satten GA. MIDASim: a fast and simple simulator for realistic microbiome data. *Microbiome* 2024;**12**:135. ⚠️
39. McKay MD, Beckman RJ, Conover WJ. A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. *Technometrics* 1979;**21**:239–245.
40. Mölder F, Jablonski KP, Letcher B, *et al.* Sustainable data analysis with Snakemake. *F1000Research* 2021;**10**:33.

Dataset sources are cited in [`methods_data.md`](methods_data.md) and are not repeated here.

**⚠️ = verify before submission.** Reference 19 (MetaDICT) and reference 22 (the composite
quantile-regression method, cited in the repository as "Park 2025") are cited in the code but I could not
confirm their bibliographic details; reference 18 (DEBIAS-M) and reference 38 (MIDASim) are cited from
memory and the volume, year or journal may be wrong. Two further citations already flagged in
`methods_data.md` — Foxx & Rivers (2025) and the MIDASim entry — are the same issue.

---

## Table 1. Methods under evaluation

| Method | Track | Output | Needs phenotype | Needs per-batch controls | Notes |
|---|---|---|---|---|---|
| Raw | baseline | counts | no | no | uncorrected input; reference for every change |
| Oracle | baseline | counts | no | no | batch-free truth; simulations only, upper bound |
| limma | A | clr | yes | no | `removeBatchEffect` on CLR |
| ComBat | A | clr | yes | no | `sva::ComBat` on CLR |
| ComBat-seq | A | counts | yes | no | `sva::ComBat_seq` |
| MMUPHin | A | relabund | yes | no | `adjust_batch` |
| ConQuR (tuned) | A | counts | yes | no | authors' tuning grid; >50% of total compute |
| ConQuR (default) | A | counts | yes | no | untuned, reference = first batch level |
| CQR | A | counts | yes | no | partial reimplementation; NB step omitted |
| PLSDA-batch | A | clr | yes | no | `balance = TRUE` |
| PLSDA-batch (weighted) | A | clr | yes | no | `balance = FALSE` |
| Percentile norm. | A | percentile | yes | **yes** | binary phenotype only; controls = `phenotype == 0` |
| MetaDICT | A | counts | yes | no | uses taxonomy or tree when present |
| DEBIAS-M | A | relabund | yes | no | binary phenotype only |
| Harmony | B | embedding | no | no | top 20 CLR principal components |
| fastMNN | B | embedding | no | no | `batchelor::fastMNN`, d = 20 |
| scANVI | B | embedding | yes | no | scvi-tools SCVI → SCANVI, defaults |
| ConQuR (permuted cov.) | control | counts | yes | no | negative control, not a candidate method |
| RUV-III-NB | — | — | — | — | wrapper kept; unavailable under current Bioconductor |
