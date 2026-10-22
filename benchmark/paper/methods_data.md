# Methods — Benchmark datasets

## Overview

Evaluating batch-correction methods requires three kinds of data, because no single kind can answer every
question. Real case-control studies show whether a method preserves a biological contrast while removing
study-to-study variation, but they have no ground truth. Designed protocol comparisons, in which the same
biological material is processed by several laboratory protocols, isolate a named technical factor and —
when the material is a defined mock community — provide the true composition. Simulations provide complete
ground truth over a controlled grid, but only inherit realism from the templates they are fitted to. We
therefore assembled all three.

The benchmark comprises **7,687 biological samples in 13 datasets** — ten entering the method comparison
(7,115 samples) and three reserved for calibration (572 samples) — together with 213 simulation scenarios.
Seven datasets are 16S rRNA amplicon surveys and six are shotgun metagenomes. Supplementary Table S1
summarises every dataset; the sections below describe how each was obtained and processed.

Two of the thirteen were reconstructed from raw sequencing data by us because no processed feature table
exists for them, and these are described in the greatest detail. All processing scripts, the assembled
tables and the exact parameters are available in the project repository (see *Data and code availability*).

## Datasets reanalysed from curated resources

Four groups of datasets were derived from resources in which the original authors had already reprocessed
every constituent study through a single pipeline. This is a prerequisite: where studies are reprocessed
uniformly, the residual between-study variation is attributable to the wet-laboratory and sequencing
differences we wish to study, rather than to bioinformatic choices.

**Shotgun metagenomes (`crc_mgx`, `ibd_mgx`, `ici_mgx`).** Species-level profiles were obtained from
curatedMetagenomicData 3.14.0 (Pasolli et al., 2017), in which raw whole-genome shotgun reads from every
study are processed with MetaPhlAn 3 against the mpa_v30 ChocoPhlAn marker database. We requested counts
(`counts = TRUE`), which curatedMetagenomicData derives as the relative abundance of each species
multiplied by the sample's read count. Cohorts were selected on the curated `study_condition` field for
colorectal cancer and inflammatory bowel disease and on objective response rate (`ORR`) for the
immune-checkpoint-inhibitor cohorts. Where a subject contributed several samples, the earliest
(`days_from_first_collection`) was retained, so that batch structure is not confounded with repeated
measures. Sequencing platform and DNA extraction kit were retained as covariates.

**16S case-control studies (`crc_16s`, `ibd_16s`, `hiv_16s`).** Operational taxonomic unit tables were taken
from MicrobiomeHD (Duvallet et al., 2017), which reprocesses 28 case-control studies by quality filtering,
dereplication, de novo clustering at 100% identity and RDP classification (Wang et al., 2007). Because OTUs
are clustered de novo within each study, their identifiers are not comparable across studies; we therefore
aggregated counts to genus using the RDP lineage. We adopted the original authors' per-study sample filters
and their floor of 100 reads per sample. Case and control definitions follow the MicrobiomeHD conventions,
with colorectal adenoma samples excluded from the cancer contrast.

**HIV re-analysis consortium (`hivrc_16s`).** Genus-level profiles for 17 of the 22 studies assembled by the
HIV Re-analysis Consortium (Tuddenham et al., 2020) were obtained from Synapse (syn18406854). The consortium
annotated all studies with a single Resphera Insight pipeline and published a genus table rarefied to 20,000
reads per sample, which we joined to the sample manifest on sequencing identifier; all 1,032 samples
matched. Age, sex, men-who-have-sex-with-men status, antiretroviral therapy use, viral load and CD4 count
were retained as covariates.

**Seed microbiome (`seed_16s`).** Amplicon sequence variant counts from eight seed-microbiome studies
processed in QIIME2 were taken from the published tables of Foxx and Rivers (2025). The studies target
different hypervariable regions (V4, V3–V4, V5–V7, V6–V8), so their ASVs scarcely overlap; applying our
standard prevalence filter across all samples would have retained only two of nine studies. We therefore
required an ASV to reach 10% prevalence **within at least one study**. This dataset is retained precisely
because its host-plant phenotype is strongly confounded with study, representing the confounded designs that
are common in practice.

## Datasets reconstructed from raw sequencing data

**Multi-laboratory amplicon data (`mbqc`).** The Microbiome Quality Control project (Sinha et al., 2017)
distributed aliquots of a common set of specimens to multiple handling laboratories, each of which performed
its own DNA extraction and 16S V4 sequencing. The consortium's integrated OTU table is no longer retrievable:
the server that hosted it is offline and no archived copy exists. We therefore rebuilt the dataset from the
raw reads (SRA SRP047083).

The experimental design was recovered from the sequence archive metadata, in which the library name encodes
specimen and handling laboratory, yielding 13 laboratories and 22 specimens; specimen attributes were taken
from the corresponding BioSample records. Because the laboratories used read configurations ranging from
single-end to 2 × 294 nt, reverse reads could not be overlapped consistently and only forward reads were
used (2,117 files).

Each laboratory was denoised separately with DADA2 1.22.0 (Callahan et al., 2016) under R 4.1.3, following
DADA2's own recommendation that error models be learned per sequencing run. Truncation length was selected
per laboratory as the last cycle at which median quality remained at or above Q25, capped at the 10th
percentile of observed read length; the cap is necessary because `filterAndTrim` discards any read shorter
than the truncation length, and the quality matrix pads shorter reads with missing values so that trailing
cycles are informed only by the longest reads. The expected-error threshold was held at `maxEE = 2` for every
laboratory; where this retained less than half of a laboratory's reads, the truncation length was shortened
in 20% steps rather than relaxing the error budget, so that the error criterion remains identical across
batches. Six FastQ files were malformed at source and were skipped.

Amplicon sequence variants were inferred and chimeras removed within each laboratory. Because truncation
lengths differ, ASVs are not comparable as sequences between laboratories, and counts were therefore
aggregated to genus, as for the MicrobiomeHD datasets. Taxonomy was assigned once over the union of 49,765
ASVs using a DADA2-formatted training set that we built from the SILVA 138.2 SSU Ref NR99 release (Quast
et al., 2013); the release taxonomy map assigns every bacterial and archaeal sequence a path of exactly six
ranks, which is the format `assignTaxonomy` expects, so 451,706 sequences were used without further
curation. The final dataset comprises 2,034 samples and 181 genera.

This dataset is **fully crossed**: each of the 22 specimens was processed by all 13 laboratories, so batch
and biological variable are completely unconfounded (Cramér's V = 0). It therefore provides the zero point
of the confounding gradient, and 100,960 within-specimen sample pairs serve as a replicate-based criterion
("replicates of one specimen should agree") in place of a known composition.

**Protocol comparison with defined mocks (`jmbc_extraction`, `jmbc_library`).** Tourlousse et al. (2021)
compared DNA extraction and library construction protocols using defined mock communities, but published no
feature table; the only supplementary file is a document, and the archive metadata contains no protocol
information. We therefore reconstructed the data from raw reads (SRA PRJNA650228).

The study design was extracted from the supplementary document: the mock community composition, comprising
20 strains with assembly accessions, genome sizes, Gram type and nominal proportions, and the assignment of
each sequencing run to an experiment and protocol. We retained the two phase-I protocol comparisons (164
runs); the phase-II experiments vary operator, reagent lot and participating laboratory instead.

Because the community membership is known, quantification required no metagenomic reference database. The
20 assemblies were retrieved from the European Nucleotide Archive and indexed with Bowtie 2 2.5.5 (Langmead
and Salzberg, 2012), giving a 75 Mb reference. Likewise, full sequencing depth was unnecessary: the archive
honours HTTP range requests, so only the leading 200,000 reads of each run were retrieved, reducing the
transfer from approximately 100 GB to 2.3 GB. Reads were aligned to the mock reference and counted per
strain using SAMtools 1.21 (Danecek et al., 2021), retaining only alignments with mapping quality of at
least 10; the community contains two *Bifidobacterium longum* subspecies whose reads cannot be assigned
unambiguously, and discarding them is preferable to resolving ties arbitrarily. One run returned no data and
was excluded.

The two experiments were kept as separate datasets because their batch variables are different quantities.
`jmbc_extraction` (79 samples, 9 extraction protocols) contrasts cell mock with faecal material and is the
only dataset in the benchmark whose batch variable is the DNA extraction method, which is widely regarded as
the dominant source of protocol bias. `jmbc_library` (84 samples, 28 library preparations of a single DNA
mock) contains no biological variation whatsoever and is therefore a null dataset: any structure a method
leaves behind is batch, and any difference it creates is a false signal.

## Mock communities with known composition

Three datasets provide a known composition against which measurement bias can be quantified directly, and
are used to calibrate the simulation rather than to rank methods.

`brooks` comprises 80 samples of seven vaginal species mixed in known proportions across two sequencing
plates, as distributed with the metacal package (Brooks et al., 2015; McLaren et al., 2019). `salim_pig`
comprises 408 samples of pig faeces and sewage held under 19 storage conditions with 125 technical replicate
groups and 24 pure mock-community samples; genus-level counts generated with MGmapper were obtained from the
authors' repository (Poulsen et al., 2021), and the known cell concentrations of the spiked community serve
as truth. `jmbc_library` is described above.

For both mock datasets the truth table records **cell proportions**, so that genome-size effects remain
part of the measured bias — which is what the multiplicative bias model of McLaren et al. (2019) is intended
to describe.

As a check on the reconstruction, observed compositions were compared with the nominal compositions. In
`jmbc_extraction`, per-taxon bias spanned 12-fold (log2 −2.09 to +1.50) and ordered taxa by the difficulty
of cell lysis, with *Bifidobacterium*, *Lactobacillus* and *Staphylococcus* strongly under-detected. In
`salim_pig` the corresponding span was approximately 360-fold. These magnitudes were used to bound the
per-taxon bias parameter of the simulation.

## Uniform preprocessing

All datasets were brought to a common format: a sample-by-taxon matrix of integer counts, a sample table
carrying a `batch` and a `phenotype` column plus any dataset-specific covariates, and, where available, a
taxonomy table and a truth table. Amplicon data were aggregated to genus and shotgun data to species, except
for the defined mocks, which are reported at strain level. Unless stated otherwise, taxa present in fewer
than 10% of samples were removed and samples with no remaining counts were dropped. No rarefaction,
transformation or compositional normalisation was applied at this stage; such choices belong to the methods
under evaluation and are applied by their respective wrappers.

## Simulated data

Simulated datasets were generated with MIDASim 2.0 in parametric mode, fitted separately to a 16S template
and a shotgun template so that sparsity and taxon-taxon correlation reflect real data of each type.
Differentially abundant taxa were implanted at known proportions and effect sizes, and batch effects were
superimposed under five models: per-taxon multiplicative bias, log-additive and scale shifts, batch-specific
zero inflation, non-systematic effects restricted to a subset of taxa, and sequencing-depth differences.
Batch-phenotype association was varied across a gradient of Cramér's V to make confounding an explicit
experimental factor rather than an incidental property. The grid comprises 213 scenarios, each replicated
20 times with fixed seeds. Null scenarios containing no biological signal were included to quantify the rate
at which each method manufactures one.

## Known limitations of the assembled data

Several properties of the source data constrain interpretation and are stated here rather than in the
discussion of individual results.

The curatedMetagenomicData counts are not sequencing counts but are derived from relative abundances, which
is an approximation for methods that assume a discrete count distribution. In `hivrc_16s` the published
genus table is rarefied, so between-sample depth variation has already been removed and methods that model
depth receive rarefied input. In `crc_mgx` the DNA extraction kit is perfectly nested within study, so kit
and study effects cannot be separated, although the kit can be used to group studies as an alternative batch
definition. The metacal distribution of the Brooks mock contains two plates rather than the full study,
limiting the precision with which per-batch bias can be estimated from it. For MBQC, the composition of the
two artificial-community specimens is published only in approximate terms and no mapping from public
specimen identifiers to consortium identifiers is available, so this dataset contributes replicate-based
rather than composition-based truth.

## Data and code availability

All datasets are derived from public sources; accessions are listed in Supplementary Table S1. The HIV
Re-analysis Consortium data require a Synapse account. The scripts that assemble every dataset, the
assembled tables, and a machine-readable dataset index are available at
`https://github.com/delect-aj/Benchmark_batch_effect`. Datasets reconstructed from raw reads can be rebuilt
from the retained intermediate files without repeating the downloads.

## Supplementary Table S1

**Table S1. Datasets used in the benchmark.** One row per dataset (`paper/table_s1_datasets.tsv`). *Role*
distinguishes datasets entering the method comparison from those reserved for calibration; `jmbc_library`
is further marked as a null dataset because it contains no biological variation. *Ground truth* states what
each dataset can be scored against: a known community composition, agreement between technical replicates of
one specimen, or neither, in which case only proxy criteria (variance attributable to batch and to the
biological variable) apply. Batch and biological variables are named as they are modelled here; several
datasets carry additional covariates, listed in the repository index. Sample counts are after the
preprocessing described above. Simulated data are not listed: they comprise 213 scenarios, each replicated
20 times, generated from a 16S and a shotgun template as described in *Simulated data*.

## References cited in this section

Brooks JP, et al. *BMC Microbiol* 2015;15:66. · Callahan BJ, et al. *Nat Methods* 2016;13:581–583. ·
Danecek P, et al. *GigaScience* 2021;10:giab008. · Duvallet C, et al. *Nat Commun* 2017;8:1784. ·
Foxx AJ, Rivers A. Mendeley Data 2025, doi:10.17632/5xrfg5dym6.1. · Langmead B, Salzberg SL. *Nat Methods*
2012;9:357–359. · McLaren MR, Willis AD, Callahan BJ. *eLife* 2019;8:e46923. · Pasolli E, et al. *Nat
Methods* 2017;14:1023–1024. · Poulsen CS, et al. *Microbiol Spectr* 2021;9:e01387-21. · Quast C, et al.
*Nucleic Acids Res* 2013;41:D590–D596. · Sinha R, et al. *Nat Biotechnol* 2017;35:1077–1086. · Tourlousse
DM, et al. *Microbiome* 2021;9:95. · Tuddenham SA, et al. *Clin Infect Dis* 2020;70:615–627. · Wang Q,
et al. *Appl Environ Microbiol* 2007;73:5261–5267.
