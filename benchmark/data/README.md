# Datasets

Every dataset the benchmark uses, where it came from, and what was done to it here.

- **[`DATASETS.tsv`](DATASETS.tsv)** is the structured index: one row per dataset with its path, size, batch
  and phenotype variables, covariates, and caveats. Regenerate it with `data/inventory.sh` after rebuilding
  anything — shape is measured from the files, provenance is curated inside the script.
- This file explains *how each one was made* and *what to watch out for when reading results from it*.
- Everything rebuilds with `snakemake data`, except the three sets whose raw reads we processed ourselves
  (see [Rebuilding](#rebuilding)).

As of 2026-10-21: **10 evaluated datasets, 3 calibration sets, 3 simulation inputs, 8,670 samples.**

---

## Conventions

Each dataset is one directory under `results/data/real/<name>/` (goes through the method × dataset grid) or
`results/data/mock/<name>/` (calibration only, not in the run list).

| File | Contents |
|---|---|
| `counts.tsv` | rows = samples, first column `sample_id`, remaining columns = taxa, **integer read counts** |
| `meta.tsv` | first column `sample_id` in the same order; required columns `batch` and `phenotype`; any further columns are covariates |
| `taxonomy.tsv` | optional; rows = taxa matching the count columns, columns = ranks, missing ranks filled `unclassified` |
| `truth.tsv` | optional; rows = mock samples, columns = taxa, values = known proportions |

Shared filter unless a dataset's section says otherwise: **taxon prevalence ≥ 10%**, drop all-zero samples.
16S data is aggregated to **genus**, shotgun metagenomes to **species** (or strain, where the community is a
defined mock). Raw downloads live under `results/raw/<name>/` and are gitignored.

To add a dataset to the method evaluation, put its name in the `real:` list in `config/config.yaml`.

---

## Curated shotgun metagenomes — `crc_mgx`, `ibd_mgx`, `ici_mgx`

**Upstream** — curatedMetagenomicData 3.14.0 reprocesses every study through one pipeline, which is what
makes them comparable: raw WGS reads → **MetaPhlAn 3** (mpa_v30 ChocoPhlAn markers) → species-level relative
abundance, with metadata standardised by the cMD curators (`study_condition`, `subject_id`,
`sequencing_platform`, `DNA_extraction_kit`). HUMAnN 3 functional profiles exist but are not used here.

**Here** (`data/cmd.R`)

1. `counts = TRUE`, i.e. cMD's relative abundance ÷ 100 × the sample's read count, rounded.
2. Cohort and sample selection: CRC and IBD by `study_condition` (case vs control), ICI by `ORR`
   (responder vs non-responder).
3. One sample per subject — the earliest `days_from_first_collection`, which longitudinal cohorts like HMP2
   require.
4. `platform` and `extraction_kit` are kept as covariates; DEBIAS-M found the extraction kit explains most
   processing bias. Taxonomy comes from `rowData`.

> **These are not sequencing counts.** They are back-calculated from relative abundance, which is a known
> approximation for methods that assume a true count distribution (ComBat-seq, ConQuR, RUV-III-NB). Say so
> in the paper.

> In `crc_mgx` the **extraction kit is perfectly nested within study** (Gnome ×3 cohorts, MoBio ×2,
> Qiagen ×2, YachidaS_2019 unrecorded), so kit and study effects cannot be separated. The kit can still be
> used to *group* studies as an alternative batch definition, which tests DEBIAS-M's central claim that the
> bias is driven by the extraction protocol.

---

## Curated 16S — `crc_16s`, `ibd_16s`, `hiv_16s`

**Upstream** — MicrobiomeHD (Duvallet et al. 2017, Zenodo 569601) reprocesses 28 case-control studies through
one pipeline: quality filtering and dereplication → **100% identity de novo OTU clustering**
(`otu_table.100.denovo`) → **RDP classifier** assignment, with the full lineage as the row name.

**Here** (`data/microbiomehd.R`)

1. Takes `*.otu_table.100.denovo.rdp_assigned` per study (some studies also ship a dbOTU table; unused).
2. Per-study sample filters copied from the authors' own `dataset_info.yaml`: `hiv_lozupone` keeps
   `time_point = 1`, `hiv_noguerajulian` keeps `cohort ∈ {BCN0, STK}`; the `ibd_gevers_2014` Zenodo table
   already holds stool only.
3. Case/control mapping: CRC vs H (adenoma `nonCRC` dropped); CD+UC vs H/nonIBD; HIV vs H.
4. The authors' floor of ≥ 100 reads per sample.
5. **Aggregated to genus via the RDP lineage** — OTUs are clustered de novo *per study*, so their numbering
   means nothing across studies and only genus is comparable.

> The studies differ in 16S region, primers and platform (`hiv_dinh` is 454 V3-V5, `crc_baxter` is MiSeq V4).
> That is precisely the batch effect under test.

---

## `hivrc_16s` — 17 HIV studies

**Upstream** — the HIV Re-analysis Consortium (Tuddenham et al. 2020 CID; Synapse syn18406854) reprocessed
**17 of 22** HIV gut microbiome studies through one pipeline: raw 16S reads → **Resphera Insight** for
species-level assignment. Two tables were published: species-level raw counts (`insight.merged_otus.txt`) and
a genus table (`insight.merged_genus.txt`) **rarefied to exactly 20,000 reads per sample**.

Synapse needs an account, so the three source files were downloaded by hand into `results/raw/hivrc/`:
`insight.merged_genus.txt` (genus counts, 1,026 × 1,032), `insight.merged_otus.txt` (species-level raw
counts, 38,124 rows, unused), `SupplementaryMaterial.xlsx` (sample manifest and per-study accessions).

**Here** (`data/hivrc_16s.R`) joins the genus table to the manifest on `SeqID` (all 1,032 samples match),
sets `batch` = study, `phenotype` = `hivstatus`, and keeps Age / gender / msm / ARTuse / vl / cd4 as
covariates. Result: **1,032 samples × 160 genera across 17 studies**.

> **Naturally confounded in part.** Perez-Santiago (0/13) and Villar-Garcia (0/24) have **no HIV-negative
> controls** and Serrano-Villar 2017 has two. Percentile Normalization, which needs controls in every batch,
> fails here — a real instance of the confounded-design problem from the literature.

> The genus table is rarefied, so sequencing-depth variation is already gone; methods that assume raw counts
> receive rarefied ones. The unrarefied species table is on the server if depth needs to be preserved.

> Strongest uncorrected batch effect of all our real data: batch R² (Aitchison) = **0.305**, matching the
> 0.12–0.30 published for HIVRC.

---

## `seed_16s` — the confounded case

**Upstream** — Foxx & Rivers 2025 (Mendeley Data 10.17632/5xrfg5dym6.1): 8 seed-microbiome studies run
through QIIME2 to an **ASV table**, published as `F&R_original_count.csv` plus metadata.

**Here** (`data/seed_16s.R`) uses the published count table, with `batch` = study, `phenotype` = host plant
(**multi-class**, not case-control), and `region` kept as a covariate.

> Prevalence is filtered at ≥ 10% **within at least one study** rather than across all samples: the studies
> target different regions (V4 / V3-V4 / V5-V7 / V6-V8) so their ASVs barely overlap, and the uniform rule
> would leave only 2 of 9 studies standing.

> ASVs are sequence hashes with no taxonomy, so there is no `taxonomy.tsv`. Barret 2015 and Liu 2019 are
> absent from the published count table.

---

## `mbqc` — 13 handling labs × 22 specimens

**Upstream** — MBQC (Microbiome Quality Control, Sinha et al. 2017) split one set of specimens and sent
aliquots to several "handling labs", each doing its own DNA extraction and 16S V4 sequencing. **The
consortium's integrated OTU table is gone** — the server that held it is offline and the Wayback Machine has
no copy — so only the raw reads survive (SRA SRP047083). We rebuilt the dataset from them.

**Here** (`data/mbqc_manifest.sh` → `data/mbqc_download.sh` → `data/mbqc.R`)

1. The design comes back out of ENA: `library_name` is `<specimen>.<handling lab>.<id>`, which gives 13 labs
   × 22 specimens. Specimen attributes (Fresh / freeze-dried / Robogut / artificial colony, health status,
   sex, age) come from the BioSample records.
2. **Forward reads only** (2,117 files, 13 GB). The labs ran anything from single-end to 2×294, so reverse
   reads cannot be overlapped consistently.
3. **Each lab is denoised separately**, which is also DADA2's own rule for separate runs. `truncLen` is
   chosen from the quality profile (last cycle whose median Q ≥ 25) and **capped at the 10th percentile of
   read length**; `maxEE = 2` everywhere, and a lab keeping under 50% of its reads is retried with a shorter
   `truncLen` rather than a looser error budget.
4. Per-lab `learnErrors` / `dada` / chimera removal. ASVs differ in length between labs and cannot be merged
   as sequences, so everything is **aggregated to genus**, as for MicrobiomeHD.
5. Taxonomy is assigned once over the union of 49,765 ASVs with our SILVA 138.2 training set (~1 hour).

Result: **2,034 samples × 181 genera**, median 24,171 reads per sample, 259 artificial-community samples.

| Lab | read length | truncLen | reads passing filter | samples |
|---|---|---|---|---|
| agoodman | 251 | 251 | 83% | 76 |
| dgevers | 175 | 174 | 81% | 19 |
| dlittman | 151 | 151 | 99% | 90 |
| dmills | 243 | 243 | 96% | 172 |
| ggloor | 253 | **160** (shortened twice) | 53% | 87 |
| jpetrosino | 250 | 250 | 64% | 462 |
| jravel | 294 | 293 | 63% | 66 |
| kjones | 151 | 151 | 97% | 184 |
| pschloss | 251 | 251 | 84% | 270 |
| pturnbaugh | 150 | 100 | 74% | 92 |
| rburk | 251 | 249 | 97% | 90 |
| rflores | 251 | 251 | 89% | 179 |
| rknight | 152 | 152 | 98% | 247 |

> **Fully crossed**: every one of the 22 specimens went to all 13 labs, so batch and phenotype are completely
> unconfounded (Cramér's V = 0). This is the zero point of the confounding gradient, and the only dataset
> here that provides one. Batch R² (Aitchison) = **0.099**, specimen R² = **0.481**; mean within-specimen
> distance 26.0 against 38.3 between specimens, over 100,960 within-specimen pairs.

> **No `truth.tsv`.** The two artificial-community specimens (sample103 oral, sample104 faecal) have a
> published composition given only as "approximate loopfuls", not quantitative proportions, and the public
> metadata has no mapping from specimen numbers to MBQC's internal IDs. The 259 artificial-community samples
> are flagged `artificial = 1` for anyone who wants them.

> Six FastQ files are malformed at source and were skipped: SRR1591126, SRR1590533, SRR1590746, SRR1590112,
> SRR1591234, SRR1589983.

---

## `jmbc_extraction` and `jmbc_library` — protocol comparisons

**Upstream** — Tourlousse et al. 2021, *Microbiome* 9:95 (Japanese Microbiome Consortium, SRA
PRJNA650228). The paper ships **one 54-page PDF and nothing else**: no feature table anywhere, and ENA's
metadata is no help either (`library_name` is a random code, the BioSample attributes are placeholders). The
MOSAIC Standards Challenge platform, where the data was also deposited, no longer resolves. Rebuilt from raw
reads.

**Here** (`data/jmbc_tables.py` → `jmbc_ref.sh` → `jmbc_download.sh` → `jmbc_quantify.sh` → `jmbc.R`)

1. `pdftotext -layout` pulls two tables out of the PDF; both are committed under `data/jmbc/`, so the cluster
   needs no poppler:
   - **Table S1** → the 20 mock strains with GenBank assembly accessions, genome sizes, Gram type and **true
     proportions** (DNA mock even at 5%; cell mock uneven, with *Megamonas funiformis* at 0).
   - **Table S8** → which experiment and protocol each run belongs to. pdftotext occasionally splits a
     two-character protocol code across the column edge (`C protocol L` = CL); the parser puts it back.
2. Only the two phase-I protocol comparisons are kept (164 runs). Phase II varies operator, reagent lot and
   participating lab instead.
3. **Reference**: whole assemblies pulled from ENA by GCA accession (chromosomes and plasmids), 20 genomes
   totalling ~75 Mb, indexed with bowtie2. The community's membership is known, so no MetaPhlAn or GTDB
   database is needed.
4. **Only the first 200,000 reads of each run** are downloaded — ENA honours HTTP Range, which turns 100 GB
   into **2.3 GB** and is still far more depth than a 20-strain mock needs.
5. Only reads with **MAPQ ≥ 10** are counted: the two *B. longum* subspecies produce alignments that cannot
   be told apart, and dropping them beats letting bowtie2 break the tie at random.

| Dataset | Size | Batch | Role |
|---|---|---|---|
| `jmbc_extraction` | 79 samples × 20 strains (54 cell mock, 25 faeces) | **9 DNA extraction protocols** | evaluated |
| `jmbc_library` | 84 samples × 20 strains | 28 library preparations (11 kits × 3 input amounts) | calibration |

> **`jmbc_extraction` is the only dataset in the benchmark whose batch variable is the DNA extraction
> method**, which the literature treats as the dominant source of protocol bias. Against the known mock,
> per-taxon log2 bias spans **−2.09 to +1.50 (~12-fold)** in the direction of how hard each strain is to
> lyse: *Bifidobacterium* (both subspecies), *Lactobacillus* and *Staphylococcus* strongly under-detected,
> *Blautia*, *Bacillus* and *Parabacteroides* over-detected. Batch R² = **0.143**, material R² = **0.976**.

> **A result that shaped the simulation design.** The bias vector agrees closely across seven of the nine
> extraction protocols (r ≈ 0.92–0.98) but **protocols O and R do not** (r ≈ 0.56–0.84) — whereas
> `salim_pig`, which only varies library preparation, agrees at r = 0.992–0.993. So the McLaren
> one-multiplicative-factor-per-taxon model is solid at the library step but **mostly, not universally**,
> true at the extraction step. That is the empirical basis for simulating both systematic and non-systematic
> batch effects (PLAN.md §2.1).

> `jmbc_library` is 84 preparations of **one** DNA mock: nothing biological varies, so it is a pure null set.
> Any structure a method leaves behind is batch, and any difference it creates is a false signal. Calibration
> only.

> 163 of 164 runs succeeded; SRR12996048 came back empty on every retry and was skipped.

---

## `salim_pig` — storage and handling (calibration)

**Upstream** — Poulsen et al. 2021, *Microbiol Spectr* 9:e01387-21 (DTU, Pamp group; ENA PRJEB31650). Four
materials (pig faeces ×2, sewage ×2) were aliquoted and held under 19 time × temperature combinations
(0 h–8,640 h at −80/−20/5/22 °C) with 0–4 freeze-thaw cycles and 2–3 technical replicates each; half the
aliquots were spiked with an 8-strain mock, plus 24 pure-mock and 30 blank extractions. Sequenced on HiSeq
4000 with **3 library kits** (Nextflex / Kappa / Nextera) across 12 flow-cell lanes. BBduk2 → **MGmapper**
(BWA-MEM against 11 reference databases) → counts halved (proper pairs) and aggregated to genus, published on
figshare 10.6084/m9.figshare.12010989.

Fachrul et al. 2022 *Sci Rep* — the "Salim 2022" entry in earlier notes — is a **methods paper reusing this
same data**, which is why the raw reads looked unavoidable at first. They are not: the count table is 2.6 MB.

Sources, all fetched automatically by the build script:

| Source | File | Size | Contents |
|---|---|---|---|
| figshare | `Item_A_TaxonomyRaw.txt` | 2.6 MB | genus-level **raw integer counts**, 2,640 × 438 |
| figshare | `Item_C_FeatureShort.txt` | 229 KB | Domain…Genus lineage per row |
| Europe PMC | `spectrum.01387-21-s0001.xlsx` (Table S1) | 185 KB | per-sample conditions |

Table S1 has to come from **Europe PMC**: PMC itself answers scripted requests with a proof-of-work page and
ASM returns 403.

**Here** (`data/salim_pig.R`) joins the count table to Table S1 on `Sample`, drops the 30 blank extractions
and one sample from a later run absent from Table S1, keeps **Bacteria + Archaea plus the mock's two
eukaryotes** (*Saccharomyces*, *Cryptosporidium*) since MGmapper also maps human, animal and plant
references. `batch` = storage condition, `phenotype` = material, `replicate` = replicate group, with
`spiked`, `library_prep`, `lane`, `time_h`, `temp_c` and `freeze_thaw` kept so any of them can serve as the
batch instead. Result: **408 samples × 1,537 genera, 19 storage conditions, 125 replicate groups**.

`truth.tsv` holds the known **cell proportions** of the 24 pure-mock samples (Table S3; Pig_feces_1 received
half as much *Propionibacterium* and *Staphylococcus* as the other three materials). As with `brooks`,
genome-size effects stay inside the measured bias, which is exactly what the McLaren model describes.

> **This is what validates the simulation's batch model.** Per-taxon log2 bias in the pure-mock samples:
>
> | Over-detected | | Under-detected | |
> |---|---|---|---|
> | *Salmonella* | +2.34 | *Cryptosporidium* | −3.41 |
> | *Bacteroides* | +1.95 | *Staphylococcus* | −4.00 |
> | *Escherichia* | +1.61 | *Saccharomyces* | −5.88 |
> | *Fusobacterium* | +1.58 | *Propionibacterium* | −6.19 |
>
> A ~360-fold span, in the direction of lysis difficulty. The key point: **the bias vector is near-identical
> across the three library kits (r = 0.992–0.993)**, so one multiplicative factor per taxon holds, and the
> bias is set by extraction and lysis rather than by library preparation — DEBIAS-M's central claim. Genome
> size accounts for only about one log2 unit of the spread.

> Measured: material R² = 0.640 (Aitchison) / 0.692 (Bray) — pig faeces and sewage are entirely different
> communities, so over-correction has something obvious to destroy. Technical batches are milder: storage
> 0.068, library prep 0.032, lane 0.039. Within-replicate-group Aitchison distance 35.4 against 67.2 between
> groups, so "replicates should agree" is usable as truth here.

> These are MGmapper reference-database counts, not MetaPhlAn marker abundances, so they are not directly
> comparable with the cMD datasets. Calibration only. The authors also publish a genome-size-normalised,
> total-sum-scaled table (figshare `Item_B`) which we do not use.

---

## `brooks` — vaginal mock (calibration)

**Upstream** — Brooks et al. 2015 cell mocks, bundled with the `metacal` package by McLaren et al. 2019
(eLife): seven vaginal species mixed in **known proportions**.

**Here** (`data/brooks_mock.R`): `counts.tsv` = observed reads, `truth.tsv` = actual proportions, `batch` =
plate, `phenotype` = how many species went into the mixture (1/2/3/4/7), with `species_list` alongside.

> **Two plates and 80 samples only** — metacal ships the cell-mixture subset, not all of Brooks 2015. (Earlier
> notes in this repo said six plates; that was the paper's scale, not the bundled data.) `Mixture_type` is
> "Cells" throughout and is dropped. With only two batches the per-batch bias estimate is weak; `salim_pig`
> (3 library kits × 24 pure mocks) is the stronger version of the same check on the shotgun side.

---

## Simulation inputs

| Input | Built by | Contents |
|---|---|---|
| `data/templates/16s.tsv` | `simulate/make_templates.R` | 146 × 614, from HMP/AGP healthy gut |
| `data/templates/mgx.tsv` | `simulate/make_templates.R` | 624 × 224, from the healthy samples of `crc_mgx` |
| `config/scenarios.tsv` | `simulate/design.R` | 213 scenarios × 12 parameters |

MIDASim is fitted to each template once (`<template>.midasim.rds`); the `sim` rule then generates
`results/data/sim/<scenario>/rep<N>/` with `counts.tsv`, `meta.tsv`, `oracle.tsv`, `truth_da.txt`,
`truth_fc.tsv` and `params.json`. Templates are gitignored because they are derived.

---

## Rebuilding

`snakemake data` rebuilds everything whose inputs are fetched automatically. Three datasets need a step
first, because they start from raw reads or a manual download:

| Dataset | Prerequisite |
|---|---|
| `hivrc_16s` | the three Synapse files in `results/raw/hivrc/` (account required) |
| `mbqc` | `data/mbqc_manifest.sh` then `data/mbqc_download.sh` (13 GB, several hours) |
| `jmbc_*` | `data/jmbc_ref.sh`, `data/jmbc_download.sh`, `data/jmbc_quantify.sh` (2.3 GB) |

`mbqc` and `jmbc` are therefore not part of the default `data` target, so a stray `snakemake data` cannot
kick off a multi-hour download. The SILVA 138.2 DADA2 training set that `mbqc` needs is built once by
`data/silva_dada2_trainset.sh` from the release files on `/beegfs`.

The MBQC build keeps per-lab checkpoints (`results/raw/mbqc/checkpoints/`) and holds a lock file, so it can
be resumed and cannot race a second copy of itself on another node.

---

## What is kept on the cluster

The prepared datasets are committed to this repository (`results/data/real/` and `results/data/mock/`,
6.4 MB), so the cluster no longer holds the only copy. Its `results/raw/` was cut from 25 GB to 113 MB on
2026-10-21; everything removed is reproducible from the scripts here.

| Kept | Size | Why it is not deleted |
|---|---|---|
| `raw/hivrc/` | 94 MB | Synapse needs an account and these were downloaded by hand; re-fetching needs a person, not a script |
| `raw/seed/` | 11 MB | Small enough that deleting it buys nothing |
| `raw/mbqc/checkpoints/` + manifests | 4.3 MB | **Rebuilds MBQC in 17 seconds instead of a 13 GB download and eleven hours of DADA2** |
| `raw/salim/` (3 source files) | 3.0 MB | Table S1 was uploaded by hand because Europe PMC is slow and flaky |
| `raw/jmbc/counts_long.tsv` | 104 KB | Rebuilds both JMBC datasets without re-downloading 2.5 GB |

Deleted: `raw/mbqc/fastq` (12 GB), `raw/mbqc/filtered` (8.3 GB, a DADA2 intermediate), `raw/jmbc/reads`
(2.5 GB), `raw/microbiomehd` (1.4 GB, re-downloads from Zenodo), `raw/jmbc/ref` (245 MB, rebuilt by
`jmbc_ref.sh`), and assorted leftovers from abandoned attempts.

The trade is the point: **4.4 MB of intermediates replaces 23 GB of raw reads and still rebuilds in
seconds.** Both paths were tested after the deletion — MBQC and both JMBC datasets rebuild byte-identical
from what remains.

> The MBQC build holds `results/raw/mbqc/checkpoints/.lock` while it runs and releases it on completion. A
> build killed partway leaves the lock behind — R has no reliable at-exit hook for a script stopped by an
> error — so it is ignored once three hours old, or delete it by hand.

---

## What we could not get

Nothing outstanding. Two entries that were open for a long time are closed:

- **MBQC** — the integrated OTU table is permanently gone, but the raw data turned out to be small (2,117
  runs, 28 GB total, already demultiplexed, lab encoded in `library_name`), so it was rebuilt.
- **Tourlousse JMBC** — no feature table exists anywhere, but the design is recoverable from the
  supplementary PDF and partial downloads keep the cost at 2.3 GB, so it was rebuilt.

**Measured bandwidth from the compute nodes** (2026-09-27), for sizing any future download: ENA https
≈ 242 KB/s, ENA FTP ≈ 153 KB/s, Zenodo ≈ 80 KB/s, NCBI S3 ≈ 27 KB/s. At the best of those, a full
reprocessing of a 461 GB study would take about 22 days — which is why partial downloads and published
tables are preferred throughout.
