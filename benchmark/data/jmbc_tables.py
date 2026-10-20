# Pull the two tables we need out of the JMBC supplementary PDF (Tourlousse et al. 2021, Microbiome 9:95).
# Nothing else in that paper is machine readable: the only additional file is this 54-page PDF, and ENA's own
# metadata is useless here (library_name is a random code, the BioSample attributes carry no protocol).
#   Table S1 -> the mock community: 20 strains with assembly accessions and their true DNA/cell proportions
#   Table S8 -> which protocol each run belongs to, merged with ENA's file URLs and sizes
# The outputs are committed, so the cluster does not need poppler to rebuild the dataset.
# Usage: python3 data/jmbc_tables.py <supplementary pdf> <outdir>     (needs pdftotext from poppler)
import io
import re
import subprocess
import sys

pdf, outdir = sys.argv[1], sys.argv[2]
txt = subprocess.run(["pdftotext", "-layout", "-f", "38", "-l", "54", pdf, "-"],
                     capture_output=True, text=True, check=True).stdout

# --- Table S1: one row per strain, with continuation lines for extra accessions ---
s1 = txt[txt.index("Supplementary Table S1"):txt.index("Supplementary Table S2")]
mock = []
for line in s1.split("\n"):
    m = re.match(r"\s{2,}(\S.*?)\s{2,}(\S.*?)\s{2,}(\S+)\s{2,}(GCA_\S+)\s{2,}([\d,]+)\s{2,}([\d.]+)\s{2,}"
                 r"(\S+)\s{2,}([\d.]+)\s{2,}([\d.]+)", line)
    if m:
        sp, strain, _, gca, size, gc, gram, dna, cell = m.groups()
        mock.append([sp.strip(), strain.strip(), gca, size.replace(",", ""), gc, gram, dna, cell])
hdr = ["species", "strain", "assembly", "genome_size", "gc", "gram", "pct_dna_mock", "pct_cell_mock"]
io.open(f"{outdir}/jmbc_mock.tsv", "w", encoding="utf-8").write(
    "\n".join("\t".join(r) for r in [hdr] + mock) + "\n")
print(f"Table S1: {len(mock)} strains")

# --- Table S8: SRA accession, free-text metafields, read pairs; grouped by the experiment above them ---
s8 = txt[txt.index("Supplementary Table S8"):]
runs, exp = {}, None
for chunk in s8.split("Supplementary Table S8"):
    m = re.search(r"Experiment description:\s*(.+?)(?:\n\s*\n|\n\s*SRA)", chunk, re.S)
    if m:
        exp = re.sub(r"\s+", " ", m.group(1)).strip().rstrip(".")
    for line in chunk.split("\n"):
        mm = re.match(r"\s*([SED]RR\d+)\s+(.+?)\s+([\d,]+)\s*[a-z]?\s*$", line)
        if mm:
            runs[mm.group(1)] = (re.sub(r"\s{2,}", " | ", mm.group(2).strip()), exp)
print(f"Table S8: {len(runs)} runs")

api = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA650228&result=read_run"
       "&fields=run_accession,sample_title,fastq_ftp,fastq_bytes,read_count")
# curl, not urllib: the macOS python this is usually run from has no usable CA bundle
report = subprocess.run(["curl", "-sS", "--retry", "3", api], capture_output=True, text=True, check=True).stdout
ena = {}
for i, line in enumerate(report.split("\n")):
    f = line.rstrip("\r").split("\t")
    if i and len(f) >= 5:
        ena[f[0]] = f

# Only the two phase-I comparisons: the same mock material put through every protocol, which is the whole
# point of this dataset. The phase-II runs vary operator, lot and participating lab instead.
KEEP = {"phase I, comparison of protocols for DNA extraction": "extraction",
        "phase I, comparison of protocols for sequencing library construction": "library"}
out = [["run", "experiment", "protocol", "material", "detail", "read_pairs", "url", "bytes"]]
for r, (meta, exp) in sorted(runs.items()):
    if exp not in KEEP or r not in ena:
        continue
    # pdftotext occasionally splits a two-character code across the column edge ("C protocol L" = CL)
    meta = re.sub(r"\b([A-Z]) protocol ([A-Z0-9])\b", r"protocol \1\2", meta)
    p = re.search(r"protocol ([A-Z][A-Z0-9]?)\b", meta)
    title = ena[r][1]
    material = ("dna_mock" if "DNA mock" in title else
                "cell_mock" if "cell mock" in title else "feces")
    url = next((f for f in ena[r][2].split(";") if f.endswith("_1.fastq.gz")), "")
    size = next((b for f, b in zip(ena[r][2].split(";"), ena[r][3].split(";"))
                 if f.endswith("_1.fastq.gz")), "0")
    if not url:
        continue
    out.append([r, KEEP[exp], p.group(1) if p else "NA", material, meta,
                ena[r][4], "https://" + url, size])
io.open(f"{outdir}/jmbc_runs.tsv", "w", encoding="utf-8").write(
    "\n".join("\t".join(map(str, r)) for r in out) + "\n")
print(f"kept {len(out) - 1} phase-I runs, "
      f"{sum(int(r[7]) for r in out[1:]) / 1e9:.1f} GB of forward reads at full depth")
