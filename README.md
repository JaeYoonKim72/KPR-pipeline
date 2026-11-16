# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

Analysis workflows and reproducible methods for the Korean Pangenome
Reference (KPR) study: genome assembly, annotations, pangenome
construction, graph-variant comparison, variant benchmarking, AMY
structural haplotypes and flanking-region population genetics.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

The numbered stages below follow the project's planned manuscript workflow.
**Stages 2 and 5 are placeholders** (not implemented; their scripts
explicitly refuse `all` until developed). Stage 7 uses
`Population_genetics_in_AMY_flanking_locus.sh`.

| # | Script | Description | Documentation | Status |
|---|---|---|---|---|
| 1 | `genome_assembly_and_statistics.sh` | Long-read and Hi-C QC; haplotype-resolved assembly, polishing, QC, repeats, SVs | [docs](genome_assembly_and_statistics.md) | Implemented |
| 2 | `Ensemble_gene_annotation.sh` | Reserved Ensembl/GENCODE gene annotation stage | [docs](Ensemble_gene_annotation.md) | Placeholder |
| 3 | `pangenome_construction_and_statistics.sh` | Graph construction, graph metrics, bubble and nonreference sequence analyses | [docs](pangenome_construction_and_statistics.md) | Implemented |
| 4 | `pangenome_variant_comparison.sh` | SV/small-variant deconstruction, allele merging and KPR–HPRC comparison | [docs](pangenome_variant_comparison.md) | Implemented |
| 5 | `Variant_benchmark.sh` | Reserved variant benchmarking against truth sets | [docs](Variant_benchmark.md) | Placeholder |
| 6 | `AMY_locus_analysis_in_pagenomegraph.sh` | Structural AMY haplotype architectures and tree in the HPRC+KPR graph; no copy-specific population SNP statistics | [docs](AMY_locus_analysis_in_pagenomegraph.md) | Implemented |
| 7 | `Population_genetics_in_AMY_flanking_locus.sh` | Distal/proximal locus-level π, θW, Tajima D, LD, nSL, FST, XP-nSL; optional individual and genomic CIs | [docs](Population_genetics_in_AMY_flanking_locus.md) | Focal + SNP-based optional CIs; incomplete normalized haplotype/background metrics |

The existing [pangenome short-read mapping](pangenome_short_read_mapping.md)
workflow (`pangenome_short_read_mapping.sh`) remains available as
a **supporting standalone pipeline** rather than a numbered stage in
this revised seven-step overview.

**Scope of the AMY revision:** The graph-based stage 6 describes
structural variation. The population-genetic stage 7 analyzes
unique distal/proximal flanking sequences at the locus level. It does
not calculate per-copy SNP π/FST/LD for the paralogous AMY genes.

## Quick start

All scripts read paths and options from environment variables (or
the listed command-line arguments). Create/activate the environment
appropriate to each analysis before starting. Example run sequence
for **implemented** stages:

```bash
export DATA_DIR=/path/to/fastq
export REF_DIR=/path/to/reference
export GRAPH_DIR=/path/to/graphs
export OUT_DIR=/path/to/results
export LIST_DIR=/path/to/lists
export THREADS=32

bash genome_assembly_and_statistics.sh all
# Step 2 Ensemble_gene_annotation.sh: planned, not implemented
bash pangenome_construction_and_statistics.sh all
bash pangenome_variant_comparison.sh all
# Step 5 Variant_benchmark.sh: planned, not implemented
bash AMY_locus_analysis_in_pagenomegraph.sh all
```

The AMY population-genetics stage takes **a joint, phased, biallelic
SNP VCF and disjoint population sample-ID lists**, not copy-specific
gene alignments. It can be run independently of graph structural steps:

```bash
# Recommended: activate the environment providing selscan v3.0.1
mamba create -n amy_popgen -c conda-forge -c bioconda python bcftools selscan=3.0.1
conda activate amy_popgen

export VCF=/path/to/joint.phased.biallelic.snps.vcf.gz
export POPS="KOR=/path/to/KOR.txt CHN=/path/to/CHN.txt"
export CALLABLE_BED=/path/to/GRCh38.StrictMask.nonN.bed.gz
export OUT_DIR=/path/to/results/amy_flanking
bash Population_genetics_in_AMY_flanking_locus.sh all
```

With **one** population, only within-population statistics are
performed. With **two or more**, the first population is the query
for Hudson FST and phased XP-nSL comparisons. If selscan v3.0.1 is
missing, nSL and XP-nSL are **skipped with a warning** while other
analyses continue.

Individual bootstrap/jackknife and genomic block bootstrap default
to **OFF**. To enable fast exploratory CIs:

```bash
INDIVIDUAL_CI=1 IND_BOOT_REPS=100 IND_JACKKNIFE_FRACTION=0.30 \
GENOMIC_BOOTSTRAP=1 GENOMIC_CHROMS=chr1 GENOMIC_BLOCK_MB=50 \
GENOMIC_BOOT_REPS=100 GENOMIC_WINDOW_STEP=1000 \
bash Population_genetics_in_AMY_flanking_locus.sh all
```

Our manuscript instead used **10,000** bootstrap replicates,
full leave-one-person-out jackknife, and **chr1–22 / 5-Mb
physical chromosome-stratified block bootstrap (10,000 draws)**.
Those settings, estimator equations, window definitions and currently
unsupported metrics are documented in [stage 7](Population_genetics_in_AMY_flanking_locus.md).
Use verified manuscript tables rather than automatically equating
exploratory public-pipeline values with the final X1/X2/X3 numbers.

## Requirements

- Unix/Bash, Python 3.9+, and the tools noted in each stage's
  documentation. Variant extraction requires `bcftools` and a
  coordinate-sorted/indexed VCF; structural graph workflows require
  graph tools (minigraph, gfabase, odgi, etc.).
- To obtain **nSL and XP-nSL**, activate a mamba/conda environment
  containing **selscan v3.0.1**. Missing selscan skips only those
  statistics, never installs it automatically.
- Stage 7 needs GRCh38 flanking SNP calls and a
  `StrictMask(20160622) ∩ non-N` callable BED for appropriate
  denominator lengths. It does not download private cohort data.
- Stage 7's optional genomewide feature needs an indexed
  joint all-chromosome VCF (chr1 for the default pilot, chr1–22
  for the manuscript-scale configuration).
- Stages 2 and 5 are placeholders; they exit with
  `NOT IMPLEMENTED` rather than report false success.

### Helper scripts

`amy_flank_popgen.py` (focal π/θW/Tajima D/LD/FST),
`amy_flank_population_stats.py` (matched phased SNPs for XP-nSL),
`amy_flank_individual_ci.py` (optional individual CIs),
`amy_flank_genomic_background.py` (optional windowed genomic percentiles),
`amy_haplotype_tags.py` and other helper scripts are kept in the
repository root. Legacy `hudson_fst.py` and `msa_to_vcf.py` are
retained for provenance but **not invoked** for paralog-specific
population genetic statistics by stage 6 or 7.

## Data availability

Raw sequencing data are deposited in the European Genome-phenome Archive (EGA)
under study accession EGAS50000002074 and dataset accession EGAD50000002977
(https://ega-archive.org/datasets/EGAD50000002977), and in the Korea BioData
Station (K-BDS) under accession KAP241783 (https://kbds.re.kr/).

The pangenome graph files, together with further detailed code, are available
at https://kgr.appex.kr/.

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jaeyoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
