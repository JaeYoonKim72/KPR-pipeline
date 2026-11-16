# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

This repository contains the analysis workflows, computational scripts, and
methodological documentation accompanying the study *Constructing a Korean
pangenome reveals distinct haplotypes at the amylase locus*. The study aims
to improve the representation of Korean genomic diversity by integrating
high-quality long-read genome assemblies into graph-based pangenome
references. Long-read sequencing datasets from 22 Korean individuals
(19 newly sequenced and three from a previous study) were processed through
read quality control, de novo assembly, polishing, assembly quality
assessment, repeat and gene annotation, and structural variant discovery.
After quality filtering, 28 haploid assemblies from 14 individuals were
used to construct the Korean Pangenome Reference (KPR) and an augmented
HPRC v1.0+KPR pangenome graph.

The workflows describe pangenome graph construction and indexing,
graph-derived variant identification, population-panel comparisons,
pangenome growth, and the characterization of non-reference sequences.
They also cover short-read alignment and variant calling using
pangenome graphs and the linear GRCh38 reference for comparative analyses.
At the structurally complex amylase (AMY) locus, the study further examines
gene copy number, gene order, and structural haplotype diversity, followed
by locus-level population-genetic analyses of the unique distal and
proximal flanking regions. The accompanying documentation specifies
inputs, software dependencies, analytical procedures, and outputs to
support reproducible use and interpretation of the implemented workflows.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

The primary analysis is organized into eight numbered stages. Each stage
has an executable Shell entry point in the repository root and its
method documentation under `src/`. Stages 2 and 5 are explicit
placeholders and do not implement analyses.

| Stage | Main script | Purpose | Documentation | Version |
|:--:|---|---|---|---|
| 1 | `1.genome_assembly_and_statistics.sh` | Long-read assembly, quality assessment, annotation and structural variants | [Documentation](src/genome_assembly_and_statistics.md) | v1.0 |
| 2 | `2.Ensemble_gene_annotation.sh` | Reserved Ensembl/GENCODE annotation workflow | [Documentation](src/Ensemble_gene_annotation.md) | v1.0 |
| 3 | `3.pangenome_construction_and_statistics.sh` | KPR and HPRC+KPR graph construction and statistics | [Documentation](src/pangenome_construction_and_statistics.md) | v1.0 |
| 4 | `4.pangenome_variant_comparison.sh` | KPR/HPRC graph variant decomposition and allele comparison | [Documentation](src/pangenome_variant_comparison.md) | v1.0 |
| 5 | `5.Variant_benchmark.sh` | Reserved truth-set benchmarking workflow | [Documentation](src/Variant_benchmark.md) | v1.0 |
| 6 | `6.pangenome_short_read_mapping.sh` | Executable short-read mapping, graph variant calling and GRCh38 baseline | [Documentation](src/pangenome_short_read_mapping.md) | v1.0 |
| 7 | `7.AMY_locus_analysis_in_pagenomegraph.sh` | AMY locus structural haplotypes and phylogeny | [Documentation](src/AMY_locus_analysis_in_pagenomegraph.md) | v1.0 |
| 8 | `8.Population_genetics_in_AMY_flanking_locus.sh` | Distal and proximal AMY flanking-region population genetics | [Documentation](src/Population_genetics_in_AMY_flanking_locus.md) | v1.0 |

**Helper scripts called by the pipelines:** `decon_split.sh`,
`typecount.py`, `wline_sample_nodes.py`, `amy_haplotype_tags.py`,
`msa_to_vcf.py`, `hudson_fst.py`, `plot_haplotypes.R`,
`sv_merge_alleles.py`, `dedup_sample.py`, `count_groups.py`,
and `summarize_variant_comparison.py`. These tools are in `src/`;
they are invoked by the main pipelines where applicable and can
also be run independently using the documented arguments.

The [pangenome short-read mapping workflow](src/pangenome_short_read_mapping.md)
is provided as a dedicated executable stage because graph alignment,
surjection, variant calling and linear-reference comparison can be
performed separately from de novo graph construction. It is run with
`6.pangenome_short_read_mapping.sh` (stage 6), not merely described as
a workflow without executable code.

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

bash 1.genome_assembly_and_statistics.sh all
bash 3.pangenome_construction_and_statistics.sh all
bash 4.pangenome_variant_comparison.sh all
bash 6.pangenome_short_read_mapping.sh all
bash 7.AMY_locus_analysis_in_pagenomegraph.sh all
```

The AMY population-genetics stage takes **a joint, phased, biallelic
SNP VCF and disjoint population sample-ID lists**, not copy-specific
gene alignments. It can be run independently of graph structural steps:

```bash
mamba create -n amy_popgen -c conda-forge -c bioconda python bcftools bedtools htslib selscan=3.0.1
conda activate amy_popgen

export VCF=/path/to/joint.phased.biallelic.snps.vcf.gz
export POPS="KOR=/path/to/KOR.txt CHN=/path/to/CHN.txt"
export CALLABLE_BED=/path/to/GRCh38.StrictMask.nonN.bed.gz
export OUT_DIR=/path/to/results/amy_flanking
bash 8.Population_genetics_in_AMY_flanking_locus.sh all
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
bash 8.Population_genetics_in_AMY_flanking_locus.sh all
```

Our manuscript instead used **10,000** bootstrap replicates,
full leave-one-person-out jackknife, and **chr1–22 / 5-Mb
physical chromosome-stratified block bootstrap (10,000 draws)**.
Those settings, estimator equations, window definitions and currently
unsupported metrics are documented in [stage 8](src/Population_genetics_in_AMY_flanking_locus.md).
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
- Stage 8 needs GRCh38 flanking SNP calls and a
  `StrictMask(20160622) ∩ non-N` callable BED for appropriate
  denominator lengths. It does not download private cohort data.
- Stage 8's optional genomewide feature needs an indexed
  joint all-chromosome VCF (chr1 for the default pilot, chr1–22
  for the manuscript-scale configuration).
- Stages 2 and 5 are placeholders; they exit with
  `NOT IMPLEMENTED` rather than report false success.

## Data availability

**Newly generated Korean long-read and Hi-C sequencing data.** Long-read
sequencing (LRS) data were generated for **19 Korean individuals**:
**12 individuals** were sequenced using PacBio Revio (HiFi), and
**seven individuals** using Oxford Nanopore PromethION.
Hi-C sequencing was additionally performed for **two of the
Revio-sequenced individuals**, one of whom was also sequenced
using PromethION. These raw sequencing data are deposited in
the European Genome-phenome Archive (EGA) under study
[EGAS50000002074](https://ega-archive.org/studies/EGAS50000002074)
and dataset
[EGAD50000002977](https://ega-archive.org/datasets/EGAD50000002977).

**Previously generated Korean PacBio Sequel II data.** Long-read
sequencing data for **three additional Korean individuals**,
generated using the **PacBio Sequel II** platform in our previous
study, are available under EGA study
[EGAS50000000375](https://ega-archive.org/studies/EGAS50000000375)
and dataset
[EGAD50000000552](https://ega-archive.org/datasets/EGAD50000000552).
Together, the newly generated and previously available long-read
datasets represent **22 Korean individuals**.

**Korean BioData Station (K-BDS).** The **same newly generated
sequencing data from 19 Korean individuals** (12 PacBio Revio
and seven Oxford Nanopore PromethION datasets, with additional
Hi-C sequencing for two individuals) are also deposited at
[K-BDS](https://kbds.re.kr/) under accession **KAP241783**.
The K-BDS deposit also includes short-read sequencing data
used in this study.

**Data access committee.** Access to controlled EGA data is
governed by the KRIBB Genomic Medicine Research Center
Data Access Committee (DAC),
[EGAC50000001142](https://ega-archive.org/dacs/EGAC50000001142).

**Pangenome references and code.** GFA graph files and associated
resources are available through the
[Korean Genome Research project](https://kgr.appex.kr/).
Executable workflow scripts and supporting analysis methods
are maintained in this repository.

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jae-Yoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
