# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

Analysis workflows and reproducible methods for the Korean Pangenome
Reference (KPR) study: genome assembly, annotations, pangenome
construction, graph-variant comparison, variant benchmarking, AMY
structural haplotypes and flanking-region population genetics.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

The primary analysis is organized into eight numbered stages. Each stage
has an executable Shell entry point in the repository root and its
method documentation under `src/`. Stages 2 and 5 are explicit
placeholders and do not implement analyses.

| Stage | Main script | Purpose | Documentation | Status |
|:--:|---|---|---|---|
| 1 | `1.genome_assembly_and_statistics.sh` | Long-read assembly, quality assessment, annotation and structural variants | [Documentation](src/genome_assembly_and_statistics.md) | Implemented |
| 2 | `2.Ensemble_gene_annotation.sh` | Reserved Ensembl/GENCODE annotation workflow | [Documentation](src/Ensemble_gene_annotation.md) | Placeholder |
| 3 | `3.pangenome_construction_and_statistics.sh` | KPR and HPRC+KPR graph construction and statistics | [Documentation](src/pangenome_construction_and_statistics.md) | Implemented |
| 4 | `4.pangenome_variant_comparison.sh` | KPR/HPRC graph variant decomposition and allele comparison | [Documentation](src/pangenome_variant_comparison.md) | Implemented |
| 5 | `5.Variant_benchmark.sh` | Reserved truth-set benchmarking workflow | [Documentation](src/Variant_benchmark.md) | Placeholder |
| 6 | `6.pangenome_short_read_mapping.sh` | Executable short-read mapping, graph variant calling and GRCh38 baseline | [Documentation](src/pangenome_short_read_mapping.md) | Implemented |
| 7 | `7.AMY_locus_analysis_in_pagenomegraph.sh` | AMY locus structural haplotypes and phylogeny | [Documentation](src/AMY_locus_analysis_in_pagenomegraph.md) | Implemented |
| 8 | `8.Population_genetics_in_AMY_flanking_locus.sh` | Distal and proximal AMY flanking-region population genetics | [Documentation](src/Population_genetics_in_AMY_flanking_locus.md) | SNP summaries; additional normalization needed |

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

**Newly generated Korean sequencing data.** The EGA study
[EGAS50000002074](https://ega-archive.org/studies/EGAS50000002074)
and dataset
[EGAD50000002977](https://ega-archive.org/datasets/EGAD50000002977)
contain the newly generated long-read data described in the
manuscript. The dataset was released by EGA on 19 September 2026;
the metadata are publicly visible and requests for controlled
access can be submitted through EGA.

**Previously generated Korean long-read data.** The three additional
Korean individuals are registered under the released EGA study
[EGAS50000000375](https://ega-archive.org/studies/EGAS50000000375).
The long-read data used in this study are in dataset
[EGAD50000000552](https://ega-archive.org/datasets/EGAD50000000552).

**Data access committee.** The KRIBB Genomic Medicine Research
Center DAC is listed at
[EGAC50000001142](https://ega-archive.org/dacs/EGAC50000001142).
Release of an EGA record does not imply unrestricted download
of participant-level sequencing data; the applicable access
procedure still applies.

**Korean BioData Station (K-BDS).** The deposited project has accession
`KAP241783` at [K-BDS](https://kbds.re.kr/). Reviewers and editors
were given a separate reviewer-access link; the regular
project record is expected to become publicly searchable
upon publication. The reviewer-specific access URL is not
embedded in this public repository.

**Pangenome references and code.** GFA graph files and associated
resources are provided through the
[Korean Genome Research project](https://kgr.appex.kr/).
Executable workflow scripts and supporting analysis
methods are maintained in this repository.

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jaeyoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
