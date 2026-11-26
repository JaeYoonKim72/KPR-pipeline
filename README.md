# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

Version 1.0.0

This repository provides computational pipelines and analysis procedures for the construction and characterization of the Korean Pangenome Reference (KPR). We assembled genomes from long-read sequencing data for 22 Korean individuals and selected 28 haplotype assemblies from 14 individuals for pangenome construction.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

| Stage | Script | Analysis | Documentation | Supplementary tables |
|:--:|---|---|---|---|
| 1 | `1.genome_assembly_and_statistics.sh` | Sample selection, genome assembly, quality assessment, repeat annotation and structural variant detection | [Methods](src/genome_assembly_and_statistics.md) | S1–S12, S18 |
| 2 | `2.Ensemble_gene_annotation.sh` | Ensemble gene annotation and copy-number gain analysis | [Methods](src/Ensemble_gene_annotation.md) | S13–S17 |
| 3 | `3.pangenome_construction_and_statistics.sh` | Pangenome construction, graph statistics, growth curves and non-reference sequence analysis | [Methods](src/pangenome_construction_and_statistics.md) | S19–S20, S22–S23 |
| 4 | `4.compare_pangenome_variant_counts.sh` | Locus- and allele-level variant comparisons | [Methods](src/compare_pangenome_variant_counts.md) | S21 |
| 5 | `5.Variant_benchmark.sh` | Variant benchmarking and Mendelian genotype inheritance inconsistency | [Methods](src/Variant_benchmark.md) | S24–S26 |
| 6 | `6.pangenome_short_read_mapping.sh` | Short-read mapping and variant calling | [Methods](src/pangenome_short_read_mapping.md) | S24–S27 |
| 7 | `7.AMY_locus_analysis_in_pagenomegraph.sh` | AMY structural haplotype analysis | [Methods](src/AMY_locus_analysis_in_pagenomegraph.md) | S28–S29 |
| 8 | `8.Population_genetics_in_AMY_flanking_locus.sh` | AMY flanking-region population genetics | [Methods](src/Population_genetics_in_AMY_flanking_locus.md) | S30–S32 |

## Usage

Each pipeline can be executed independently with dataset-specific inputs and analysis parameters. Software requirements, input preparation, execution commands, output files and additional analysis procedures are provided in the corresponding documentation.

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

**Pangenome references and code.** Pangenome graph files and associated
resources are available through the
[Korean Genome Research project](https://kgr.appex.kr/).
Executable workflow scripts and supporting analysis methods
are maintained in this repository.

## Citation

See [CITATION.cff](CITATION.cff).

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jae-Yoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
