# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

This repository provides computational workflows and analysis procedures for the Korean Pangenome Reference (KPR). We generated or incorporated long-read sequencing data from 22 Korean individuals and retained 28 haplotype assemblies from 14 individuals after assembly quality assessment. We constructed KPR and the augmented HPRC+KPR pangenome, compared graph variants at locus and allele levels, evaluated variant calling, and characterized structural haplotypes and flanking-region population genetics at the AMY locus.

The scripts accept user-supplied references, assemblies, reads, variant calls and population manifests. The accompanying methods describe the study-specific configurations and additional analyses performed using the pipeline outputs and external analysis tools.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

| Stage | Script | Analysis | Documentation | Study tables |
|:--:|---|---|---|---|
| 1 | `1.genome_assembly_and_statistics.sh` | Sample selection, assembly, quality assessment, repeats and consensus SVs | [Methods](src/genome_assembly_and_statistics.md) | S1–S12, S18 |
| 2 | `2.Ensemble_gene_annotation.sh` | Ensemble annotation and copy-number gains | [Methods](src/Ensemble_gene_annotation.md) | S13–S17 |
| 3 | `3.pangenome_construction_and_statistics.sh` | Graph construction, graph statistics, growth and non-reference sequences | [Methods](src/pangenome_construction_and_statistics.md) | S19–S20, S22–S23 |
| 4 | `4.compare_pangenome_variant_counts.sh` | Locus-level and allele-level variant comparisons | [Methods](src/compare_pangenome_variant_counts.md) | S21 |
| 5 | `5.Variant_benchmark.sh` | GIAB and KOR benchmarking and trio inheritance inconsistency | [Methods](src/Variant_benchmark.md) | S24–S26 |
| 6 | `6.pangenome_short_read_mapping.sh` | Short-read mapping and variant calling | [Methods](src/pangenome_short_read_mapping.md) | S24–S27 |
| 7 | `7.AMY_locus_analysis_in_pagenomegraph.sh` | AMY structural haplotypes and phylogeny | [Methods](src/AMY_locus_analysis_in_pagenomegraph.md) | S28–S29 |
| 8 | `8.Population_genetics_in_AMY_flanking_locus.sh` | AMY flanking-region population genetics | [Methods](src/Population_genetics_in_AMY_flanking_locus.md) | S30–S32 |

## Requirements

Bash, Python 3.9+ and the analysis software listed in each methods document. Variant benchmarking with Truvari v5.5.0 requires Python ≥3.12. Install the specified versions and prepare reference indexes and analysis-specific input files before running each stage.

## Run

Each stage has a Shell entry point and a methods document with setup, execution, step descriptions, configuration, outputs, study analysis procedures, references and software versions. Stages can be executed independently using their documented inputs.

~~~bash
export DATA_DIR=/path/to/reads
export REF_DIR=/path/to/reference
export OUT_DIR=/path/to/results
export LIST_DIR=/path/to/lists
export THREADS=32

bash 1.genome_assembly_and_statistics.sh all
bash 2.Ensemble_gene_annotation.sh all
bash 3.pangenome_construction_and_statistics.sh all
bash 4.compare_pangenome_variant_counts.sh all
bash 6.pangenome_short_read_mapping.sh all
bash 5.Variant_benchmark.sh all
bash 7.AMY_locus_analysis_in_pagenomegraph.sh all
bash 8.Population_genetics_in_AMY_flanking_locus.sh all
~~~

Configure each stage's additional required variables according to its methods document. In particular, supply the stage 1 acrocentric BED for misjoin filtering, the stage 2 assembly/annotation manifest, stage 3 ordered assembly seqfiles, stage 4 population map, stage 5 truth and query sets with confidence BEDs, stage 6 matching graph indexes, stage 7 structural-label map and stage 8 SNP inputs, callable mask and population lists.

The study used 302 Korean individuals and eight 1000 Genomes populations for AMY flanking analyses. Diversity, LD, Hudson's FST and phased XP-nSL were assessed in the distal and proximal flanks; the structurally variable region was excluded. See stage 8 for the study's 10,000-replicate individual and 5-Mb genomic block resampling configurations, separate LD background calculations and phased XP-nSL normalization and subsampling.

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

## Citation

See [CITATION.cff](CITATION.cff).

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jae-Yoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
