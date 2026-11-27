# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

Version 1.0.0

This repository provides computational pipelines and analysis procedures for the construction and characterization of the Korean Pangenome Reference (KPR). We assembled genomes from long-read sequencing data for 22 Korean individuals and selected 28 haplotype assemblies from 14 individuals for pangenome construction.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

| Stage | Script | Analysis | Documentation |
|:--:|---|---|---|
| 1 | `1.genome_assembly_and_statistics.sh` | Sample selection, genome assembly, quality assessment, repeat annotation and structural variant detection | [Methods](src/genome_assembly_and_statistics.md) |
| 2 | `2.Ensemble_gene_annotation.sh` | Ensemble gene annotation and copy-number gain analysis | [Methods](src/Ensemble_gene_annotation.md) |
| 3 | `3.pangenome_construction_and_statistics.sh` | Pangenome construction, graph statistics, growth curves and non-reference sequence analysis | [Methods](src/pangenome_construction_and_statistics.md) |
| 4 | `4.compare_pangenome_variant_counts.sh` | Locus- and allele-level variant comparisons | [Methods](src/compare_pangenome_variant_counts.md) |
| 5 | `5.Variant_benchmark.sh` | Variant benchmarking and Mendelian genotype inheritance inconsistency | [Methods](src/Variant_benchmark.md) |
| 6 | `6.pangenome_short_read_mapping.sh` | Short-read mapping and variant calling | [Methods](src/pangenome_short_read_mapping.md) |
| 7 | `7.AMY_locus_analysis_in_pagenomegraph.sh` | AMY structural haplotype analysis | [Methods](src/AMY_locus_analysis_in_pagenomegraph.md) |
| 8 | `8.Population_genetics_in_AMY_flanking_locus.sh` | AMY flanking-region population genetics | [Methods](src/Population_genetics_in_AMY_flanking_locus.md) |

## Usage

Each pipeline can be executed independently with dataset-specific inputs and analysis parameters. Software requirements, input preparation, execution commands, output files and additional analysis procedures are provided in the corresponding documentation.

## Data availability

**Newly generated Korean long-read and Hi-C sequencing data.** Long-read sequencing (LRS) data were generated for **19 Korean individuals**: **12 individuals** were sequenced using PacBio Revio (HiFi), and **seven individuals** using Oxford Nanopore PromethION. Hi-C sequencing was additionally performed for **two of the Revio-sequenced individuals**, one of whom was also sequenced using PromethION. These raw sequencing data are deposited in the European Genome-phenome Archive (EGA) under study [EGAS50000002074](https://ega-archive.org/studies/EGAS50000002074) and dataset [EGAD50000002977](https://ega-archive.org/datasets/EGAD50000002977).

**Previously generated Korean PacBio Sequel II data.** Long-read sequencing data for **three additional Korean individuals**, generated using the **PacBio Sequel II** platform in our previous study, are available under EGA study [EGAS50000000375](https://ega-archive.org/studies/EGAS50000000375) and dataset [EGAD50000000552](https://ega-archive.org/datasets/EGAD50000000552). Together, the newly generated and previously available long-read datasets represent **22 Korean individuals**.

**Genome assemblies.** All 37 polished genome assembly files (FASTA) from **22 Korean individuals** — including the three sequenced on PacBio Sequel II in the previous study — are deposited together with the newly generated raw sequencing data under EGA study [EGAS50000002074](https://ega-archive.org/studies/EGAS50000002074) and dataset [EGAD50000002977](https://ega-archive.org/datasets/EGAD50000002977). Of these, 28 high-quality haploid assemblies (`KOR-001` to `KOR-014`, two haplotypes each) were used to build the KPR.

**Korean BioData Station (K-BDS).** The **same newly generated sequencing data from 19 Korean individuals** (12 PacBio Revio and seven Oxford Nanopore PromethION datasets, with additional Hi-C sequencing for two individuals) are also deposited at [K-BDS](https://kbds.re.kr/) under accession **KAP241783**. The K-BDS deposit also includes short-read sequencing data used in this study. Genome assemblies are not held at K-BDS.

**Pangenome graph files.** The KPR and HPRC+KPR graphs in GFA format are available under controlled access from the [Korean Genome Research data download page](https://kgr.appex.kr/resource/data-download/). Because the graphs contain individual-level haplotype sequences, they are governed by the same access terms as the KOR genome assemblies. Researchers with current DAC approval covering the KOR genome assemblies can obtain the GFA files by submitting a short form with their name, institutional affiliation, approved EGA account email address and requested graph file. **A download link is sent to the approved email address within two working days following verification** and is valid for up to seven days, or until DAC approval expires or is revoked, whichever occurs earlier. This is a verification of existing access permissions, **not a second research review**. Downloading the FASTA files is not a prerequisite. The graphs can also be reconstructed from the deposited assemblies and corresponding reference resources using stage 3 of this repository.

**Data access committee.** Access to data under EGAS50000002074 and the KPR and HPRC+KPR graph files is governed by the KRIBB Genomic Medicine Research Center Data Access Committee (DAC), [EGAC50000001142](https://ega-archive.org/dacs/EGAC50000001142). The access policy — the categories of research that will be approved, the uses that will not, the review procedure, and the expected decision time — is published at [KRIBB Korean Pangenome Project Data Access Policy](https://github.com/JaeYoonKim72/Policy_KPRv1.0). Applications are assessed only against the published eligibility criteria; they are not evaluated on scientific merit, novelty, research expertise or overlap with the data producers' own work. The earlier PacBio Sequel II raw-read dataset EGAD50000000552 remains subject to its own EGA access policy and DAC. Access to K-BDS deposits follows the conditions applicable to those deposits.

**Openly available resources.** The analysis pipelines and methods documentation in this repository (MIT License) and the supplementary tables published with the primary study are available without application to the DAC.

## Citation

See [CITATION.cff](CITATION.cff).

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jae-Yoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
