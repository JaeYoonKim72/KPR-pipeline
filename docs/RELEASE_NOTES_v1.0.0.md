# KPR-pipeline v1.0.0

Source code and analysis procedures for the Korean Pangenome Reference (KPR) workflows accompanying the study *Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus*.

## Included workflows

1. Genome assembly, quality assessment, repeat annotation and structural variant detection
2. Ensemble gene annotation and copy-number analysis
3. Pangenome graph construction and characterization
4. Locus- and allele-level pangenome variant comparisons
5. Variant benchmarking and Mendelian genotype inheritance inconsistency
6. Short-read mapping and variant calling
7. AMY structural haplotype analysis
8. AMY flanking-region population genetics

Each workflow includes an executable Shell script and documentation describing its inputs, analysis parameters, commands and outputs. Additional study analyses are described in the corresponding methods documents.

## Validation and data

Validation includes Shell and Python syntax checks and synthetic regression tests. Sequencing data, reference genomes, truth sets, confidence BEDs, pangenome graph files and third-party software are obtained separately. See the [README](../README.md) and pipeline documentation.

## Citation and license

See [CITATION.cff](../CITATION.cff) and [LICENSE](../LICENSE) (MIT).
