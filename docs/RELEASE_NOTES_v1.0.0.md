# KPR-pipeline v1.0.0

This release archives the source code and methodological documentation for the Korean Pangenome Reference (KPR) workflows accompanying the study *Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus*.

## Included workflows

1. Long-read genome assembly, quality assessment, and structural variant discovery
2. Ensemble gene annotation and copy-number analysis
3. Pangenome graph construction and characterization
4. Pangenome variant-count comparisons
5. Small-variant and SV benchmarking with GIAB trio genotype inheritance inconsistency analysis
6. Graph and GRCh38 short-read mapping and variant calling
7. AMY locus structural haplotype analysis
8. Population-genetic analysis of the AMY distal and proximal flanking regions

The executable entry points are numbered `1.*.sh` through `8.*.sh` in the repository root, with supporting methods under `src/`.

## Reproducibility and licensing

The repository is licensed under MIT. Sequencing data, reference genomes, confidence BEDs, truth sets, pangenome graph files and third-party executables are not bundled. Refer to the stage-specific guides for requirements and limitations. Passing synthetic smoke tests does not establish reproduction of all manuscript-scale results. A machine-readable citation is in `CITATION.cff`.
