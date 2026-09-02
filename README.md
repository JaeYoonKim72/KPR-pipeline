# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

Analysis pipelines used in the Korean Pangenome Reference (KPR) study, covering
long-read genome assembly, pangenome graph construction, pangenome-based
short-read mapping and variant calling, and the analysis of the amylase (AMY)
locus.

<!-- Add the Zenodo badge here after the first release is archived:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

Run in the order below; each script consumes the output of the previous one.

| # | Script | Description | Documentation |
|---|---|---|---|
| 1 | `genome_assembly_and_statistics.sh` | Long-read and Hi-C read QC, haplotype-resolved *de novo* assembly, polishing, assembly statistics, repeat and gene annotation, structural variant detection | [docs](docs/genome_assembly_and_statistics.md) |
| 2 | `pangenome_construction_and_statistics.sh` | Pangenome graph construction from the quality-controlled assemblies and graph statistics | [docs](docs/pangenome_construction_and_statistics.md) |
| 3 | `pangenome_short_read_mapping.sh` | Short-read mapping to the pangenome graph and graph-based small variant and structural variant calling | [docs](docs/pangenome_short_read_mapping.md) |
| 4 | `AMY_locus_analysis.sh` | Copy number and haplotype structure analysis at the amylase locus | [docs](docs/AMY_locus_analysis.md) |

## Quick start

Every script reads its paths from the environment; nothing is hard-coded.

```bash
export DATA_DIR=/path/to/fastq        # input sequencing data
export REF_DIR=/path/to/reference     # reference files
export OUT_DIR=/path/to/results       # output directory
export LIST_DIR=/path/to/lists        # sample lists
export THREADS=32

bash genome_assembly_and_statistics.sh all
bash pangenome_construction_and_statistics.sh all
bash pangenome_short_read_mapping.sh all
bash AMY_locus_analysis.sh all
```

Each script also accepts a single step name instead of `all`, which is
recommended on a cluster where individual steps are submitted as array jobs.
Run a script without arguments to list its steps. See the documentation linked
above for the inputs, reference files and outputs of each pipeline.

## Requirements

Unix shell environment with the third-party tools listed in the documentation
of each pipeline available on `PATH`. Reference genomes and annotation files
are downloaded separately; the download URLs are given in the configuration
block at the top of each script.

## Data availability

Raw sequencing data are deposited in the European Genome-phenome Archive (EGA)
under study accession EGAS50000002074 and dataset accession EGAD50000002977
(https://ega-archive.org/datasets/EGAD50000002977). Genome assemblies and the
pangenome graph are available at https://kgr.appex.kr/.

## Citation

If you use this code, please cite the paper and the archived software release.
Machine-readable metadata is provided in `CITATION.cff`.

## License

Released under the MIT License. See [LICENSE](LICENSE).
