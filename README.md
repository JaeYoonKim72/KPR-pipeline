# Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus

Analysis pipelines used in the Korean Pangenome Reference (KPR) study, covering
long-read genome assembly, pangenome graph construction, pangenome-based
short-read mapping and variant calling, and the analysis of the amylase (AMY)
locus.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Pipelines

Run in the order below; each script consumes the output of the previous one.

| # | Script | Description | Documentation |
|---|---|---|---|
| 1 | `genome_assembly_and_statistics.sh` | Long-read and Hi-C read QC, haplotype-resolved *de novo* assembly, polishing, assembly statistics, repeat and gene annotation, structural variant detection | [docs](genome_assembly_and_statistics.md) |
| 2 | `pangenome_construction_and_statistics.sh` | Pangenome graph construction from the quality-controlled assemblies, graph statistics, graph-derived variants, growth curves, bubbles and non-reference sequence | [docs](pangenome_construction_and_statistics.md) |
| 2-1 | `pangenome_variant_comparison.sh` | SV and small-variant decomposition of the graph-derived variants and comparison between the KPR and HPRC haplotypes at the locus and allele level | [docs](pangenome_variant_comparison.md) |
| 3 | `pangenome_short_read_mapping.sh` | Short-read mapping to the pangenome graph and graph-based small variant and structural variant calling | [docs](pangenome_short_read_mapping.md) |
| 4 | `AMY_locus_analysis.sh` | Structural haplotype analysis of the AMY locus (no copy-based SNP population statistics) | [docs](AMY_locus_analysis.md) |
| 5 | `5.Population_genetics_in_AMY_flanking_locus.sh` | Within- and between-population statistics for distal/proximal AMY flanks (focal estimates; optional selscan) | [docs](5.Population_genetics_in_AMY_flanking_locus.md) |

Helper scripts called by the pipelines: `decon_split.sh`, `typecount.py`,
`wline_sample_nodes.py`, `amy_haplotype_tags.py`, `msa_to_vcf.py`,
`hudson_fst.py`, `plot_haplotypes.R`, `sv_merge_alleles.py`,
`dedup_sample.py`, `count_groups.py` and `summarize_variant_comparison.py`.
The new step 5 uses `amy_flank_popgen.py` and `amy_flank_population_stats.py`.
They live next to the pipelines and can also be run on their own; see the
documentation of each pipeline.

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
bash pangenome_variant_comparison.sh all
bash pangenome_short_read_mapping.sh all
bash AMY_locus_analysis.sh all

# Separately, after preparing a joint phased SNP VCF and population lists:
# export VCF=/path/to/joint.phased.biallelic.snps.vcf.gz
# export POPS="KOR=/path/to/KOR.samples CHN=/path/to/CHN.samples"
# export CALLABLE_BED=/path/to/GRCh38.StrictMask.bed.gz
# bash 5.Population_genetics_in_AMY_flanking_locus.sh all
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
(https://ega-archive.org/datasets/EGAD50000002977), and in the Korea BioData
Station (K-BDS) under accession KAP241783 (https://kbds.re.kr/).

The pangenome graph files, together with further detailed code, are available
at https://kgr.appex.kr/.

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Contact

Jaeyoon Kim, Korea Research Institute of Bioscience and Biotechnology (KRIBB) —
jaeyoonkim@kribb.re.kr
