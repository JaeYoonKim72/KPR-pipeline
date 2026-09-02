# genome_assembly_and_statistics.sh

Haplotype-resolved *de novo* genome assembly, assembly quality assessment,
repeat and gene annotation, and structural variant detection from long-read
sequencing data.

Part of the [KPR-pipeline](../README.md) repository.

## Quick start

```bash
export DATA_DIR=/path/to/fastq        # required
export REF_DIR=/path/to/reference     # required
export OUT_DIR=/path/to/results       # default: ./results
export LIST_DIR=/path/to/lists        # default: ./lists
export THREADS=32

bash genome_assembly_and_statistics.sh all
```

Steps can also be run one at a time, which is recommended on a cluster:

```bash
bash genome_assembly_and_statistics.sh read_qc
bash genome_assembly_and_statistics.sh assembly
```

## Steps

| Step | Description |
|---|---|
| `read_qc` | FastQC, HiFiAdapterFilt, Porechop, Trimmomatic, raw read metrics |
| `assembly` | HiFiasm, HiFiasm Hi-C mode, Flye, contaminant removal |
| `polishing` | Inspector evaluation and correction |
| `asm_stats` | Inspector summary, QUAST-LG, NGx curves, Minigraph misjoins, HMM-Flagger |
| `annotation` | RepeatMasker, etrf, sdust, Liftoff |
| `sv_calling` | pbmm2, pbsv, Sniffles, SVIM, SVIM-asm, SURVIVOR merging |
| `all` | All of the above, in order |

## Input

FASTQ files in `DATA_DIR`, named by sample ID:

```
SAMPLE.hifi.fastq.gz                            PacBio HiFi
SAMPLE.ont.fastq.gz                             ONT PromethION
SAMPLE.hic_1.fastq.gz  SAMPLE.hic_2.fastq.gz    Hi-C
```

Sample lists in `LIST_DIR`, one ID per line: `hifi.list`, `ont.list`,
`hic.list`. A sample in both `hifi.list` and `hic.list` is assembled in Hi-C
phasing mode. `qc_pass.list` is written by `asm_stats` and lists the assemblies
meeting N50 >= 20 Mb and <= 1,500 contigs.

## Configuration

| Variable | Default | Description |
|---|---|---|
| `DATA_DIR` | required | input FASTQ directory |
| `REF_DIR` | required | reference file directory |
| `OUT_DIR` | `./results` | output directory |
| `LIST_DIR` | `./lists` | sample list directory |
| `THREADS` | `32` | threads per job |
| `GENOME_SIZE` | `3200000000` | genome size used for raw depth |
| `MIN_N50` | `20000000` | minimum contig N50 |
| `MAX_CONTIG` | `1500` | maximum number of contigs |

Individual reference files can be overridden as well, for example
`GRCH38=/path/to/GRCh38.fa`.

## Reference data

Download URLs are given in the configuration block at the top of the script.
Required files: GRCh38 primary assembly, Ensembl Regulatory Build GFF,
T2T-CHM13 v2.0, GENCODE v38, an augmented CHM13 satellite library, and a
Kraken2 database.

## Output

```
results/
├── read_qc/       FastQC reports, trimmed reads, read metrics
├── assembly/      raw and contaminant-filtered assemblies
├── polishing/     polished assemblies, Inspector summaries
├── statistics/    QUAST-LG reports, NGx curves, misjoins, HMM-Flagger labels
├── annotation/    repeat, VNTR, low-complexity and gene annotations
├── sv/            per-sample and shared structural variant call sets
└── logs/
```

Main tables:

- `statistics/inspector/inspector_summary.tsv` — quality value, self-mapping
  rate, self-mapping depth and small-scale error counts per haplotype
- `statistics/assembly_qc_summary.tsv` — contig number, N50 and PASS/FAIL
  status against the contiguity criteria
- `sv/merged/shared_SV.vcf` — structural variants of at least 51 bp detected in
  every sample retained after assembly quality control

The polished assemblies in `polishing/corrected/` are the input for
`pangenome_construction_and_statistics.sh`.

## Software

FastQC 0.12.1 · HiFiAdapterFilt 3.0 · Porechop 0.2.4 · Trimmomatic 0.39 ·
HiFiasm 0.24 (0.25 for ONT mode) · Flye 2.9.5 · gfatools · Kraken2 · seqkit ·
Inspector 1.3 · QUAST-LG 5.3 · Minigraph 0.21 · HMM-Flagger 1.0.0 ·
RepeatMasker 4.1.6 · etrf · sdust · Liftoff 1.6.3 · pbmm2 1.17 · pbsv 2.9.0 ·
Sniffles 2.2 · SVIM 2.0.0 · SVIM-asm 1.0.3 · SURVIVOR · samtools · minimap2 ·
jq

All executables must be available on `PATH`.
