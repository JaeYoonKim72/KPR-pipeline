# genome_assembly_and_statistics.sh

Haplotype-resolved *de novo* genome assembly, assembly quality assessment,
repeat and gene annotation, and structural variant detection from long-read
sequencing data.

Part of the [KPR-pipeline](README.md) repository.

---

## 1. Set up

### 1.1 Software

Install the tools listed under [Software](#software) and make sure they are on
`PATH`.

### 1.2 Reference files

Download the reference files into one directory. The URLs are given in the
configuration block at the top of the script.

```bash
mkdir -p reference && cd reference

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip -d chm13v2.0.fa.gz && samtools faidx chm13v2.0.fa

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gzip -d gencode.v38.annotation.gtf.gz

wget https://ftp.ensembl.org/pub/release-108/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz
```

An augmented CHM13 satellite library and a Kraken2 database are also needed;
see the script for the sources.

### 1.3 Input reads

FASTQ files go in `DATA_DIR` and are named by sample ID:

```
SAMPLE.hifi.fastq.gz                            PacBio HiFi
SAMPLE.ont.fastq.gz                             ONT PromethION
SAMPLE.hic_1.fastq.gz  SAMPLE.hic_2.fastq.gz    Hi-C
```

### 1.4 Sample lists

Plain text files in `LIST_DIR`, one sample ID per line:

```
lists/hifi.list      samples sequenced with PacBio HiFi
lists/ont.list       samples sequenced with ONT PromethION
lists/hic.list       samples with Hi-C data
```

A sample listed in both `hifi.list` and `hic.list` is assembled in Hi-C phasing
mode; if a trimmed ONT FASTQ also exists for it, those reads are passed to
HiFiasm through `--ul`. `lists/qc_pass.list` is written by the `asm_stats` step
and holds the samples whose assemblies meet N50 ≥ 20 Mb and ≤ 1,500 contigs;
the later steps read it.

### 1.5 Environment

```bash
export DATA_DIR=/path/to/fastq        # required
export REF_DIR=/path/to/reference     # required
export OUT_DIR=/path/to/results       # default: ./results
export LIST_DIR=/path/to/lists        # default: ./lists
export THREADS=32
```

---

## 2. Run

Print the available steps:

```bash
bash genome_assembly_and_statistics.sh
```

Run everything in order:

```bash
bash genome_assembly_and_statistics.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash genome_assembly_and_statistics.sh read_qc
bash genome_assembly_and_statistics.sh assembly
bash genome_assembly_and_statistics.sh polishing
bash genome_assembly_and_statistics.sh asm_stats
bash genome_assembly_and_statistics.sh annotation
bash genome_assembly_and_statistics.sh sv_calling
```

Assembly and polishing are long-running; detach them from the terminal:

```bash
nohup bash genome_assembly_and_statistics.sh assembly > assembly.log 2>&1 &
```

Override individual settings on the command line:

```bash
THREADS=64 bash genome_assembly_and_statistics.sh assembly
MIN_N50=15000000 bash genome_assembly_and_statistics.sh asm_stats
```

On a cluster, the per-sample loops inside each step are the natural unit of
parallelism: split the sample list and submit one job per subset, for example

```bash
split -l 1 lists/hifi.list lists/part_
for P in lists/part_*; do
  qsub -N "$(basename "$P")" -pe smp 32 -V -cwd \
    -b y "HIFI_LIST=$P bash genome_assembly_and_statistics.sh assembly"
done
```

Every step writes to `$OUT_DIR/logs/pipeline.log`. Steps are independent, so
one can be rerun after fixing a problem without repeating the others.

---

## 3. Steps

| Step | Description |
|---|---|
| `read_qc` | FastQC, HiFiAdapterFilt, Porechop, Trimmomatic, raw read metrics |
| `assembly` | HiFiasm, HiFiasm Hi-C mode, Flye, HiFiasm ONT mode, contaminant removal |
| `polishing` | Inspector evaluation, correction and re-evaluation |
| `asm_stats` | Inspector summary, QUAST-LG, NGx curves, contiguity filtering, Minigraph misjoins, HMM-Flagger |
| `annotation` | RepeatMasker (two rounds), etrf, sdust, Liftoff |
| `sv_calling` | pbmm2, pbsv, Sniffles, SVIM, SVIM-asm, SURVIVOR merging |
| `all` | All of the above, in order |

Dependencies between steps: `assembly` needs `read_qc`; `polishing` needs
`assembly`; `asm_stats` needs `polishing` and writes `qc_pass.list`;
`annotation` and `sv_calling` need `polishing`, and `sv_calling` also needs the
`qc_pass.list` written by `asm_stats`.

---

## 4. Configuration

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

---

## 5. Output

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

---

## Software

FastQC 0.12.1 · HiFiAdapterFilt 3.0 · Porechop 0.2.4 · Trimmomatic 0.39 ·
HiFiasm 0.24 (0.25 for ONT mode) · Flye 2.9.5 · gfatools · Kraken2 · seqkit ·
Inspector 1.3 · QUAST-LG 5.3 · Minigraph 0.21 · HMM-Flagger 1.0.0 ·
RepeatMasker 4.1.6 · etrf · sdust · Liftoff 1.6.3 · pbmm2 1.17 · pbsv 2.9.0 ·
Sniffles 2.2 · SVIM 2.0.0 · SVIM-asm 1.0.3 · SURVIVOR · samtools · minimap2 ·
jq
