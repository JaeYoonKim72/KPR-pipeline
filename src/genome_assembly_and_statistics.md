# 1.genome_assembly_and_statistics.sh

## Scope and rationale

A high-quality Korean pangenome depends on accurate, haplotype-resolved
de novo assemblies. Short reads provide reliable genotypes across much of
the genome, but cannot consistently traverse long repeats, highly
homologous segmental duplications or complex structural variants.
PacBio high-fidelity (HiFi) reads, complemented by long Oxford Nanopore
(ONT) reads and chromatin-conformation data where available, provide
sequence continuity and haplotype information needed to resolve
structures that are absent or collapsed in GRCh38. Reference alignment
is used for evaluation, rather than as the basis of the assembly
itself [1–3].

This workflow processes PacBio Sequel II/Revio, ONT PromethION and
Hi-C reads through platform-specific quality control, assembly,
polishing, reference-based and reference-free validation, repeat/gene
annotation and structural variant discovery. The distinction between
`read_qc`, `assembly`, `polishing`, and `asm_stats` is deliberate:
high raw-read accuracy or high contiguity alone does not establish that
a genome assembly is complete, correctly phased or free of collapsed
repeats, contaminant sequence and interchromosomal misjoins.

In the KPR manuscript, long-read datasets were collected for
22 Korean individuals: 19 newly sequenced individuals (12 Revio,
seven PromethION) and three previously generated PacBio Sequel II
datasets. Assemblies were retained for pangenome construction
only when contig N50 was at least 20 Mb and the number of contigs
was at most 1,500. The retained panel consisted of **28 haploid
assemblies from 14 individuals**. All seven PromethION-derived
assemblies in this particular dataset failed the contig-count
criterion; this outcome should not be interpreted as a general
technology comparison. The workflow reports the underlying QC
metrics so that exclusions and retained haplotype assemblies can
be traced [1].

This stage produces the polished FASTA sequences and validation
evidence for [pangenome construction](pangenome_construction_and_statistics.md).
Population-level comparisons are not calculated by this stage.

Part of the [KPR-pipeline](../README.md) repository.

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
bash 1.genome_assembly_and_statistics.sh
```

Run everything in order:

```bash
bash 1.genome_assembly_and_statistics.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash 1.genome_assembly_and_statistics.sh read_qc
bash 1.genome_assembly_and_statistics.sh assembly
bash 1.genome_assembly_and_statistics.sh polishing
bash 1.genome_assembly_and_statistics.sh asm_stats
bash 1.genome_assembly_and_statistics.sh annotation
bash 1.genome_assembly_and_statistics.sh sv_calling
```

Assembly and polishing are long-running; detach them from the terminal:

```bash
nohup bash 1.genome_assembly_and_statistics.sh assembly > assembly.log 2>&1 &
```

Override individual settings on the command line:

```bash
THREADS=64 bash 1.genome_assembly_and_statistics.sh assembly
MIN_N50=15000000 bash 1.genome_assembly_and_statistics.sh asm_stats
```

On a cluster, the per-sample loops inside each step are the natural unit of
parallelism: split the sample list and submit one job per subset, for example

```bash
split -l 1 lists/hifi.list lists/part_
for P in lists/part_*; do
  qsub -N "$(basename "$P")" -pe smp 32 -V -cwd \
    -b y "HIFI_LIST=$P bash 1.genome_assembly_and_statistics.sh assembly"
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

### 3.1 Read quality control and platform-specific assembly

The PacBio HiFi branch removes residual adapter sequence with
HiFiAdapterFilt and assembles long accurate reads using HiFiasm.
For samples with Hi-C links, HiFiasm's Hi-C mode distinguishes
maternal and paternal haplotype paths; the underlying per-sample
read identifiers must remain consistent across technologies.
The ONT branch removes adapters with Porechop, trims low-quality
bases and assembles with Flye; an alternative HiFiasm ONT
mode is available to assess assembler-dependent contiguity.
Where a suitable trimmed ONT readset belongs to a Hi-C-enabled
HiFi sample, the wrapper can pass those reads to HiFiasm
with `--ul`. This does not imply that every ONT sample has
Hi-C data or belongs in a haplotype-resolved final panel.

### 3.2 Contaminant screening and polishing

Assemblies of lymphoblastoid cell-line-derived genomes can
contain viral or other non-human sequence. Kraken2-based
classification and contaminant removal precede downstream
annotation. Inspector assesses local assembly errors, applies
polishing and re-evaluates the polished sequences; retain both
pre- and post-polishing summaries for audit. Polishing is not
a substitute for checking long-range structural accuracy.

### 3.3 Contiguity, completeness and structural integrity

Each haplotype assembly is evaluated across complementary
dimensions:

| Measure | Interpretation | Manuscript relationship |
|---|---|---|
| Contig count and N50 | Sequence fragmentation and long-range continuity | PASS requires N50 >=20 Mb and <=1,500 contigs |
| Total assembled length | Missing or excess sequence | Compare to ungapped GRCh38/CHM13; neither is a complete truth set |
| Inspector self-alignment and QV | Local consensus quality and coverage | Retain self-depth and small-error summaries |
| QUAST-LG aligned/unaligned sequence | Reference representation and possible non-reference sequence | Evaluate independently against GRCh38 and CHM13 |
| Duplication ratio | Potential duplicated or collapsed content | Investigate discordance with read depth |
| Minigraph misjoins | Candidate interchromosomal connections | Review acrocentric chromosome alignments separately |
| HMM-Flagger states | Reliable versus collapsed, duplicated or erroneous sequence | Preserve flagged bases as a QC limitation |

The manuscript retained assemblies with a mean contig N50 of
50.34 Mb and a mean assembled length of 3.02 Gb.
Mean coverage of GRCh38 and CHM13 was 96.06% and 93.20%,
respectively, and HMM-Flagger identified 43.89 Mb (1.46%)
of unreliable assembly sequence on average. These are
**observed manuscript results**, not universal software
thresholds. Acrocentric alignments, including apparent
interchromosomal joins, require particular caution because
homologous satellite and segmental-duplication sequence
can produce complex mappings.

### 3.4 Repeats, gene models and copy-number gains

RepeatMasker annotates interspersed elements, satellites and
other repeat families. Separate etrf and sdust passes add VNTR
and low-complexity sequence after accounting for previously
masked intervals. Liftoff projects GENCODE gene/transcript
models to the assembly. Report the percentage of models
represented and the additional copies of each protein-coding
gene relative to the selected reference. An inferred extra
gene model is a copy-number-gain candidate, not an
independently validated expression or functional result.
The manuscript compared these quantities among the
KPR, HPRC-EAS and HPRC-nEAS assemblies.

### 3.5 Structural-variant discovery and integration

The SV workflow combines read- and assembly-derived
evidence. PacBio reads are aligned with pbmm2 and
interrogated by pbsv, Sniffles and SVIM; assembly-based
candidates are identified by SVIM-asm. The wrapper
merges near-equivalent calls with a 100-bp breakpoint
tolerance and reports shared >=51-bp events across
retained samples. Such shared events should be
interpreted relative to the selected reference,
calling algorithms and available samples, not
as variants fixed in all Koreans. Keep insertion,
deletion, duplication and inversion types separate
where the underlying callsets support them.

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
`3.pangenome_construction_and_statistics.sh`.

The downstream graph must use **the actual two retained haplotype
assemblies per selected individual**, not duplicate FASTA records
or both raw and polished versions. Check sequence identifiers,
sample provenance and haplotype FASTA indexing before constructing
seqfiles. Record versions of the retained FASTA files, the
`qc_pass.list`, contamination-screen reports and key
quality summary tables together with graph input manifests.

---

## References

1. Kim JY, et al. *Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus.* KPR manuscript, revision 2026; sections "Genome assembly and quality assessment" and "Repeat and gene annotation".
2. Cheng H, et al. (2024). Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph. *Nature Methods*. https://doi.org/10.1038/s41592-024-02269-8
3. Kolmogorov M, et al. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*. https://doi.org/10.1038/s41587-019-0072-8
4. Chen Y, et al. (2021). Accurate long-read de novo assembly evaluation with Inspector. *Genome Biology*. https://doi.org/10.1186/s13059-021-02527-4
5. Gurevich A, et al. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btt086

## Software and versions

The versions below follow the manuscript Methods. Tools that implement
auxiliary command-line operations without a version in the manuscript
are listed separately; their precise versions should be captured in
the execution log rather than assumed.

| Analysis stage | Software (manuscript version) |
|---|---|
| Raw-read QC and trimming | FastQC 0.12.1; HiFiAdapterFilt 3.0; Porechop 0.2.4; Trimmomatic 0.39; Guppy 6.5.7 (ONT basecalling) |
| Genome assembly | HiFiasm 0.24 (HiFi/Hi-C); HiFiasm 0.25 (ONT evaluation); Flye 2.9.5 |
| Polishing and assembly assessment | Inspector 1.3; QUAST 5.3; HMM-Flagger 1.0.0; Minigraph 0.21 |
| Repeats and annotation | RepeatMasker 4.1.6; Liftoff 1.6.3; etrf; sdust |
| Read-based structural variants | pbmm2 1.17; pbsv 2.9.0; Sniffles 2.2; SVIM 2.0.0 |
| Assembly-based structural variants | SVIM-asm 1.0.3 |
| Supporting tools | Kraken2; seqkit; samtools; minimap2; gfatools; SURVIVOR; jq |

The manuscript's GATK-based short-read sample-selection workflow is
separate from this long-read assembly stage; it used bwa-mem2 2.2.1,
GATK 4.3.0.0, SHAPEIT5 5.1.1, fastSTRUCTURE 1.0,
PLINK 2.0 and GenoCore 1.0 [1].
