# 1.genome_assembly_and_statistics.sh

Long-read assembly, quality assessment, repeat annotation and structural variant discovery in Korean genomes (Figures 1, 2a, S1–S4; Tables S1–S12, S18).

Part of the [KPR-pipeline](../README.md) repository.

## 1. Set up

### 1.1 Software

Install the tools listed under [Software and versions](#software-and-versions) and make sure they are on
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
HiFiasm through `--ul`. `lists/qc_pass.list is written by the `asm_stats` step.
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


### Acrocentric short-arm intervals

Supply a zero-based, half-open CHM13 BED containing the short arms of chr13, chr14, chr15, chr21 and chr22. Set `ACROCENTRIC_BED` to this file for the study's misjoin filtering. `PAFTOOLS` specifies the Minimap2 v2.30 `paftools.js` executable.

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


## 3. Steps

| Step | Description |
|---|---|
| `read_qc` | FastQC, HiFiAdapterFilt, Porechop, Trimmomatic, raw read metrics |
| `assembly` | HiFiasm, HiFiasm Hi-C mode, Flye, HiFiasm ONT mode, contaminant removal |
| `polishing` | Inspector evaluation, correction and re-evaluation |
| `asm_stats` | Inspector summary, QUAST-LG, NGx curves, contiguity filtering, Minigraph/paftools.js misjoins, filtered assemblies, HMM-Flagger |
| `annotation` | RepeatMasker (two rounds), etrf, sdust, Liftoff |
| `sv_calling` | pbmm2, pbsv, Sniffles, SVIM, SVIM-asm; SVIM-asm plus at least two of the three read-based callers, followed by SURVIVOR cohort merging |
| `all` | All of the above, in order |

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
| `SV_MIN_SIZE` | `51` | minimum SV length (bp) |
| `SV_MIN_DP` / `SV_MIN_ALT_SUPPORT` | `7` / `3` | read-caller depth and ALT support |
| `SNIFFLES_MIN_GQ` / `SVIM_MIN_QUAL` | `20` / `10` | caller-specific quality thresholds |
| `TRUVARI_REFDIST` | `500` | Truvari reference distance (bp) |
| `TRUVARI_PCTSEQ` / `TRUVARI_PCTSIZE` | `0.70` / `0.70` | sequence and size similarity |
| `EXPECTED_SV_SAMPLES` | `14` | manuscript cohort size; `0` disables the count check |
| `ACROCENTRIC_BED` | unset | CHM13 acrocentric short-arm BED; required for classification and filtering |
| `PAFTOOLS` | `paftools.js` | Minimap2 v2.30 misjoin utility |
| `QC_ASM_OUT` | `$OUT_DIR/polishing/qc_assemblies` | QC-passed assemblies for downstream analysis |

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
- `sv/filtered/SAMPLE.CALLER.vcf.gz` — per-caller PASS, autosomal,
  ≥51-bp SVs with read-based quality thresholds as applicable
- `sv/filtered/SAMPLE.CALLER.filter_audit.tsv` — filtering reasons/counts,
  including missing depth or ALT-read support
- `sv/truvari/SAMPLE/{pbsv,sniffles,svim}/tp-base.vcf.gz` — the
  SVIM-asm variants matched by each independent read-based caller
- `sv/merged/SAMPLE.merged.vcf` — per-sample SVIM-asm SVs matched by
  at least two of pbsv, Sniffles and SVIM
- `sv/merged/SAMPLE.caller_support.tsv` — 0/1/2/3 matched read-caller counts
- `sv/merged/shared_SV.vcf` — ≥51-bp SVs observed in all 14 QC-passed
  individuals after SURVIVOR merging at 100 bp

The polished assemblies in `polishing/corrected/` are the input for
`3.pangenome_construction_and_statistics.sh`.

The downstream graph must use **the actual two retained haplotype
assemblies per selected individual**, not duplicate FASTA records
or both raw and polished versions. Check sequence identifiers,
sample provenance and haplotype FASTA indexing before constructing
seqfiles. Record versions of the retained FASTA files, the
`qc_pass.list`, contamination-screen reports and key
quality summary tables together with graph input manifests.

`polishing/qc_assemblies/` contains the retained haplotype FASTAs. `statistics/misjoin/interchromosomal_misjoins.tsv` records the candidate pairs, their CHM13 coordinates and classification. `lists/qc_pass.list` contains individuals with two passing haplotypes.

## 6. Selection of samples for long-read sequencing

We aligned paired-end SRS from 308 Korean individuals to GRCh38.p14 with bwa-mem2 v2.2.1, followed by GATK v4.3.0.0 duplicate removal, base quality recalibration, HaplotypeCaller and joint genotyping. Variants with MQ <30, QUAL <30, mean DP <7 or GQ <10 were filtered; VQSR used a truth sensitivity of 99.7%. Biallelic SNPs with <5% missingness and MAF ≥0.05 were phased with SHAPEIT5 v5.1.1 and imputed with IMPUTE5. The pre-MAF-filter QC SNP dataset was retained for AMY population genetics.

The merged KOR/1000 Genomes panel contained 2,856 individuals (308 + 2,548). We ran fastSTRUCTURE v1.0 for K=2–5 and visualized K=4 in Figure 1a. PLINK v2.0 calculated PC1–PC30; the explained variance of PC i was its eigenvalue divided by the trace of the genomic relationship matrix (sum of all 2,856 diagonal elements), multiplied by 100. The study selected 16 core individuals using GenoCore v1.0 (-cv 99 -d 0.001). SRS structural variation was characterized using Manta INS/DUP (QUAL ≥30) and Delly PASS DEL (Tables S1–S4).

## 7. Genome assembly and quality assessment

PacBio Sequel II/Revio reads were adapter-trimmed with HiFiAdapterFilt v3.0 and assembled with HiFiasm v0.24; Hi-C data were used for Hi-C phasing. ONT reads were trimmed with Porechop v0.2.4 and Trimmomatic v0.39 and assembled using Flye v2.9.5; HiFiasm v0.25 ONT-mode assemblies provided an assembler comparison. Non-human contigs, including EBV sequences, were removed. Inspector v1.3 was used for polishing and self-alignment assessment. We retained both haplotypes of individuals passing N50 ≥20 Mb and ≤1,500 contigs. QUAST v5.3 assessed the assemblies against GRCh38 and CHM13 (Tables S5–S9).

## 8. Interchromosomal misjoins

Minigraph v0.21 (-cxasm) aligned the QC-passed haplotypes to CHM13. Minimap2 v2.30 paftools.js misjoin used a minimum alignment block length of 1 Mb. Both alignment intervals were intersected with the CHM13 acrocentric short-arm BED. Pairs for which both intervals overlapped acrocentric short arms were retained; all other interchromosomal candidate pairs were classified as assembly errors and the corresponding contigs excluded. The -e option outputs the alignment pairs for classification.

~~~bash
export ACROCENTRIC_BED=/path/to/CHM13.acrocentric_short_arms.bed
bash 1.genome_assembly_and_statistics.sh asm_stats
~~~

Use the resulting QC_ASM_OUT FASTAs for annotation, SVIM-asm and the stage 3 assembly manifests. HMM-Flagger v1.0.0 assessed read coverage of the filtered assemblies. Erroneous, duplicated, collapsed and unknown sequence lengths were summed for the unreliable fraction; haploid regions provided the reliable annotation intervals (Figure 1f–g; Tables S10–S11).

## 9. Repeat annotation

For the assembly repeat analysis we used RepeatMasker v4.1.6 with the human repeat library, followed by an augmented CHM13 satellite library. Previously assigned repeats were hard-masked before the next annotation round. etrf identified VNTRs with repeat units ≥7 bp; sdust detected low-complexity regions after hard-masking the earlier repeat classes. The combined repeat intervals were applied to the original assembly using soft masking to prepare the stage 2 input. The annotation step supplies individual repeat annotations and a Liftoff comparison; combine the repeat intervals separately for the final soft-masked FASTA. Count repeat sites and annotated bases by Alu, L1, L2, SVA, ERV, other TE, VNTR, STR, LCR and satellite classes (Figure 2a; Table S12).

## 10. Structural variant detection

PacBio LRS was aligned to GRCh38 using pbmm2 v1.17. We called SVs with pbsv v2.9.0, Sniffles v2.2, SVIM v2.0.0 and SVIM-asm v1.0.3. All calls were autosomal PASS SVs ≥51 bp. pbsv required DP ≥7 and ALT reads ≥3; Sniffles required GQ ≥20, DR+DV ≥7 and DV ≥3; SVIM required QUAL ≥10, DP ≥7 and ALT reads ≥3. SVIM-asm calls matching at least two read-based callers were selected by Truvari bench (breakpoint distance ≤500 bp; size and sequence similarity ≥70%). SURVIVOR merged selected calls with a 100-bp breakpoint distance. Variants carried by all 14 individuals were counted in Table S18. Stage 5 uses a separate, more stringent truth-set procedure for seven benchmark samples.

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
| Supporting tools | Kraken2; seqkit; samtools; minimap2; gfatools; SURVIVOR; Truvari bench; bcftools; Python 3; jq |


Minimap2 v2.30 (`paftools.js`), fastSTRUCTURE v1.0, PLINK v2.0, GenoCore v1.0, SHAPEIT5 v5.1.1, IMPUTE5.
