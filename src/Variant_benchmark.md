# 5.Variant_benchmark.sh

Variant-calling accuracy was assessed by comparing SNPs, small insertions and
deletions, and structural insertions and deletions from graph-based and
GRCh38-based call sets against reference truth sets. The workflow includes
the GIAB v4.2.1 small-variant truth set, T2T-HG002-Q100 v1.1 structural-variant
truth set, and long-read-derived truth sets for seven KOR individuals (Figures 4a–b and 4d–e; Tables a1, a3 and S29–S30)
(Wagner et al., 2022; English et al., 2022).

The pipeline uses hap.py with RTG vcfeval for small variants and Truvari
bench for structural variants. Matching alleles and genotypes are evaluated
separately. Mendelian genotype inheritance inconsistency is assessed in
the Ashkenazi Jewish (HG002–HG004) and Han Chinese (HG005–HG007) GIAB trios.
All comparisons use the same reference assembly and are limited to autosomes.
The workflow is part of the [KPR-pipeline](../README.md).

---

## 1. Set up

### 1.1 Software

Install bcftools, samtools, bgzip/tabix, RTG Tools, hap.py, Truvari,
Python and pysam. **Truvari v5.5.0 requires Python ≥3.12.** The
version of each executable should be recorded with the benchmark output.
The required benchmark settings are listed under
[Software and versions](#11-software-and-versions).

### 1.2 Cohorts and reference data

The seven KOR samples with short-read, PacBio long-read and assembly
data are KOR-001, KOR-002, KOR-006, KOR-007, KOR-008, KOR-009 and
KOR-010. A plain-text sample list contains one ID per row:

~~~text
KOR-001
KOR-002
KOR-006
KOR-007
KOR-008
KOR-009
KOR-010
~~~

GIAB analyses use HG002 (Ashkenazi Jewish child) with the
small-variant truth set **v4.2.1** and structural-variant truth set
**T2T-HG002-Q100 v1.1**, each mapped to the GRCh38 coordinate system
used for comparison. HG005 is the child of the Han Chinese trio.
There is no corresponding HG005 SV truth set for the comparison
described in this manuscript.

The inputs require an indexed GRCh38 FASTA and an RTG sequence
data file built from the same reference:

~~~bash
samtools faidx /path/to/GRCh38.fa
rtg format -o /path/to/GRCh38.sdf /path/to/GRCh38.fa
~~~

### 1.3 Truth sets and confidence BED files

The KOR truth-set inputs follow the existing main-branch pipeline layout:

~~~text
truth_input/
├── hifi_deepvariant/
│   ├── SAMPLE.hifi_dv.biSNP.addfilter.VAF0208.vcf.gz
│   └── SAMPLE.hifi_dv.Small.vcf.gz
├── long_read_sv/
│   ├── svim_asm_SAMPLE.PASS.SVLEN51.chr.vcf.gz
│   ├── pbsv_SAMPLE.PASS.SVLEN51.chr.vcf.gz
│   ├── sniffles_SAMPLE.PASS.SVLEN51.chr.vcf.gz
│   └── svim_SAMPLE.PASS.SVLEN51.chr.vcf.gz
└── confidence_beds/
    ├── SAMPLE.small.bed
    └── SAMPLE.sv.bed
~~~

**The confidence BEDs must come from the actual truth-set
construction or evaluation design.** The script does not infer
high-confidence regions from the presence of variants. It requires
a separate BED for each sample and truth class. For GIAB, use the
matching official benchmark confidence BED and place it under
`CONF_BED_DIR` as `HG002.small.bed` or `HG002.sv.bed`.
The SV BED should cover only the eligible GRCh38 regions from the
SV truth set. The `--includebed` option requires SV intervals to
be contained within an eligible region, rather than merely touching it.

Using the entire autosomal reference as the confidence interval is
**not** equivalent to a GIAB or KOR high-confidence BED. An explicit
`ALLOW_AUTOSOME_AS_CONFIDENCE=1` option exists for exploratory runs;
it is off by default and must not be used to label exploratory values
as confidence-restricted manuscript results.

### 1.4 Query input

~~~text
query_input/
├── kpr/         kpr.SAMPLE.snp.vcf.gz       kpr.SAMPLE.indel.vcf.gz
│               kpr.SAMPLE.sv.DEL.vcf.gz    kpr.SAMPLE.sv.INS.vcf.gz
├── hprc_kpr/    analogous files
├── hprc/        analogous files
├── deepvariant/ deepvariant.SAMPLE.snp.vcf.gz
│               deepvariant.SAMPLE.indel.vcf.gz
├── manta/       manta.SAMPLE.sv.DEL.vcf.gz manta.SAMPLE.sv.INS.vcf.gz
└── delly/       delly.SAMPLE.sv.DEL.vcf.gz
~~~

Graph-based short-read calls are produced in GRCh38 coordinates by
[stage 6](pangenome_short_read_mapping.md); the baseline uses
DeepVariant for small variants, Manta for INS and DEL, and Delly
for DEL where available. For KOR comparisons, the construction of
KPR-containing graphs excludes the focal individual.

### 1.5 Environment

~~~bash
export REF_DIR=/path/to/reference
export TRUTH_DIR=/path/to/truth_input
export QUERY_DIR=/path/to/query_input
export LIST_DIR=/path/to/lists
export SAMPLE_LIST="$LIST_DIR/benchmark_samples.list"
export CONF_BED_DIR="$TRUTH_DIR/confidence_beds"
export OUT_DIR=/path/to/results
export JOBS=8
export THREADS=4
export TRUTH_MODE=kor
~~~

---

## 2. Run

### 2.1 Check and run

~~~bash
bash 5.Variant_benchmark.sh --help
bash 5.Variant_benchmark.sh check
bash 5.Variant_benchmark.sh all
~~~

The KOR mode builds the three small-variant truth classes and three
SV support levels before benchmarking. A single analysis stage may
be called when its prerequisites are already complete:

~~~bash
bash 5.Variant_benchmark.sh small_truth
bash 5.Variant_benchmark.sh sv_truth_filter
bash 5.Variant_benchmark.sh sv_truth_support
bash 5.Variant_benchmark.sh sv_truth_readsupport
bash 5.Variant_benchmark.sh small_query
bash 5.Variant_benchmark.sh small_bench
bash 5.Variant_benchmark.sh sv_query
bash 5.Variant_benchmark.sh sv_bench
bash 5.Variant_benchmark.sh fn_classify
bash 5.Variant_benchmark.sh summary
~~~

### 2.2 GIAB HG002 benchmark

GIAB mode reads supplied, GRCh38-compatible GIAB truth VCFs
rather than deriving a long-read truth set from the seven KOR
samples. Set `SAMPLE_LIST` to a file containing only `HG002`.
Specify the v4.2.1 small-variant and v1.1 SV truth VCFs, and
their matching confidence BEDs.

~~~bash
export TRUTH_MODE=giab
export SAMPLE_LIST=/path/to/lists/HG002.list
export GIAB_SMALL_VCF=/path/to/GIAB_HG002_v4.2.1.GRCh38.vcf.gz
export GIAB_SV_VCF=/path/to/T2T-HG002-Q100_v1.1.GRCh38.vcf.gz
export SV_SUPPORT_LEVELS=ge2
export CONF_BED_DIR=/path/to/giab_confidence_beds
bash 5.Variant_benchmark.sh all
~~~

The same query names are expected for HG002 as for the KOR samples.
Keep the KOR and GIAB runs in **different `OUT_DIR` locations**,
so that their truth sets and summaries cannot overwrite each other.

### 2.3 GIAB trio Mendelian analysis

The Mendelian stage uses separate genotype VCFs from the child,
father and mother, including *genuinely called* homozygous-reference
genotypes. The three DeepVariant gVCFs should first be represented
as sample-specific genotype-called VCFs with the required genotype
information; missing genotypes must not be imputed as reference.
For the graph SV comparison use the `vg call` genotype VCFs;
for GRCh38 SVs, set the INS input directory to Manta and
the DEL input directory to Delly.

~~~bash
export MEND_TRIO=AJ
export MEND_ARM=graph
export MEND_REFERENCE=kpr
export MEND_VCF_DIR=/path/to/graph_trio_genotype_vcfs
export MEND_TRUTH_SMALL=/path/to/HG002_small_truth.vcf.gz
export MEND_CONF_SMALL=/path/to/HG002_small_truth.bed
export MEND_TRUTH_SV=/path/to/HG002_sv_truth.vcf.gz
export MEND_CONF_SV=/path/to/HG002_sv_truth.bed
bash 5.Variant_benchmark.sh mendelian
~~~

For the Han Chinese trio, set `MEND_TRIO=HC` and supply the
HG005 small-variant truth. Do **not** set `MEND_TRUTH_SV` or
`MEND_CONF_SV`; the observed SV inconsistency counts can be
reported, but SV truth-based adjudication is unavailable.

For GRCh38-based calls, specify `MEND_ARM=linear`,
`MEND_SMALL_DIR`, `MEND_INS_DIR` and `MEND_DEL_DIR`
to point to the corresponding trio VCF folders.
The referenced helper is
[`src/trio_mendelian_inconsistency.sh`](trio_mendelian_inconsistency.sh),
which calls RTG `mendelian` and the child-truth adjudication helper.

---

## 3. Steps

| Step | Command | Description |
|---|---|---|
| 1 | `small_truth` | KOR HiFi DeepVariant truth SNP, INS and DEL sets |
| 2 | `sv_truth_filter` | PASS/autosome/51-bp and caller-specific quality filtering |
| 3 | `sv_truth_support` | SVIM-asm variants matched by ≥1, ≥2 or 3 read-based callers |
| 4 | `sv_truth_readsupport` | Require ALT supporting reads ≥5 and DP ≥10 in one matching caller |
| 5 | `small_query` | Normalize, split and classify small query variants |
| 6 | `small_bench` | hap.py with RTG vcfeval, annotated VCF and summary.csv |
| 7 | `sv_query` | Normalize and classify query INS/DEL of 51–50,000 bp |
| 8 | `sv_bench` | Truvari v5.5.0 genotype-aware matching within confidence BED |
| 9 | `fn_classify` | Detect FN causes and audit unresolved GBZ evidence |
| 10 | `summary` | Per-sample metrics and sample-level means |
| Alternative | `giab_truth` | Prepare supplied HG002 GIAB truth VCFs |
| Separate | `mendelian` | GIAB trio genotype inheritance inconsistency |

The SV support filter and benchmarking steps use different
Truvari thresholds, as described below. The Mendelian stage
is separate from the seven-sample benchmarking average.

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `TRUTH_MODE` | `kor` | `kor` or `giab` |
| `SAMPLE_LIST` | `lists/benchmark_samples.list` | Cohort sample IDs |
| `GRAPHS` | `kpr hprc_kpr hprc` | Graph call sets |
| `LINEAR_SMALL` | `deepvariant` | GRCh38 small-variant query |
| `SV_SUPPORT_LEVELS` | `ge1 ge2 ge3` | KOR SV truth support levels; primary `ge2` |
| `TRUTH_MIN_GQ` / `TRUTH_MIN_DP` | 20 / 7 | HiFi DeepVariant truth genotype |
| `HET_VAF_MIN` / `HET_VAF_MAX` | 0.30 / 0.70 | Heterozygous truth VAF |
| `HOM_VAF_MIN` | 0.90 | Homozygous ALT truth VAF |
| `SUPPORT_REFDIST` | 500 | SVIM-asm/read-based caller Truvari distance (bp) |
| `SUPPORT_PCTSIZE` / `SUPPORT_PCTSEQ` | 0.70 / 0.70 | SV truth construction matching |
| `READSUP_MIN_ALT` / `READSUP_MIN_DP` | 5 / 10 | Strong support in at least one matching read caller |
| `TV_REFDIST` | 500 | Benchmark Truvari distance (bp) |
| `TV_PCTSIZE` / `TV_PCTSEQ` / `TV_PCTOVL` | 0.80 / 0.80 / 0.80 | Benchmark allele similarities |
| `TV_SIZEMIN` / `TV_SIZEFILT` / `TV_SIZEMAX` | 51 / 30 / 50000 | Benchmark SV size handling |
| `CONF_BED_DIR` | `TRUTH_DIR/confidence_beds` | Per-sample `SAMPLE.small.bed`, `SAMPLE.sv.bed` |
| `ALLOW_AUTOSOME_AS_CONFIDENCE` | 0 | Exploratory whole-autosome comparison, not manuscript confidence |
| `GBZ_EVIDENCE_DIR` | Unset | Independently prepared Truvari results for personalized-GBZ allele representation |
| `JOBS` / `THREADS` | 8 / 4 | Parallel jobs and hap.py threads |
| `MEND_TRIO` / `MEND_ARM` | AJ / graph | Trio and reference arm |
| `MEND_MIN_GQ` / `MEND_MIN_DP` | 20 / 10 | Mendelian genotype filter; linear SV DP unavailable |

---

## 5. Output

~~~text
results/benchmarking/
├── autosomes.bed
├── truth/
│   ├── small/SAMPLE/                   SNP/DEL/INS truth VCFs
│   └── sv/SAMPLE/                      caller-filtered and supported SV sets
│       └── ge1/, ge2/, ge3/           read-support evidence and INS/DEL truth
├── query/
│   ├── small/SAMPLE/                   preprocessed query VCFs
│   └── sv/SAMPLE/                      preprocessed SV query VCFs
├── bench/
│   ├── small/SAMPLE/                   hap.py *.summary.csv and annotated *.vcf.gz
│   └── sv/SAMPLE/DATASET/LEVEL/TYPE/   Truvari summary.json, tp-base, tp-comp, fp, fn
├── summary/
│   ├── kor_small_detail.tsv
│   ├── kor_small_mean.tsv
│   ├── kor_sv_detail.tsv
│   ├── kor_sv_mean.tsv
│   └── fn_failure_classes.tsv
└── mendelian/REFERENCE.TRIO/
    ├── SNP.annot.vcf.gz, Indel.annot.vcf.gz, INS.annot.vcf.gz, DEL.annot.vcf.gz
    ├── CLASS.mcv.vcf.gz
    ├── CLASS.child_mcv.vcf.gz
    ├── CLASS.truth.vcf.gz or CLASS.truvari/
    ├── mendelian_summary.tsv
    └── mendelian_variant_adjudication.tsv
~~~

The `mendelian_summary.tsv` retains total, child non-reference and
testable variants, observed Mendelian-consistent and -inconsistent
counts and rates, adjudicated true/false inconsistent variants,
and the false-inconsistent rate. An unresolved or unsupported
truth comparison is recorded separately, not silently counted false.

---

## 6. Construction of KOR truth sets

### 6.1 Small variants

PacBio LRS is aligned to GRCh38 and variants are called with
DeepVariant v1.6. SNPs and 1–50 bp simple insertions/deletions
must meet `FILTER=PASS`, `GQ ≥20`, `DP ≥7` and the VAF
requirements: **0.30–0.70** for heterozygotes and **≥0.90**
for homozygous ALT genotypes. Multiallelic variants are split
and normalized with `bcftools norm -m -any`
(Danecek et al., 2021).

### 6.2 SV callers and support

SVIM-asm is the assembly-based anchor. Each of pbsv, Sniffles
and SVIM must pass the same initial filters described in
[stage 1](genome_assembly_and_statistics.md): autosomes,
`FILTER=PASS` and **51–50,000 bp**, with the read-based
quality filters (pbsv: DP ≥7, ALT reads ≥3; Sniffles: GQ ≥20,
DR+DV ≥7, DV ≥3; SVIM: QUAL ≥10, DP ≥7, ALT reads ≥3).

Each read-based caller is compared separately with SVIM-asm
using Truvari bench, with a 500-bp maximum breakpoint
difference, **≥70%** size similarity and **≥70%** sequence
similarity. The `ge1`, `ge2`, and `ge3` sets contain
SVIM-asm variants supported by at least one, at least two,
or all three read-based callers. The **ge2** set is used
for the main SV truth benchmark (Table S30).

Only INS and DEL are included. At least one read caller that
supports the SVIM-asm allele must have both **ALT-supporting
reads ≥5** and **DP ≥10**. This support condition is
additional to the caller-level screening. The support table
and per-variant evidence are preserved for auditing.

---

## 7. Variant benchmarking and metrics

### 7.1 Variant preprocessing and matching

Only autosomal `FILTER=PASS` query records are retained.
Multiallelic records are decomposed using `bcftools norm
-m -any`. Variants are classified as SNPs, small INS/DEL
(1–50 bp), or SV INS/DEL (51–50,000 bp).

Small-variant comparisons use hap.py with the RTG vcfeval
engine, which allows haplotype-equivalent representations
(Cleary et al., 2015). Truvari bench v5.5.0
(English et al., 2022) compares SVs within the relevant
truth-set high-confidence BED, using:

~~~bash
--includebed CONFIDENCE.bed --refdist 500 --pctsize 0.8 \
--pctseq 0.8 --pctovl 0.8 --pick single \
--sizemin 51 --sizefilt 30 --sizemax 50000 --passonly
~~~

Truvari produces `summary.json`. hap.py produces
`summary.csv` and, with `-V`, an annotated VCF.

### 7.2 Genotype-aware precision, recall and F1

A true positive requires both allele and genotype concordance.

For small variants, use the `METRIC.Precision`,
`METRIC.Recall`, and `METRIC.F1_Score` values from
hap.py (`TRUTH.TP` for genotype-concordant truth calls).

For SVs, Truvari distinguishes matching alleles from
matching genotypes. Let `B` denote `base cnt` (truth
variants) and `C` denote `comp cnt` (query variants):

\[
\mathrm{Recall}_{GT}=\frac{\mathrm{TP\text{-}base\_TP\text{-}gt}}{B}
\]

\[
\mathrm{Precision}_{GT}=\frac{\mathrm{TP\text{-}comp\_TP\text{-}gt}}{C}
\]

\[
F_1=\frac{2\,\mathrm{Precision}_{GT}\,\mathrm{Recall}_{GT}}
{\mathrm{Precision}_{GT}+\mathrm{Recall}_{GT}}
\]

A genotype-discordant SV match contributes one truth-side
FN (`TP-base_FP-gt`) and one query-side FP
(`TP-comp_FP-gt`). These are added to the unmatched
truth FN and query FP when reporting GT-aware errors.
The allele-level `TP-base` and `TP-comp` counts must
not be substituted for GT-aware TP.

### 7.3 Presence and genotype concordance

**Presence** represents the fraction of truth variants
with a matching query allele, independent of genotype:

\[
\mathrm{Presence}_{small}
=\frac{\mathrm{TRUTH.TP}+\mathrm{BK{=}am}}{\mathrm{TRUTH.TOTAL}}
\]

\[
\mathrm{Presence}_{SV}
=\frac{\mathrm{TP\text{-}base}}{\mathrm{base\ cnt}}
\]

Here, `BK=am` identifies genotype-only mismatches in
the hap.py annotated truth records. **GT concordance**
is the fraction of present alleles with concordant
genotypes: `TRUTH.TP/(TRUTH.TP + BK=am)` for small
variants, and `gt_concordance` for SVs.

Values are calculated per sample and then averaged
across the seven KOR individuals. The pipeline does not
perform a Wilcoxon test.

---

## 8. False-negative classification

### 8.1 Small variants

Hap.py's annotated VCF is inspected among truth-side
false negatives. A missing corresponding query call
(`BK=.`) is a **detection failure**. Allele (`BK=lm`)
or genotype (`BK=am`) discordance is a **matching
failure**. Haplotype-aware vcfeval matching is retained;
equivalent allelic representations are not counted
as mismatches solely due to VCF record layout.

### 8.2 Structural variants

Truvari's `fn.vcf.gz` contains unmatched truth SVs.
For each, a query call within ±500 bp is sought.
An unmatched SV without a query in this window is
classified as a detection failure. Nearby calls that
disagree in allele or genotype can indicate a matching
failure. Genotype-discordant `tp-base` records are
also matching failures in GT-aware analysis.

For graph-based analyses, the sample-specific personalized
GBZ used in k-mer-based haplotype sampling and mapping must
additionally be examined for representation of the truth
allele, using the same Truvari matching criteria.
**If the personalized-GBZ allele comparison is not
available, a nearby query call alone does not establish
that the allele was present in the graph.** The helper
reports these cases as `unresolved` rather than
assigning a definitive detection/matching cause.

The optional `GBZ_EVIDENCE_DIR` accepts previously
validated Truvari `tp-base.vcf.gz` evidence under
`SAMPLE/DATASET.TYPE/`, derived by comparing the truth
FN allele with the personalized-GBZ allele inventory.
The script does not construct a personalized GBZ or
assert that graph allele presence was checked when
such evidence is absent.

---

## 9. Mendelian genotype inheritance inconsistency

### 9.1 GIAB trios and initial variant filtering

Two GIAB trios are considered:

| Trio | Child | Father | Mother |
|---|---|---|---|
| Ashkenazi Jewish | HG002 | HG003 | HG004 |
| Han Chinese | HG005 | HG006 | HG007 |

For small variants, each DeepVariant trio genotype
must satisfy **GQ ≥20** and **DP ≥10**. For SVs,
GRCh38-based INS and DEL calls are obtained from
Manta and Delly, respectively, and graph-based SV
calls from `vg call`. SVs are restricted to
autosomal `FILTER=PASS` INS/DEL **51–50,000 bp**
with **GQ ≥20**. Delly calls additionally require
precise breakpoints. Graph-based SV calls additionally
require **DP ≥10**; this cannot be applied to the
GRCh38 Manta/Delly callsets where DP is unavailable.

The three sample VCFs are normalized separately.
`bcftools merge` combines trio calls, preserving
called `0/0` genotypes. Missing genotypes remain
missing; they are not imputed to reference.

### 9.2 Observed inconsistency

RTG Tools `mendelian` identifies variants that
cannot be inherited from the called parental
genotypes. The annotated VCF's `INFO/MCV` marks
an observed Mendelian-inconsistent variant.
The analysis is **child-centered**, restricted
to autosomal child non-reference variants with
testable diploid genotypes for all three members.

\[
\mathrm{Observed\ inconsistency\ rate\ (\%)}
=\frac{\mathrm{Observed\ inconsistent\ child\ variants}}
{\mathrm{Testable\ child\ nonreference\ variants}}\times100
\]

This rate is a genotype inheritance incompatibility
measure, not proof of calling errors or of de novo
mutation.

### 9.3 Truth-set adjudication of inconsistent variants

Each observed inconsistent variant is checked
against the *child's* matching truth set within
its high-confidence region. Small variants are
compared with hap.py/RTG vcfeval; SVs are compared
using the same Truvari criteria in Section 7.

- **True inconsistent:** the observed child call
  is allele- and genotype-concordant with the truth.
- **False inconsistent:** the observed child call
  is unsupported by the truth or has a discordant
  genotype within an adjudicable high-confidence region.
- **Unresolved:** the truth comparison did not provide
  a conclusive classification.
- **Outside confidence:** the candidate is outside
  the corresponding truth-set confidence region.

\[
\mathrm{False\ inconsistent\ rate\ (\%)}
=\frac{\mathrm{False\ inconsistent\ variants}}
{\mathrm{True\ inconsistent}+\mathrm{False\ inconsistent}}\times100
\]

The denominator contains only adjudicated
inconsistent variants from the truth-confidence
intervals. An observed inconsistency outside those
regions is **not** counted as false. Because parental
genotypes are not independently re-verified against
their own truth sets, the terms true and false here
describe **truth-concordance of the child's call**,
not verified biological truth of the trio inheritance
event (Figure 4c; Table a2).

Figure 4c contains the SNP, small-indel, INS and DEL
results for HG002–HG004. HG005–HG007 results are
not included in that figure because no corresponding
SV truth set was available. GRCh38-based SV calls
are also excluded from that figure because their
breakpoints and ALT sequences are not directly
comparable to the graph-based calls.

---

## 10. References

1. Wagner J, Olson ND, Harris L, et al. (2022). Benchmarking challenging small variants with linked and long reads. *Cell Genomics*. https://doi.org/10.1016/j.xgen.2022.100128
2. English AC, Menon VK, Gibbs RA, Metcalf GA, Sedlazeck FJ. (2022). Truvari: refined structural variant comparison preserves allelic diversity. *Genome Biology*. https://doi.org/10.1186/s13059-022-02840-6
3. Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*. https://doi.org/10.1093/gigascience/giab008
4. Cleary JG, Braithwaite R, Gaastra K, et al. (2015). Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. *bioRxiv*. https://doi.org/10.1101/023754
5. Wu D, Yang C, Chen Q, et al. (2026). Hidden genetic diversity in 320 nearly-complete East Asian genome assemblies. *bioRxiv*. https://doi.org/10.64898/2026.07.27.740875
6. Suzuki Y, Owa C, Kobayashi H, et al. (2026). Accessing medically relevant complex regions with a pangenome graph of 20 near-complete Japanese haplotypes. *Nature Communications*. https://doi.org/10.1038/s41467-026-73461-x

---

## 11. Software and versions

| Tool | Version or requirement | Use |
|---|---|---|
| DeepVariant | 1.6 | Long-read small-variant truth generation; short-read baseline |
| bcftools | 1.24 (manuscript) | VCF normalization, genotype filtering and merging |
| hap.py | Record installed version | Small variant vcfeval benchmarking and annotated VCF |
| RTG Tools | 3.13 | `vcfeval` and `mendelian` |
| Truvari bench | **5.5.0** | SV benchmarking and allele/genotype concordance |
| Python | **≥3.12** for Truvari 5.5.0 | Pipeline helper execution |
| pysam | Record installed version | SV truth-set preparation and allele matching |
| samtools, bgzip, tabix | Record installed versions | FASTA and VCF indexing |

The pipeline reproduces the stated comparison criteria but
does not bundle the underlying GIAB and KOR truth VCFs,
confidence BEDs, personalized GBZ files or raw reads.
Those inputs must be supplied separately. For a reproducible
run, record the full software environment, the reference
build, sample list, truth-set versions and input checksums.
