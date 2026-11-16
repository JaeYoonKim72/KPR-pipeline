# Population_genetics_in_AMY_flanking_locus.sh

The human amylase (AMY) locus on chromosome 1 contains the salivary
AMY1 and pancreatic AMY2 gene families. Tandem and segmental
duplications, inversions, deletions, variable gene copy number and
extensive sequence homology have produced a diversity of structural
haplotypes that cannot be adequately represented by a single linear
reference genome. Bolognini et al. (2024) referred to the
copy-number-variable, structurally complex interval as the
**structurally variable region (SVR)**. Its arrangement and length
vary between haplotypes, even when they carry the same gene copy
numbers [1]. The KPR study resolved 31 structural architectures
among 124 KPR/HPRC/reference haplotype assemblies; this graph-based
structural analysis is documented separately in
[AMY_locus_analysis_in_pagenomegraph.md](AMY_locus_analysis_in_pagenomegraph.md).

The SVR presents several obstacles to conventional SNP-based
population genetics. Highly similar paralogous gene copies generate
ambiguous short-read mappings, copy-specific variant assignment and
genotyping uncertainty, while interlocus gene conversion can
exchange sequence between paralogs and obscure their genealogical
relationships [13,14]. Standard diploid-locus estimates of
nucleotide diversity, FST, linkage disequilibrium and
haplotype-based statistics require well-defined allelic
observations at homologous positions. More importantly, treating
every duplicated AMY copy as a separate independent person-like
sample constitutes **pseudoreplication**: the number of analyzed
copies is determined by the copy-number genotype, and multiple
copies from the same individual share a genetic, population and
pedigree context. An increased count of gene copies therefore does
not increase the number of independent sampled individuals or
justify the usual sampling variance and confidence intervals of
diploid population analyses. These problems cannot be resolved
simply by assigning additional sample IDs to individual copies.

An alternative is to characterize **single-copy, non-duplicated
sequences flanking the SVR**. Bolognini et al. identified the unique
distal and proximal regions as principal bundles 0 and 1, in which
SNPs can be genotyped reliably, and used their haplotype sequences
to investigate the genealogical background of the intervening
structural architectures [1]. They observed substantial
cross-SVR linkage disequilibrium, approximately 7- to 20-fold
higher than that for similarly separated SNP pairs elsewhere on
chromosome 1, together with reduced recombination in the region.
Consequently, flanking haplotypes can retain informative
associations with linked SVR structures and their demographic
and evolutionary history, although recombination and recurrent
structural mutations mean that a flank does not uniquely identify
an AMY copy-number state.

Here, we examine the distal and proximal AMY flanks at the
**locus level**, rather than applying diploid SNP estimators to
individual AMY paralogs. Within-population diversity (pi,
Watterson's theta and Tajima's D), allele-frequency spectra,
within- and across-flank LD, and nSL describe sequence diversity
and linked haplotype structure. Between-population Hudson FST and
XP-nSL evaluate differentiation and relative haplotype
homozygosity. Each observed statistic is compared, where
implemented, with an appropriate empirical genomic background,
including accessible-base and physical-distance matching.
Patterns consistent with linked selection, background selection,
drift, demographic history or reduced recombination are
interpreted as **indirect signatures**; none of these flanking
statistics alone proves selection on an AMY gene copy or
establishes a particular causal SVR architecture [1,13,14].

Analyses use GRCh38-aligned, diploid SNP genotypes from
consent-eligible KOR302 individuals and public 1000 Genomes
populations. The earlier 818-person cohort and all derived
inputs are excluded. The code supports configurable cohort
lists, focal intervals, optional uncertainty analyses, and
independently generated genomic backgrounds. Status and
implementation limits are summarized below.

**Implementation status.** The public `all` workflow calculates
focal pi, projected Watterson's theta, an exploratory projected
Tajima's D, within-flank LD, cross-flank LD, Hudson FST, and,
when selscan v3.0.1 is available, raw nSL/XP-nSL scores.
Individual-level bootstrap/jackknife is implemented for
SNP-based statistics; genomewide window/block-bootstrap
comparisons are implemented for pi, theta, Tajima's D and FST.
Normalized haplotype metrics, genomewide LD enrichment and
the complete manuscript uncertainty analyses require the
separate inputs or downstream procedures documented below.
Do not treat these public exploratory outputs as automatically
identical to the verified manuscript tables.

---

## 1. Set up

### 1.1 Software

Install the programs under [Software](#software). The non-haplotype statistics
require `bash`, `python3` and `bcftools`. `selscan v3.0.1` is required
**only** for nSL and XP-nSL.

Create and **activate** a mamba/conda environment before running:

~~~bash
mamba create -n amy_popgen -c conda-forge -c bioconda python bcftools bedtools htslib selscan=3.0.1
conda activate amy_popgen
selscan --version
~~~

If that exact version cannot be resolved on the selected channel/platform,
install the v3.0.1 binary separately and make it available on `PATH`.
The wrapper never installs software automatically. If v3.0.1 is missing,
`all` logs warnings and records `SKIPPED_MISSING_SELSCAN` for the two
haplotype analyses; all applicable SNP-based calculations continue.

### 1.2 AMY locus and Bolognini flanking coordinates

The structurally variable region (SVR) includes highly similar AMY1 and
AMY2 gene copies, frequent copy-number variation and gene conversion.
Copy-specific short-read population genetics is not treated as independent
diploid-locus variation. Instead, the unique flanking regions are used to
characterize variation and long-range LD at the **AMY locus level**.

Following **Bolognini et al. (2024)** [1], who examined 94 HPRC
haplotype-resolved assemblies and their principal AMY graph bundles,
we retain the exact GRCh38 definitions of **bundle 0 (distal)**
and the first 66 kb of **bundle 1 (bundle 1a, proximal)**:

| Region | GRCh38, 1-based inclusive | Length |
|---|---|---:|
| Distal, bundle 0 | `chr1:103,456,163–103,561,526` | 105,364 bp |
| Excluded interval (contains SVR) | `chr1:103,561,527–103,760,697` | 199,171 bp |
| Proximal, bundle 1a | `chr1:103,760,698–103,826,698` | 66,001 bp |
| Entire AMY exclusion span | `chr1:103,456,163–103,826,698` | 370,536 bp |

These are **default parameters, not hardcoded/immutable values**:
override them using `--distal` and `--proximal` (or `DISTAL` and
`PROXIMAL`). Do not include the intervening SVR in SNP-based
diversity, LD, FST or haplotype estimates.

### 1.3 GRCh38 StrictMask and non-N callable regions

Download the 20160622 GRCh38 StrictMask BED from the 1000 Genomes
Project:

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed

The file named **GRCh38.StrictMask.nonN.bed.gz** is generated locally,
not downloaded directly. It contains the reference-accessible
StrictMask intervals after excluding N bases in the same GRCh38 FASTA
used for variant calling. Both source and derived files use zero-based,
half-open BED coordinates. Install bedtools, bgzip and tabix in addition
to the programs listed in section 1.1.

~~~bash
set -euo pipefail

GRCH38_FASTA=/path/to/GRCh38.fa
MASK_DIR=inputs/masks
STRICT_URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed

mkdir -p "$MASK_DIR"
curl -fL --retry 3 "$STRICT_URL" -o "$MASK_DIR/20160622.allChr.mask.bed"

python3 - "$GRCH38_FASTA" "$MASK_DIR/GRCh38.Nruns.bed" <<'PY'
import gzip
import re
import sys

fasta_path, out_path = sys.argv[1:]
open_fasta = gzip.open if fasta_path.endswith(".gz") else open

with open_fasta(fasta_path, "rt") as fin, open(out_path, "w") as fout:
    chrom = None
    offset = 0
    for line in fin:
        if line.startswith(">"):
            chrom = line[1:].split()[0]
            offset = 0
            continue
        seq = line.strip()
        for match in re.finditer(r"[Nn]+", seq):
            fout.write(
                f"{chrom}\t{offset + match.start()}\t{offset + match.end()}\n"
            )
        offset += len(seq)
PY

LC_ALL=C sort -k1,1 -k2,2n "$MASK_DIR/GRCh38.Nruns.bed" |
  bedtools merge -i - > "$MASK_DIR/GRCh38.Nruns.merged.bed"

awk 'BEGIN{OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1,$2,$3}' \
  "$MASK_DIR/20160622.allChr.mask.bed" |
  LC_ALL=C sort -k1,1 -k2,2n > "$MASK_DIR/StrictMask.sorted.bed"

bedtools subtract \
  -a "$MASK_DIR/StrictMask.sorted.bed" \
  -b "$MASK_DIR/GRCh38.Nruns.merged.bed" |
  LC_ALL=C sort -k1,1 -k2,2n |
  bedtools merge -i - |
  bgzip -c > "$MASK_DIR/GRCh38.StrictMask.nonN.bed.gz"

tabix -f -p bed "$MASK_DIR/GRCh38.StrictMask.nonN.bed.gz"
tabix -l "$MASK_DIR/GRCh38.StrictMask.nonN.bed.gz" | grep -qx chr1
~~~

The reference FASTA, accessibility mask and VCF must use matching
contig names and the same GRCh38 assembly. This is a reference-level
accessibility mask, **not** a sample-specific genotype-callability
mask. Keep both the downloaded mask and derived file for provenance.

### 1.4 Input SNP data and population lists

`VCF` must be a **coordinate-sorted, indexed joint multi-sample VCF**
containing phased **diploid biallelic SNPs**, aligned to the same GRCh38
REF/ALT representation. Provide one text file of sample IDs per population,
one name per line. Each sample belongs to exactly one population; the
**first population** is the between-population query. Analyzing two populations
does not require equal sample counts.

Genotype-level DP≥10 and GQ≥20, sample-level consent/QC, and SHAPEIT5
phasing/IMPUTE5 imputation are **upstream inputs**, not operations carried
out by this script. The focal Python statistic engine additionally enforces
genotype call rate ≥90% **for each relevant population**. All SNPs used in
π/θW/Tajima's D and FST retain low-frequency variants (no MAF filter).
LD requires MAF≥0.05; nSL/XP-nSL core MAF is independently configurable.

For callable-base denominators, provide `CALLABLE_BED` from GRCh38
StrictMask (20160622) [2] **intersected with GRCh38 non-N regions**.
This file is BED0 (0-based half-open). No callable BED produces
a loud warning and uses the **physical length** as a denominator;
this is not equivalent to the paper's masked analysis.

In our revision the focal population was KOR302; the eight public
comparison populations are CHN, CDX, JPT, KHV, SAS, EUR, AMR and AFR.
The historic 818-person dataset and all 818-derived files are excluded.

### 1.5 Input layout

~~~text
inputs/
├── joint.phased.biallelic.snps.vcf.gz
├── joint.phased.biallelic.snps.vcf.gz.tbi   (or .csi)
├── masks/GRCh38.StrictMask.nonN.bed.gz
└── sample_lists/
    ├── KOR.txt
    ├── CHN.txt
    ├── JPT.txt
    └── ...
~~~

For `GENOMIC_BOOTSTRAP=1`, `GENOME_VCF` must contain the selected
chromosomes, not merely the AMY interval. Default `GENOME_VCF=$VCF`.

### 1.6 Environment

~~~bash
export VCF=/path/to/joint.phased.biallelic.snps.vcf.gz
export POPS="KOR=/path/KOR.txt CHN=/path/CHN.txt JPT=/path/JPT.txt"
export CALLABLE_BED=/path/GRCh38.StrictMask.nonN.bed.gz
export OUT_DIR=/path/results/AMY_popgen
export THREADS=4
~~~

---

## 2. Run

### 2.1 Quick start

~~~bash
bash Population_genetics_in_AMY_flanking_locus.sh --help
bash Population_genetics_in_AMY_flanking_locus.sh all
~~~

The same inputs can be supplied as command-line arguments:

~~~bash
bash Population_genetics_in_AMY_flanking_locus.sh all \
  --vcf /path/to/joint.phased.biallelic.snps.vcf.gz \
  --pop KOR=/path/KOR.txt --pop CHN=/path/CHN.txt \
  --mask /path/GRCh38.StrictMask.nonN.bed.gz \
  --out /path/results
~~~

A **single population** computes within-population diversity, within-flank
and across-flank LD and nSL. With **two or more populations**, these are
computed per population and Hudson FST/XP-nSL are added for the first
population versus each subsequent population.

### 2.2 Optional individual bootstrap and jackknife (default OFF)

~~~bash
INDIVIDUAL_CI=1 \
IND_BOOT_REPS=100 \
IND_JACKKNIFE_FRACTION=0.30 \
bash Population_genetics_in_AMY_flanking_locus.sh all
~~~

Or append `--individual-ci --bootstrap-reps 100 --jackknife-fraction 0.30`.
For the **paper's full-sampling configuration**, use
`IND_BOOT_REPS=10000 IND_JACKKNIFE_FRACTION=1.0`.

### 2.3 Optional genomewide analysis (default OFF)

Fast chr1 pilot: **50-Mb blocks**, 1-kb sliding windows and
100 bootstrap replicates:

~~~bash
GENOMIC_BOOTSTRAP=1 \
GENOMIC_CHROMS=chr1 \
GENOMIC_BLOCK_MB=50 \
GENOMIC_WINDOW_STEP=1000 \
GENOMIC_BOOT_REPS=100 \
bash Population_genetics_in_AMY_flanking_locus.sh all
~~~

The paper's **full autosomal chromosome-stratified 5-Mb analysis**
used autosomes 1–22 and **10,000** block bootstrap replicates:

~~~bash
export GENOMIC_CHROMS="$(printf 'chr%s,' {1..22} | sed 's/,$//')"
export GENOMIC_BLOCK_MB=5
export GENOMIC_BOOT_REPS=10000
export GENOMIC_WINDOW_STEP=1000
export GENOME_VCF=/path/to/all_autosomes_joint.phased.biallelic.snps.vcf.gz

GENOMIC_BOOTSTRAP=1 \
bash Population_genetics_in_AMY_flanking_locus.sh all
~~~

These genomewide settings correspond to the **published analysis design**;
the simplified public code's currently supported statistics are explicitly
listed under `9.6`.

---

## 3. Steps

| Step | Primary statistic | Required population count | Main result |
|---|---|---|---|
| `diversity` | π, projected θW, projected Tajima's D | ≥1 | `stats/diversity.tsv` |
| `ld` | Mean within-flank corrected dosage r² | ≥1 | `stats/within_ld.tsv` |
| `across_ld` | Mean distal × proximal corrected dosage r² | ≥1 | `stats/across_ld.tsv` |
| `fst` | Hudson FST | ≥2 | `stats/hudson_fst.tsv` |
| `nsl` | Unphased nSL raw core scores | ≥1, selscan | `selscan/` |
| `xpnsl` | Phased XP-nSL raw core scores | ≥2, selscan | `selscan/` |
| `individual_ci` | Bootstrap and jackknife for SNP-derived focal statistics | ≥1 | `individual_ci/individual_ci.tsv` |
| `genomic` | 1-kb genomic background windows, block-bootstrap genomic percentile | ≥1; genomewide VCF | `genomic/genomic_percentile_ci.tsv` |
| `all` | Above focal steps; optional CI and genomic if ON | ≥1 | `status.tsv` |

`all` ALWAYS attempts focal within and eligible between statistics.
`INDIVIDUAL_CI` and `GENOMIC_BOOTSTRAP` are OFF by default, even
when the run mode is `all`. Use the corresponding environment variables
or switches to enable them. `nsl`/`xpnsl` are skipped with a warning if
selscan v3.0.1 is absent. No scheduler or automatic parallel job launch
is required.

---

## 4. Configuration

| Variable / option | Default | Description |
|---|---|---|
| `VCF` / `--vcf` | required | Joint phased diploid bi-allelic SNP VCF and index |
| `POPS` / repeated `--pop` | required | `NAME=sample-list`; first group is query |
| `OUT_DIR` / `--out` | `./results/amy_flanking_population_genetics` | Output root |
| `DISTAL` / `--distal` | `chr1:103456163-103561526` | Distal bundle 0 |
| `PROXIMAL` / `--proximal` | `chr1:103760698-103826698` | Proximal bundle 1a |
| `CALLABLE_BED` / `--mask` | unset | GRCh38 StrictMask ∩ non-N callable BED |
| `PROJECTION_N` / `--projection` | `168` | Projection allele copies, n₀ |
| `CALL_RATE` / `--call-rate` | `0.90` | Per-population genotype call rate |
| `LD_MAF` / `--ld-maf` | `0.05` | MAF for LD SNPs |
| `LD_DISTAL_MAX` / `--ld-distal-max` | `100000` | Distal within-LD distance must be less than value |
| `LD_PROXIMAL_MAX` / `--ld-proximal-max` | `60000` | Proximal within-LD distance must be less than value |
| `NSL_CORE_MAF` / `--core-maf` | `0.05` | Core SNP MAF cutoff for selscan |
| `NSL_MAX_EXTEND` / `--max-extend-nsl` | `50` | Maximum SNP markers for nSL/XP-nSL extension |
| `THREADS` | `4` | selscan threads (no outer parallelism) |
| `INDIVIDUAL_CI` / `--individual-ci` | `0` (OFF) | Enable individual-level resampling CIs |
| `IND_BOOT_REPS` / `--bootstrap-reps` | `100` | Individual bootstrap replicates; paper used 10,000 |
| `IND_JACKKNIFE_FRACTION` / `--jackknife-fraction` | `0.30` | ceil(30% N) deleted-person replicates/group; paper used 1.0 |
| `GENOMIC_BOOTSTRAP` / `--genomic` | `0` (OFF) | Enable genomewide windows and block bootstrap |
| `GENOMIC_CHROMS` / `--chromosomes` | `chr1` | Pilot one chromosome; paper used chr1–chr22 |
| `GENOMIC_BLOCK_MB` / `--block-mb` | `50` | Pilot physical block; paper used 5 Mb |
| `GENOMIC_WINDOW_STEP` / `--window-step` | `1000` | Sliding-window step in bp; paper also used 1 kb |
| `GENOMIC_BOOT_REPS` / `--genomic-reps` | `100` | Genome block draws; paper used 10,000 |
| `GENOMIC_FST_CALLABLE_MATCH` / `--fst-callable-match` | `0.10` | FST background callable-bp ±10%; -1 disables |
| `GENOME_VCF` / `--genome-vcf` | `VCF` | Genomewide joint VCF for configured chromosomes |
| `GENOME_CHROM_SIZES` / `--chrom-sizes` | built-in GRCh38 | Optional reference `.fai` to override chromosome lengths |
| `SEED` / `--seed` | `20261117` | Deterministic resampling seed |

Important: `IND_JACKKNIFE_FRACTION=0.30` is a
**partial delete-one jackknife approximation**, not the exact complete
leave-one-out procedure. To reproduce the original method set it to
`1.0`.

---

## 5. Output

~~~text
results/amy_flanking_population_genetics/
├── stats/
│   ├── diversity.tsv             π, θW, exploratory Tajima's D
│   ├── within_ld.tsv             within-flank corrected r²
│   ├── across_ld.tsv             distal × proximal corrected r²
│   └── hudson_fst.tsv            first-vs-other Hudson FST
├── individual_ci/                only when INDIVIDUAL_CI=1
│   ├── individual_ci.tsv         percentile/shifted bootstrap and JK intervals
│   └── individual_replicates.tsv person-resampling replicate values
├── genomic/                      only when GENOMIC_BOOTSTRAP=1
│   └── genomic_percentile_ci.tsv upper-tail ranks and block bootstrap 95% interval
├── intermediate/                 region VCFs used for selscan
├── selscan/                      raw nSL and XP-nSL outputs (if installed)
├── logs/
└── status.tsv                    DONE, RAW_ONLY, SKIPPED, FAILED
~~~

Correspondence with paper supplementary tables:

| Table | Within/between | Current public-code focal coverage | Additional published metrics |
|---|---|---|---|
| X1 | Within | π, θW, exploratory Tajima D, raw nSL | Validated projected Tajima D, normalized nSL; all genomic ranks/CIs |
| X2 | Within | Within-LD, across-LD corrected r² | Distance-matched LD enrichment and genomewide LD rank/interval |
| X3 | Between | Hudson FST, raw XP-nSL | Normalized XP-nSL, empirical 93:93 subsampling, complete ranks/CIs |

`status.tsv` distinguishes successfully computed values from raw,
skipped, or failed stages. **Raw nSL and XP-nSL are not the signed/absolute
normalized scores required for interpretation in the final paper.**

---

## 6. Within-population statistics

### 6.1 Nucleotide diversity (π)

For a biallelic SNP `i` with alternative allele count `AC_i` and
number of called alleles `AN_i`, the sample-size-unbiased pairwise
nucleotide difference estimator is [3]:

`π_i = 2 AC_i(AN_i − AC_i) / [AN_i(AN_i − 1)]`

The per-base nucleotide diversity of a flank is:

`π = (Σ_i π_i) / L_callable`

`L_callable` is the number of StrictMask ∩ non-N bases in the flank.
SNPs with `0 < AC_i < AN_i` are included without MAF filtering.
The formula uses site-wise observed allele counts, without applying
frequency projection to π.

### 6.2 Watterson's theta (θW)

To compare populations with unequal sample sizes, the
segregating-site frequency spectrum is projected to
`n₀ = PROJECTION_N` allele copies (default 168).
At site `i` the expected projected segregating probability is:

`Pr(segregating after projection) =
1 − C(AC_i,n₀)/C(AN_i,n₀) − C(AN_i−AC_i,n₀)/C(AN_i,n₀)`.

The expected projected count is `E[S_n₀] = Σ_i Pr_i` and
Watterson's estimator [4] is:

`θW = E[S_n₀] / (a_(n₀−1) × L_callable)`

where `a_(n₀−1) = Σ_(j=1)^(n₀−1) (1/j)`.
The manuscript's choice `n₀=168` corresponds to 90% genotype
call rate among the smallest population CDX (n=93).

### 6.3 Tajima's D

Tajima's D [5] compares total pairwise differences with the
segregating-site estimator and scales their difference by
the expected variance:

`D = [π_sum − S/a₁] / sqrt(e₁S + e₂S(S−1))`

where `a₁`, `e₁` and `e₂` depend on sample size.
**Important implementation note:** the current helper substitutes
projected `E[S_n₀]` for `S` and uses an unprojected π sum. Its output
is explicitly exploratory (`Tajima_D` in `stats/diversity.tsv`), not
a calibrated demographic neutrality-test P value. Publication values
should be checked against the validated projected-frequency pipeline.

### 6.4 Within-flank linkage disequilibrium (LD)

For unphased diploid genotype dosages `0,1,2`, the squared Pearson
correlation of SNPs `i` and `j` is the raw LD statistic
`r²_raw = Corr(g_i,g_j)²`. With `n` individuals having valid
genotypes at both sites, the sample-size-corrected estimate is:

`r²_corrected = [(n−1) r²_raw − 1] / (n−2)`.

Corrected negative values are retained. Use biallelic SNPs with
`MAF ≥ 0.05` and call rate ≥90%. Pairs must have physical
separation **<100 kb** within distal or **<60 kb** within proximal.
The summary is the **arithmetic mean of corrected r² across valid
SNP pairs**, not the square of the mean r.

### 6.5 Cross-flank / across-locus LD

Across-locus LD is estimated by pairing eligible SNPs from the distal
and proximal unique flanks. The focal output `stats/across_ld.tsv`
contains the arithmetic mean of sample-size-corrected dosage r² and
its valid pair count for each population. A focal r² alone does not
quantify excess LD relative to the genome: enrichment requires a
separately estimated, distance-matched background.

For a consistent estimator and pair inclusion rule:

- `F` is the mean focal cross-flank r².
- `B_k` is the mean genomic r² for SNP pairs separated by distance
  bin `k` (1-kb bins).
- `w_k` is the fraction of eligible focal SNP pairs in bin `k`.
- `B = sum_k(w_k * B_k)` is the focal-distance-weighted background.
- `enrichment = F / B` and `delta = F - B`.

For example, F = 0.12 and B = 0.03 would give a fourfold
enrichment and a delta of 0.09. Both F and B must use identical
population, QC, MAF, callable-region and r² definitions. A
corrected-r² focal estimate cannot be divided by a raw-r² background.

**Distance-filter distinction.** This public helper uses all
eligible distal-by-proximal pairs and does not impose the
original study's additional 190–370-kb distance band. To
reproduce the completed cross-SVR analysis, impose that same
190–370-kb band on the focal and genomic pairs. The original
background used chromosome 1 distance-matched pairs, with
weights from the focal 1-kb pair-distance histogram. It is not
an all-autosome rank. Section 9.2 explains how to use the
previously derived results and compute a corrected-r² ratio.

### 6.6 nSL: within-population haplotype homozygosity

Unphased nSL is calculated with selscan v3.0.1 [6–8,12].
The core must meet MAF >= 0.05 by default; low-frequency
markers may contribute to extension under `--keep-low-freq`.
The focal extension stops at the flank boundary and uses
`--max-extend-nsl 50` by default.

#### Normalize the nSL results

Normalize each population independently with the genomewide
**unphased** nSL reference for that population. Our
normalization reference uses chromosomes 2–22 and 100
core-allele-frequency bins. Each core is standardized as
`z_i = (raw_nSL_i - mean_bin) / sd_bin`. The background
must be computed with consistent core-MAF, phasing mode,
QC filters and extension cap; the two short AMY focal files
must not be used to estimate the reference distribution.

First generate raw unphased nSL scores across chr2–22
for each population, using `selscan --nsl --unphased`
and the same options as the focal run. With these
independently derived background files in place:

~~~bash
set -euo pipefail
export OUT_DIR=/path/to/results/amy_flanking
export NSL_BG=/path/to/autosomal_unphased_nsl/KOR/chr2_22
mkdir -p "$OUT_DIR/normalization"

mapfile -t bg < <(find "$NSL_BG" -type f -name '*.nsl.out' | sort)
mapfile -t focal < <(find "$OUT_DIR/selscan" -maxdepth 1   -type f -name 'KOR.*.nsl.out' | sort)
(( ${#bg[@]} > 0 && ${#focal[@]} == 2 )) || {
  echo "Expected background files and both AMY flank files" >&2
  exit 1
}

selscan norm --nsl --bins 100   --log "$OUT_DIR/normalization/KOR.nsl.chr2_22.log"   --files "${bg[@]}"

selscan norm --nsl --bins 100   --log-input "$OUT_DIR/normalization/KOR.nsl.chr2_22.log"   --files "${focal[@]}"
~~~

The normalized `.out.norm` results can then be summarized
separately for distal and proximal flanks as mean absolute
nSL and the proportion of valid cores with absolute
normalized nSL > 2. Repeat for CHN, CDX, JPT, KHV,
SAS, EUR, AMR and AFR with each population's own
reference log. Do not pool populations. The current
public wrapper generates only raw focal nSL scores;
genomewide reference scans, normalization and ranking
remain independent downstream steps.

---

## 7. Between-population statistics

### 7.1 Hudson FST

For KOR (query) vs each reference population, the Hudson estimator
[9,10] corrects allele-frequency sampling variability at differing
numbers of called alleles. Using `p₁=AC₁/AN₁` and `p₂=AC₂/AN₂`:

`N_i = (p₁−p₂)² − p₁(1−p₁)/(AN₁−1) − p₂(1−p₂)/(AN₂−1)`

`D_i = p₁(1−p₂) + p₂(1−p₁)`

`FST_Hudson = (Σ_i N_i) / (Σ_i D_i)`.

Only shared POS/REF/ALT SNPs meeting genotype call rate ≥90%
in **both** populations and callable-mask requirements are used.
No additional MAF filter is applied. A negative estimated FST is
retained rather than truncated to zero.

The KOR302 paper's focal shared SNP counts were 166 (distal) and
128 (proximal). These numbers are dataset-dependent and are
**not hardcoded** into the script.

### 7.2 XP-nSL

Phased XP-nSL compares haplotype-homozygosity extension
between KOR (query) and each reference population using
matched POS/REF/ALT variants and selscan v3.0.1 [7,8,12].
The core-MAF default is 0.05 and extension is capped at
50 SNP markers. Positive signed XP-nSL indicates longer
haplotype homozygosity in KOR, while a negative sign
indicates longer homozygosity in the reference population.

#### Normalize the XP-nSL results

Unlike nSL, XP-nSL is standardized by the genomewide
mean and SD **without allele-frequency bins**:
`z_i = (raw_XPnSL_i - mean_genome) / sd_genome`.
Fit the normalization reference separately for each
ordered population pair, using phased genomewide
autosomal XP-nSL scores with consistent QC, sites
and extension settings. Exclude AMY from fitting the
reference distribution to avoid focal self-inclusion.

First generate genomewide raw XP-nSL files on
chr1–22 for KOR versus each reference population,
using `selscan --xpnsl` with the query/reference
order retained. Then run for each pair:

~~~bash
set -euo pipefail
export OUT_DIR=/path/to/results/amy_flanking
export PAIR=KOR_vs_CHN
export XPN_BG=/path/to/autosomal_phased_xpnsl/KOR_vs_CHN/chr1_22
mkdir -p "$OUT_DIR/normalization"

mapfile -t bg < <(find "$XPN_BG" -type f -name '*.xpnsl.out' | sort)
mapfile -t focal < <(find "$OUT_DIR/selscan" -maxdepth 1   -type f -name "${PAIR}.*.xpnsl.out" | sort)
(( ${#bg[@]} > 0 && ${#focal[@]} == 2 )) || {
  echo "Expected background files and both AMY flank files" >&2
  exit 1
}

selscan norm --xpnsl   --log "$OUT_DIR/normalization/${PAIR}.xpnsl.autosomes.log"   --files "${bg[@]}"

selscan norm --xpnsl   --log-input "$OUT_DIR/normalization/${PAIR}.xpnsl.autosomes.log"   --files "${focal[@]}"
~~~

Repeat with `PAIR=KOR_vs_CDX`, `KOR_vs_JPT`,
`KOR_vs_KHV`, `KOR_vs_SAS`, `KOR_vs_EUR`,
`KOR_vs_AMR` and `KOR_vs_AFR`, changing
`XPN_BG` to the matching population-pair directory.
Never reuse the normalization log from a different pair
or switch query/reference orientation.

The final flank-level summaries use valid normalized
cores: scored-core count, mean signed XP-nSL,
mean absolute XP-nSL, fraction > 2 and fraction < -2.
The **mean absolute** normalized score is used for
the upper-tail genomic percentile, whereas the
signed mean preserves directional information.

The full study also evaluated 100 independent,
without-replacement 93:93 individual subsamples.
CDX has exactly 93 individuals, so its entire
cohort occurs in each replicate; only the KOR
subset changes. These descriptive subsampling
intervals are not implemented in the public helper,
which currently produces raw focal XP-nSL only.

---

## 8. Individual bootstrap and jackknife confidence intervals

### 8.1 Individual bootstrap

For each population, resample **individuals**, with replacement,
to the original population size `N`. In Hudson FST comparisons,
resample the query and reference **independently**, preserving
their respective sizes. Recompute the statistic in each replicate
(`B=IND_BOOT_REPS`).

For replicate values `t₁,...,t_B` and observed `t_obs`,
the raw 95% percentile interval is:

`CI_boot_raw = [Q_0.025(t*), Q_0.975(t*)]`.

The mean shift is `b = mean(t*) − t_obs` and the
bias-recentered interval is:

`CI_boot_shifted = [Q_0.025(t*)−b, Q_0.975(t*)−b]`.

Both are written to `individual_ci/individual_ci.tsv`.
For the manuscript, `B=10,000`, and the mean-shifted FST
interval was the preferred bootstrap presentation.
The public **default when ON** is `B=100` only.

### 8.2 Delete-one jackknife

For each group, recompute the statistic after removing one
individual. Complete leave-one-out means `N` deleted-person
replicates in a single population, and `N_query + N_reference`
recalculations for a two-population FST comparison. The
standard error is based on between-replicate variability
(with the contribution of each group calculated separately
for between-population statistics). The conventional interval is:

`CI_JK = t_obs ± 1.96 × SE_JK`.

**Pilot default:** `IND_JACKKNIFE_FRACTION=0.30` deletes
`ceil(0.30×N)` randomly selected individuals **one at a time**.
This is not an exact classical jackknife. The helper uses
a finite-population variance estimate from these randomly
selected leave-one-out values, labels the result
`SUBSET_JACKKNIFE_APPROX` and reports the number used/total
per group. Setting `IND_JACKKNIFE_FRACTION=1.0` runs the
**complete** leave-one-out jackknife, as used in the paper.

### 8.3 Implemented statistics and limitations

Optional individual CI currently supports **π, projected θW,
exploratory Tajima's D, within-LD, across-LD and Hudson FST**.
It genuinely resamples individual genotype dosages and
recalculates each statistic. The nSL/XP-nSL CIs are **not**
produced because they require repeated haplotype-level selscan
runs with the original genomewide normalization; raw nSL scores
must not be treated as equivalent.

---

## 9. Genome-wide comparison and block bootstrap

### 9.1 Matched-size sliding-window backgrounds

For each flank independently, take windows of the **same physical
length** as the focal flank: distal `105,364 bp` and proximal
`66,001 bp`, with a **1,000-bp (1 kb) step** along the chosen
background chromosome(s). The entire locus span
`chr1:103,456,163–103,826,698` and every background window
overlapping this span are excluded.

Windows lacking callable bases or a valid denominator for a
statistic are omitted. Hudson FST background calculations
can additionally require callable-base length within ±10%
of the corresponding focal flank (`GENOMIC_FST_CALLABLE_MATCH=0.10`).
This is the callable-matched primary background defined for
our manuscript; unmatched windows are a sensitivity analysis.
For the physical chr1 pilot, the default still uses the callable
matching unless disabled.

The study's distance-matched genomic cross-flank LD background
preserves the two flank lengths and their relative positions,
and matches SNP-pair physical distances in 1-kb bins.
That LD background/enrichment calculation is not yet present
in this public helper.

### 9.2 Genome-wide LD enrichment

**Within-flank LD.** For each population, create an autosomal
reference distribution of sliding-window mean corrected
r², with the same physical length as the distal
(105,364 bp) or proximal (66,001 bp) flank, a
1-kb step, identical call-rate and MAF >= 0.05
filters, and the same maximum pair separation
(<100 kb distal; <60 kb proximal). Exclude windows
overlapping chr1:103456163–103826698. Compare
the focal within-LD mean to the distribution of
eligible window means; the genomic median can
serve as the denominator of an exploratory fold
enrichment, and the full distribution supports
the percentile in section 9.3. If the denominator
is not positive, report delta and percentile
rather than an unstable fold ratio.

**Cross-flank LD.** A valid genomic comparison
must account for the dependence of r² on
physical SNP-pair separation. Within each
population, partition focal SNP pairs into
1-kb distance bins. Let `w_k` be the proportion
of focal pairs in bin k, and let `B_k` be the
mean genomic r² in that distance bin, after
applying the same SNP, callable, sample and
r² filters. The matched reference is
`B = sum_k(w_k * B_k)`. Report the
observed mean `F`, the ratio `F/B`
and the difference `F-B`. This preserves
the observed focal pair-distance distribution.
A geometry-matched background may additionally
preserve both flank lengths and their relative
physical displacement; such analyses should
be labeled separately from the earlier
chromosome 1 pair-distance-matched analysis.

**Completed KOR example: original raw-r²
cross-SVR analysis.** The focal mean raw
r² was **0.469826650626** and the
chromosome 1 distance-matched mean raw r²
was **0.0198553931296**, yielding
**23.662419956-fold enrichment** and
an observed-minus-background difference
of **0.449971257496**. Those estimates
used 1-kb pair-distance-bin matching and
the original 190–370-kb restriction.
They do not represent all-autosome LD
enrichment and cannot be applied to
the corrected-r² focal file without
recomputing the corrected-r² background.

**How to use the public outputs.**
The file `stats/across_ld.tsv` contains
`population`, `n_pairs`, and
`mean_corrected_r2`. When an
independently generated distance-matched
genomic reference table contains
`population` and `matched_mean_corrected_r2`,
the two completed estimates can be joined:

~~~bash
python3 - "$OUT_DIR/stats/across_ld.tsv"   /path/to/matched_cross_flank_background.tsv <<'PY'
import csv
import math
import sys

focal_path, background_path = sys.argv[1:]
with open(background_path, newline="") as handle:
    background = {
        row["population"]: float(row["matched_mean_corrected_r2"])
        for row in csv.DictReader(handle, delimiter="\t")
    }

print("population\tfocal_r2\tmatched_r2\tfold\tdelta")
with open(focal_path, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        population = row["population"]
        focal = float(row["mean_corrected_r2"])
        reference = background[population]
        if not (math.isfinite(focal) and math.isfinite(reference)):
            raise ValueError("Invalid LD estimate for " + population)
        ratio = f"{focal/reference:.6g}" if reference > 0 else "NA"
        print(f"{population}\t{focal:.9g}\t{reference:.9g}\t"
              f"{ratio}\t{focal-reference:.9g}")
PY
~~~

This command joins **already computed** focal and
matched-genomic estimates; it does not perform genomic
LD estimation or bootstrapping. To estimate uncertainty
for an LD ratio by individual resampling, recompute
both numerator and denominator within each replicate
before dividing; do not divide confidence interval
endpoints. Genomic block bootstrap of the LD
percentile resamples background physical blocks and
is a different descriptive interval. Retain
availability flags for BCa and other CI methods:
the existing KOR raw ratio percentile interval
[19.14355, 21.69726] does not contain its
observed point estimate 23.66242 and must **not**
be described as a calibrated 95% confidence interval.

The completed cross-SVR statistics may be shown
alongside Table X2 as descriptive enrichment.
The public helper currently implements neither
genomic LD background generation nor LD
percentile/interval estimation; those modules
remain separate from `amy_flank_genomic_background.py`.

### 9.3 Upper-tail genomic percentile

For focal estimate `T_focal` and `W` eligible background windows:

`P_upper = 100 × #{w : T_w ≥ T_focal} / W`.

Thus a **smaller percentage means a more extreme upper-tail
value**. These ranks describe relative genomic positions; they
are **not selection P values**. For normalized nSL/XP-nSL,
the manuscript used the mean absolute normalized score to
calculate upper-tail rank (not the signed mean).

### 9.4 Chromosome-stratified physical block bootstrap

Split each chosen chromosome into contiguous, nonoverlapping
physical blocks of length `GENOMIC_BLOCK_MB`. Within each
chromosome, resample its original number of physical blocks
**with replacement**. Pool all eligible background windows
in the selected blocks and calculate the percentile at each
replicate. From `R=GENOMIC_BOOT_REPS` replicate percentiles,
report the mean, SD and 2.5/97.5 percentiles as a descriptive
95% interval. Report the number of physical and occupied
blocks, valid replicates and contributing windows.

Defaults **when ON**: chr1 only, **50 Mb**, **100**
block-bootstrap replicates, 1-kb window step.
**Manuscript configuration**: autosomes chr1–chr22, **5 Mb**,
**10,000** block-bootstrap replicates, same 1-kb window step.
These are settings, not distinct definitions of the estimator.

### 9.5 Example configuration

~~~bash
# Pilot
GENOMIC_BOOTSTRAP=1 GENOMIC_CHROMS=chr1 GENOMIC_BLOCK_MB=50 \
GENOMIC_BOOT_REPS=100 GENOMIC_WINDOW_STEP=1000 \
bash Population_genetics_in_AMY_flanking_locus.sh all

# Paper
GENOMIC_CHROMS="$(printf 'chr%s,' {1..22} | sed 's/,$//')" \
GENOMIC_BLOCK_MB=5 GENOMIC_BOOT_REPS=10000 \
GENOMIC_WINDOW_STEP=1000 GENOMIC_BOOTSTRAP=1 \
bash Population_genetics_in_AMY_flanking_locus.sh all
~~~

### 9.6 Current implementation coverage

`amy_flank_genomic_background.py` scans the selected chromosomes
from `GENOME_VCF`, builds per-site aggregate numerator/denominator
contributions, computes same-length 1-kb sliding-window backgrounds,
and applies chromosome-stratified block resampling for
**π, θW, exploratory Tajima D, Hudson FST**.

The genomewide ranks for within-LD, across-LD, normalized nSL
and normalized XP-nSL, the distance-matched LD enrichment,
the 93:93 empirical XP-nSL analysis, and bootstrap intervals
for selscan scores **remain outside the current implementation**.
Their definitions are recorded here to make the full manuscript
methods reproducible once those modules are supplied.

---

## 10. References

1. Bolognini D, et al. (2024). Recurrent evolution and selection shape structural diversity at the amylase locus. *Nature* 634:617–625. https://doi.org/10.1038/s41586-024-07911-1
2. Zheng-Bradley X, et al. (2017). Alignment of 1000 Genomes Project reads to reference assembly GRCh38. *GigaScience* 6:gix038. https://doi.org/10.1093/gigascience/gix038
3. Nei M, Li WH (1979). Mathematical model for studying genetic variation in terms of restriction endonucleases. *PNAS* 76:5269–5273. https://doi.org/10.1073/pnas.76.10.5269
4. Watterson GA (1975). On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology* 7:256–276. https://doi.org/10.1016/0040-5809(75)90020-9
5. Tajima F (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* 123:585–595. https://doi.org/10.1093/genetics/123.3.585
6. Ferrer-Admetlla A, et al. (2014). On detecting incomplete soft or hard selective sweeps using haplotype structure. *Molecular Biology and Evolution* 31:1275–1291. https://doi.org/10.1093/molbev/msu077
7. Szpiech ZA, Hernandez RD (2014). selscan: an efficient multithreaded program to perform EHH-based scans for positive selection. *Molecular Biology and Evolution* 31:2824–2827. https://doi.org/10.1093/molbev/msu211
8. Szpiech ZA (2024). selscan 2.0: scanning for sweeps in unphased data. *Bioinformatics* 40:btae006. https://doi.org/10.1093/bioinformatics/btae006
9. Hudson RR, Slatkin M, Maddison WP (1992). Estimation of levels of gene flow from DNA sequence data. *Genetics* 132:583–589. https://doi.org/10.1093/genetics/132.2.583
10. Bhatia G, et al. (2013). Estimating and interpreting FST: the impact of rare variants. *Genome Research* 23:1514–1521. https://doi.org/10.1101/gr.154831.113
11. Byrska-Bishop M, et al. (2022). High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell* 185:3426–3440.e19. https://doi.org/10.1016/j.cell.2022.08.004
12. Rahman A, Smith TQ, Szpiech ZA (2026). Selection scans and downstream analysis with selscan. *Human Population Genetics and Genomics* 6(2):0005. https://doi.org/10.47248/hpgg2606020005
13. Vollger MR, et al. (2023). Increased mutation and gene conversion within human segmental duplications. *Nature* 617:325–334. https://doi.org/10.1038/s41586-023-05895-y
14. Jeong H, et al. (2025). Structural polymorphism and diversity of human segmental duplications. *Nature Genetics* 57:390–401. https://doi.org/10.1038/s41588-024-02051-8

---

## Software

Bash; Python 3.9+ (standard-library focal/CI/genomic helpers);
bedtools, bgzip and tabix for reference-mask preparation;
bcftools 1.15+ for indexed VCF operations; selscan **v3.0.1**
for unphased nSL and phased XP-nSL. For final normalization,
`selscan norm` must be run on the appropriate autosomal
empirical score distribution, with an activated mamba/conda
environment supplying the exact executable version.
All helper scripts are in the repository root.
