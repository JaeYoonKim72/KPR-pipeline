# pangenome_variant_comparison.sh

Structural variants and small variants derived from the pangenome graph, and
their comparison between the Korean (KPR) and the HPRC haplotypes at the locus
and the allele level.

Part of the [KPR-pipeline](README.md) repository.

The variants of the KPR–HPRC graph, deconstructed against GRCh38, are split
into structural variants (SVs, more than 50 bp) and small variants (50 bp or
less), following the graph-based procedure of the Asian Pangenome Project
phase 1 (APGp1; Wu et al., 2026). SVs are decomposed with `vcfwave`, grouped
into loci with `vcfcreatemulti`, and similar alleles within each locus are
merged with PanSVMerger. Small variants are decomposed per sample with
`vcfwave` and `rtg vcfdecompose` and merged across samples. Every variant is
then counted as KPR-specific, shared or HPRC-specific.

---

## 1. Set up

### 1.1 Software

Install the tools listed under [Software](#software) and make sure they are on
`PATH`. Python needs `pysam`.

PanSVMerger is used from its repository, at the commit used in the study:

```bash
git clone https://github.com/Asian-Pan-Genome/PanSVMerger.git
git -C PanSVMerger checkout 9d9827f
```

`sv_merge_alleles.py` imports the clustering functions of
`PanSVMerger/scripts/SV.merge_variants.py` unchanged. Set `PANSVMERGER` if the
repository is not cloned next to the pipeline.

### 1.2 Input

The input is produced by
[`pangenome_construction_and_statistics.sh`](pangenome_construction_and_statistics.md):

```
results/
├── graph_variants/hprc_kpr/reference/hprc_kpr.GRCh38.vcf.gz   vg deconstruct -P GRCh38 -a
└── index/hprc_kpr/hprc_kpr.xg                                 GRCh38 sequence of the graph
```

The VCF must be the raw `vg deconstruct` output, with all nesting levels and
without `vcfbub`, because bubble selection is part of this pipeline.

### 1.3 Reference files

In `REF_DIR`, downloaded from UCSC:

| File | Source |
|---|---|
| `hg38.centromeres.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz |
| `hg38.gap.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz |

Centromeres and the telomere entries of the gap table are excluded.

### 1.4 Population table

`LIST_DIR/population.tsv` holds one sample per line, sample ID and population
separated by a tab. The sample IDs are the sample columns of the input VCF.
Samples of `TARGET_POP` (`KOR`) are compared with all other samples:

```
IBD-409   KOR
HG00438   CHN
HG02055   AFR
```

### 1.5 Environment

```bash
export REF_DIR=/path/to/reference     # required
export OUT_DIR=/path/to/results       # default: ./results
export LIST_DIR=/path/to/lists        # default: ./lists
export JOBS=4                         # chromosomes or samples in parallel
```

---

## 2. Run

Print the available steps:

```bash
bash pangenome_variant_comparison.sh
```

Run everything in order:

```bash
bash pangenome_variant_comparison.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash pangenome_variant_comparison.sh prepare_input
JOBS=23 bash pangenome_variant_comparison.sh sv_decompose
JOBS=4 MERGE_WORKERS=16 bash pangenome_variant_comparison.sh sv_merge
JOBS=24 bash pangenome_variant_comparison.sh small_decompose
JOBS=12 bash pangenome_variant_comparison.sh small_merge
JOBS=12 bash pangenome_variant_comparison.sh count
bash pangenome_variant_comparison.sh summary
```

Test on one chromosome first:

```bash
CHROMS=chr22 OUT_DIR=test bash pangenome_variant_comparison.sh all
```

SV steps run per chromosome, `small_decompose` per sample. `sv_decompose` runs
one single-threaded `vcfwave` per chromosome and needs up to 35 GB each;
`sv_merge` uses `MERGE_WORKERS` cores per chromosome. Each job writes its own
log to `$OUT_DIR/logs/variant_comparison/`, and every step writes to
`$OUT_DIR/logs/pipeline.log`. A step stops when any of its jobs fails.

---

## 3. Steps

| Step | Description |
|---|---|
| `prepare_input` | `vcfbub`, contig renaming, removal of the CHM13 haplotype, centromeres and telomeres, SV input, reference FASTA |
| `sv_decompose` | `vcfwave` decomposition, SV alleles of more than 50 bp, SV loci with `vcfcreatemulti` |
| `sv_merge` | PanSVMerger allele clustering within each SV locus |
| `small_decompose` | Per sample: `vcfwave`, `bcftools norm`, `rtg vcfdecompose`, small variants of 50 bp or less |
| `small_merge` | Per chromosome: sample merge, SNV and indel files |
| `count` | KPR-specific, shared and HPRC-specific variants per chromosome |
| `summary` | Summary tables at the locus and allele level |
| `all` | All of the above, in order |

Dependencies between steps: `sv_decompose` and `small_decompose` need
`prepare_input`; `sv_merge` needs `sv_decompose`; `small_merge` needs
`small_decompose`; `count` needs `sv_merge` and `small_merge`; `summary` needs
`count`. The SV branch and the small variant branch are independent and can
run at the same time.

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `REF_DIR` | required | reference file directory |
| `OUT_DIR` | `./results` | output directory, shared with the other pipelines |
| `LIST_DIR` | `./lists` | population table directory |
| `TMP_DIR` | `$OUT_DIR/tmp` | temporary files |
| `GRAPH` | `hprc_kpr` | graph whose variants are compared |
| `INPUT_VCF` | `$OUT_DIR/graph_variants/$GRAPH/reference/$GRAPH.GRCh38.vcf.gz` | deconstructed VCF |
| `GRAPH_XG` | `$OUT_DIR/index/$GRAPH/$GRAPH.xg` | graph used to extract the reference sequence |
| `POP_TABLE` | `$LIST_DIR/population.tsv` | sample and population |
| `TARGET_POP` | `KOR` | population compared with all other samples |
| `TARGET_LABEL`, `OTHER_LABEL` | `KPR`, `HPRC` | labels used in the output tables |
| `EXCLUDE_SAMPLES` | `CHM13v2` | reference haplotype stored as a sample column |
| `REF_SAMPLE` | `GRCh38` | reference path of the deconstructed VCF |
| `CHROMS` | `chr1`–`chr22`, `chrX` | chromosomes analysed |
| `VCFBUB_LEVEL` | `0` | `vcfbub -l` |
| `VCFBUB_MAX_REF` | `10000000` | `vcfbub -r` |
| `WAVE_MAXLEN` | `100000` | `vcfwave -L` |
| `SV_MIN_LEN` | `50` | SVs are longer, small variants not |
| `MERGE_IDENTITY` | `0.9` | vsearch identity used by PanSVMerger |
| `MERGE_WORKERS` | `16` | parallel vsearch processes per chromosome |
| `MULTI_MEM_GB` | `32` | memory cap of `vcfcreatemulti` |
| `RTG_MAXLEN` | `5000` | records with a longer allele are not passed to rtg |
| `RTG_MEM` | `8g` | rtg heap size |
| `SORT_MEM` | `8G` | `bcftools sort` buffer |
| `PANSVMERGER` | `./PanSVMerger/scripts/SV.merge_variants.py` | PanSVMerger script |
| `JOBS` | `4` | chromosomes or samples processed in parallel |

---

## 5. Output

```
results/variant_comparison/
├── input/            preprocessed VCF, SV input, sample list, reference FASTA
├── sv/CHR/           bi-allelic, SV, SV locus and merged SV VCFs per chromosome
├── small_variants/
│   ├── samples/      per-sample small variant VCFs
│   └── merged/CHR/   merged small variants, SNVs and indels per chromosome
├── counts/           per-chromosome counts
├── variant_comparison_summary.tsv
├── variant_comparison_by_class.tsv
└── variant_comparison_shared_loci.tsv
```

Main tables:

- `variant_comparison_summary.tsv` — KPR-specific, shared and HPRC-specific
  SVs and small variants at the locus and allele level, with the novel
  percentage
- `variant_comparison_by_class.tsv` — the same counts by variant class (SNV,
  insertion, deletion, complex). A complex locus holds alleles of more than one
  class, or multi-base substitutions
- `variant_comparison_shared_loci.tsv` — shared SV loci carrying a
  KPR-specific allele, and the number of such alleles
- `sv/CHR/CHR.sv.merged.vcf.gz` — SV loci after allele merging
- `small_variants/merged/CHR/CHR.small.vcf.gz` — merged small variants

---

## 6. What the steps do

The pipeline has two branches that share the same preprocessed input and meet
again at the counting step:

```
vg deconstruct -P GRCh38 -a  (graph_variants/hprc_kpr/reference/hprc_kpr.GRCh38.vcf.gz)
        │
        ▼
prepare_input     vcfbub -l 0 -r 10000000, CHM13 haplotype and centromeres/telomeres removed
        │
        ├──────────────────────────────────┐
        ▼                                  ▼
  SV (> 50 bp)                       Small variants (≤ 50 bp)
  all samples together, per chr      each sample separately
        │                                  │
sv_decompose   vcfwave                small_decompose   vcfwave
               norm -m-                                 norm -c s -f GRCh38 -m -any
               keep > 50 bp                             rtg vcfdecompose
               vcfcreatemulti (loci)                    keep ≤ 50 bp, merge haplotypes
        │                                  │
sv_merge       PanSVMerger             small_merge      bcftools merge (all samples)
               (vsearch id ≥ 0.9)                       SNV and indel files
        │                                  │
        └──────────────┬───────────────────┘
                       ▼
count          KPR-specific / shared / HPRC-specific, locus and allele
                       ▼
summary        tables, novel % = KPR-specific / (shared + HPRC-specific)
```

The numbered comments in the script (`1-1)`, `2-1)`, …) correspond to the
items below.

### Step 1. `prepare_input`

Input: `INPUT_VCF`, `GRAPH_XG`, the UCSC tables. Output: `variant_comparison/input/`.

- **1-1) Excluded regions.** Centromeres from the UCSC centromere table and the
  `telomere` entries of the UCSC gap table are written to
  `exclude.cen_telo.bed`.
- **1-2) Contig names.** The contigs of the deconstructed VCF are named after
  the reference path (`GRCh38#0#chr1`). A rename table maps them to `chr1`,
  so that the VCF agrees with the BED files and the reference FASTA.
- **1-3) Common preprocessing.** This is the only filtering shared by both
  branches.
  - `vcfbub -l 0 -r 10000000` keeps the top-level bubbles. A top-level bubble
    whose reference allele is longer than 10 Mb is replaced by its child
    bubbles, so that a few very large bubbles do not hide the variants nested
    inside them. This is the APGp1 setting.
  - The CHM13 haplotype, stored as a sample column, is removed. It is a
    reference, not a sample of either group.
  - Only `chr1`–`chr22` and `chrX` are kept, and sites overlapping a
    centromere or telomere are removed.
  - AC and AN are recomputed and sites no longer carried by any sample are
    removed (`bcftools view -c 1`).

  Output: `all.vcf.gz` and the sample list `samples.txt`.
- **1-4) SV input.** Sites whose reference or an alternative allele is longer
  than 50 bp are copied to `sv_input.vcf.gz`. A site whose alleles are all
  50 bp or shorter cannot contain an allele differing from the reference by
  more than 50 bp, so it cannot yield an SV after decomposition. The filter only
  shortens the `vcfwave` run of the SV branch and does not change its result.
  It is not part of APGp1.
- **1-5) Reference FASTA.** The GRCh38 sequence is extracted from the graph
  with `vg paths -S GRCh38 -F`, renamed to `chr1` … and indexed. It is used by
  `bcftools norm -f` in the small variant branch.

### Step 2. `sv_decompose` (per chromosome)

Input: `input/sv_input.vcf.gz`. Output: `sv/CHR/`.

- **2-1) Decomposition.** A record deconstructed from a bubble holds the whole
  sequence of each path through the bubble, so one long allele can hide
  several variants. `vcfwave -L 100000` aligns every alternative allele to the
  reference allele and writes each primitive difference (SNV, indel, SV,
  inversion) as its own record. For example, an alternative allele that differs
  from the reference by one substitution and one 3-bp deletion becomes two
  records, an SNV and a deletion. Alleles longer than 100 kb are not realigned
  and are kept as they are. The records are then split into bi-allelic records
  (`bcftools norm -m-`), spanning-deletion alleles (`*`) are removed, the file
  is sorted, and AC/AN are recomputed. Output: `CHR.bi.vcf.gz`.
- **2-2) SV alleles.** Bi-allelic records with |len(ALT) − len(REF)| > 50 bp
  are kept. Output: `CHR.sv.bi.vcf.gz`.
- **2-3) SV loci.** After decomposition, the same SV in different haplotypes
  is often written as several records with slightly different breakpoints or
  lengths, for example three insertions of 300, 305 and 298 bp written at the
  same position. `vcfcreatemulti` merges records with overlapping reference
  spans into one multi-allelic record. Each such record is one **SV locus**. Output:
  `CHR.sv.multi.vcf.gz`.

### Step 3. `sv_merge` (per chromosome)

Input: `sv/CHR/CHR.sv.multi.vcf.gz`. Output: `sv/CHR/CHR.sv.merged.vcf.gz`.

- **3-1) Allele clustering.** Within a locus, alleles that differ only by a
  few bases are still separate alleles, which inflates the number of distinct
  SV alleles. `sv_merge_alleles.py` applies the PanSVMerger clustering to every
  multi-allelic locus:
  - the reference and alternative sequences are clustered with
    `vsearch --cluster_fast --id 0.9 --strand both`;
  - a cluster is split when its sequences differ in length by more than 50 bp;
  - the cluster containing the reference allele becomes the reference, and
    every other cluster becomes one alternative allele, represented by its
    longest sequence;
  - genotypes are rewritten to the new allele numbers, and AC and AN are
    recomputed;
  - a locus whose alternative alleles all cluster with the reference is
    removed.

  In the example above, the 300, 305 and 298 bp insertions form one cluster
  and the locus is left with one alternative allele. Each cluster is one
  **SV allele**.

  The clustering functions are imported from PanSVMerger
  (`SV.merge_variants.py`) unchanged, with its default limits (a locus with
  more than 100 alleles of 50 kb or longer is not clustered). Compared with
  running PanSVMerger directly, the wrapper keeps the diploid genotypes as they
  are and clusters the loci of a chromosome in parallel
  (`MERGE_WORKERS`); the result is the same as serial processing.

### Step 4. `small_decompose` (per sample)

Input: `input/all.vcf.gz`, `input/GRCh38.fa`. Output:
`small_variants/samples/SAMPLE/SAMPLE.small.vcf.gz`.

Small variants are processed one sample at a time, following the graph-based
processing of loss-of-function variants in APGp1.

- **4-1) Extraction and normalisation.**
  `bcftools view -a -I -s SAMPLE -c 1:nref` keeps the sites carried by the
  sample and trims the alleles it does not carry. `vcfwave -L 100000`
  decomposes the alleles as in step 2-1, and `bcftools norm -c s -f GRCh38
  -m -any` left-aligns the variants against GRCh38 and splits them into
  bi-allelic records.
- **4-2) Breaking MNPs and complex indels.**
  `rtg vcfdecompose --break-mnps --break-indels` writes multi-base
  substitutions as separate SNVs and complex indels as an SNV plus a simple
  indel, so that a small variant has the same representation in every sample.
  Before rtg, the graph path annotation `INFO/AT` (up to megabytes per record)
  is removed, and inversions and records with an allele longer than 5 kb are
  set aside. These records are SVs; rtg splits an inversion into thousands of
  artificial SNVs, and its memory grows steeply with the allele length. This
  filter is not part of APGp1.
- **4-3) Small variants.** Records with |len(ALT) − len(REF)| ≤ 50 bp are
  kept. After decomposition, an allele carried on both haplotypes of a sample
  can be written as two records (`1|0` and `0|1`). `dedup_sample.py` merges
  them into one record (`1|1`), so that a homozygous variant is not counted as
  two heterozygous ones. The number of merged records is written to
  `SAMPLE.dedup.txt`.

### Step 5. `small_merge` (per chromosome)

Input: the per-sample VCFs. Output: `small_variants/merged/CHR/`.

- **5-1) Sample merge.** `bcftools merge -m none --missing-to-ref` merges the
  samples without combining different alleles into one record. Each sample file
  holds only the variants that sample carries, so a sample without a record at
  a site carries the reference allele there. AC, AN and AF are recomputed.
  Output: `CHR.small.vcf.gz`.
- **5-2) SNVs and indels.** `bcftools view -v snps` and `-v indels` write
  `CHR.snv.vcf.gz` and `CHR.indel.vcf.gz`.

### Step 6. `count` (per chromosome)

Input: the merged SV and small variant files and `POP_TABLE`. Output:
`counts/CHR.tsv`.

**6-1) Group comparison.** `count_groups.py` reads the genotypes of every
record. A group carries an allele when at least one haplotype of one of its
samples carries it; missing genotypes are ignored. Every locus and every allele
is then classified as carried by KPR only, by both groups (shared) or by HPRC
only. Two units are counted for each variant type:

| | Locus | Allele |
|---|---|---|
| SV | one record of `CHR.sv.merged.vcf.gz` (`--mode locus`) | one PanSVMerger allele cluster of that record |
| Small variant | records of `CHR.small.vcf.gz` whose reference spans overlap (`--mode overlap`) | one distinct CHROM/POS/REF/ALT of `CHR.snv.vcf.gz` and `CHR.indel.vcf.gz` (`--mode exact`) |

A locus is shared when both groups carry any alternative allele at it, even if
they carry different alleles. For example, an SV locus with one allele carried
only by KPR samples and another carried only by HPRC samples counts as one
shared locus, one KPR-specific allele and one HPRC-specific allele. Such shared
loci with a KPR-specific allele are also counted separately, as new alleles at
known SV sites.

Every count is also split by class: SNV, insertion or deletion by the lengths
of REF and ALT, and complex for a locus whose alleles are of more than one
class or are multi-base substitutions.

### Step 7. `summary`

Input: `counts/`. Output: the three tables in `variant_comparison/`.

**7-1)** `summarize_variant_comparison.py` sums the chromosome counts and
computes the novel percentage, the number of KPR-specific variants relative
to all variants found in HPRC:

```
novel % = KPR-specific / (shared + HPRC-specific) × 100
```

### Running the helper scripts on their own

```bash
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' chr22.sv.merged.vcf.gz \
  | python3 count_groups.py --mode locus --varset sv --pop lists/population.tsv --target KOR
python3 summarize_variant_comparison.py --counts results/variant_comparison/counts \
  --out results/variant_comparison
```

### Counts obtained in the study

KPR (n = 14) versus HPRC (n = 47), GRCh38, `variant_comparison_summary.tsv`:

| Variant | Unit | KPR-specific | Shared | HPRC-specific | Novel (%) |
|---|---|---:|---:|---:|---:|
| SV (> 50 bp) | locus | 7,778 | 45,463 | 50,660 | 8.09 |
| SV (> 50 bp) | allele | 18,118 | 59,350 | 110,424 | 10.67 |
| Small variant (≤ 50 bp) | locus | 2,059,844 | 8,956,135 | 14,635,890 | 8.73 |
| Small variant (≤ 50 bp) | allele | 2,359,766 | 10,222,061 | 17,180,351 | 8.61 |

In addition, 10,211 KPR-specific SV alleles were found in 4,783 SV loci shared
with HPRC (`variant_comparison_shared_loci.tsv`).

---

## Software

vcfbub 0.1.0 · vcflib 1.0.15 (vcfwave, vcfcreatemulti) · PanSVMerger (commit
9d9827f) · vsearch 2.32.0 · RTG Tools 3.13 · bcftools 1.24 · samtools ·
vg 1.63.1 · Python 3 with pysam 0.22.1

All executables must be available on `PATH`.
