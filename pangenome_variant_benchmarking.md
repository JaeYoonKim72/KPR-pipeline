# pangenome_variant_benchmarking.sh

Benchmarking of the short-read variant calls obtained with the pangenome graphs
and with GRCh38, against long-read truth sets of seven Korean individuals.

Part of the [KPR-pipeline](README.md) repository.

For every individual, HiFi reads and a haplotype-resolved assembly are used to
build truth sets for small variants (SNPs, small deletions and insertions) and
for structural variants (SVs, deletions and insertions). The calls of the KPR,
HPRC+KPR and HPRC graphs, and the GRCh38 calls of GATK (small variants), Manta
and Delly (SVs), are then compared with the truth sets: small variants with
hap.py (vcfeval engine) and SVs with Truvari.

---

## 1. Set up

### 1.1 Software

Install the tools listed under [Software](#software) and make sure they are on
`PATH`. Python needs `pysam`.

### 1.2 Reference files

In `REF_DIR`:

| File | Description |
|---|---|
| `GRCh38.fa` | GRCh38 with `chr` contig names, indexed with `samtools faidx` |
| `GRCh38.sdf` | RTG SDF of the same FASTA, `rtg format -o GRCh38.sdf GRCh38.fa` |

The autosome BED (`chr1`–`chr22`) used by every step is written from
`GRCh38.fa.fai`.

### 1.3 Truth-set input

In `TRUTH_DIR`, one set of files per sample (`SAMPLE`), from HiFi reads aligned
to GRCh38 and from the haplotype-resolved assembly:

```
truth_input/
├── hifi_deepvariant/
│   ├── SAMPLE.hifi_dv.biSNP.addfilter.VAF0208.vcf.gz   DeepVariant, bi-allelic SNPs
│   └── SAMPLE.hifi_dv.Small.vcf.gz                     DeepVariant, small variants
└── long_read_sv/
    ├── svim_asm_SAMPLE.PASS.SVLEN51.chr.vcf.gz         SVIM-asm, assembly-based SVs
    ├── pbsv_SAMPLE.PASS.SVLEN51.chr.vcf.gz             read-based SVs
    ├── sniffles_SAMPLE.PASS.SVLEN51.chr.vcf.gz
    └── svim_SAMPLE.PASS.SVLEN51.chr.vcf.gz
```

The SV files hold PASS calls of at least 51 bp on `chr`-named contigs.

### 1.4 Query input

In `QUERY_DIR`, one subdirectory per call set, with files named
`CALLSET.SAMPLE.KIND.vcf.gz`:

```
query_input/
├── kpr/        kpr.SAMPLE.snp.vcf.gz  kpr.SAMPLE.indel.vcf.gz  kpr.SAMPLE.sv.DEL.vcf.gz  kpr.SAMPLE.sv.INS.vcf.gz
├── hprc_kpr/   same four files
├── hprc/       same four files
├── gatk/       gatk.SAMPLE.snp.vcf.gz  gatk.SAMPLE.indel.vcf.gz
├── manta/      manta.SAMPLE.sv.DEL.vcf.gz  manta.SAMPLE.sv.INS.vcf.gz  manta.SAMPLE.sv.DUP.vcf.gz
└── delly/      delly.SAMPLE.sv.DEL.vcf.gz
```

The graph call sets are produced by
[`pangenome_short_read_mapping.sh`](pangenome_short_read_mapping.md), in
GRCh38 coordinates. For a Korean sample, the KPR and HPRC+KPR call sets come
from graphs built without that sample.

### 1.5 Sample list

`LIST_DIR/benchmark_samples.list` holds one sample ID per line, as used in the
input file names. The study used seven Korean individuals:

```
KOR-001
KOR-002
KOR-006
KOR-007
KOR-008
KOR-009
KOR-010
```

### 1.6 Environment

```bash
export REF_DIR=/path/to/reference     # required
export TRUTH_DIR=/path/to/truth_input # required
export QUERY_DIR=/path/to/query_input # required
export OUT_DIR=/path/to/results       # default: ./results
export LIST_DIR=/path/to/lists        # default: ./lists
export JOBS=8                         # samples or call sets in parallel
export THREADS=4                      # threads per hap.py run
```

---

## 2. Run

Print the available steps:

```bash
bash pangenome_variant_benchmarking.sh
```

Run everything in order:

```bash
bash pangenome_variant_benchmarking.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash pangenome_variant_benchmarking.sh small_truth
bash pangenome_variant_benchmarking.sh sv_truth_filter
bash pangenome_variant_benchmarking.sh sv_truth_support
bash pangenome_variant_benchmarking.sh sv_truth_readsupport
bash pangenome_variant_benchmarking.sh small_query
JOBS=6 THREADS=8 bash pangenome_variant_benchmarking.sh small_bench
bash pangenome_variant_benchmarking.sh sv_query
bash pangenome_variant_benchmarking.sh sv_bench
bash pangenome_variant_benchmarking.sh summary
```

Evaluate a single graph by restricting `GRAPHS`:

```bash
GRAPHS=kpr bash pangenome_variant_benchmarking.sh sv_bench
```

Every job writes its own log to `$OUT_DIR/logs/benchmarking/`, and every step
writes to `$OUT_DIR/logs/pipeline.log`. A step stops when any of its jobs
fails.

---

## 3. Steps

| Step | Description |
|---|---|
| `small_truth` | HiFi DeepVariant small variant truth sets (SNP, DEL, INS) |
| `sv_truth_filter` | Filtering of the SVIM-asm, pbsv, Sniffles and SVIM calls |
| `sv_truth_support` | SVIM-asm SVs supported by at least 1, 2 or 3 read-based callers |
| `sv_truth_readsupport` | Read-support filter, deletion and insertion truth sets |
| `small_query` | Small variant queries (SNP, DEL, INS) |
| `small_bench` | hap.py benchmark of the small variants |
| `sv_query` | SV queries (DEL, INS) |
| `sv_bench` | Truvari benchmark of the SVs |
| `summary` | Summary tables and paired Wilcoxon tests |
| `all` | All of the above, in order |

Dependencies between steps: `sv_truth_support` needs `sv_truth_filter`;
`sv_truth_readsupport` needs `sv_truth_support`; `small_bench` needs
`small_truth` and `small_query`; `sv_bench` needs `sv_truth_readsupport` and
`sv_query`; `summary` needs `small_bench` and `sv_bench`. The small variant
branch and the SV branch are independent and can run at the same time.

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `REF_DIR` | required | reference file directory |
| `TRUTH_DIR` | required | truth-set input directory |
| `QUERY_DIR` | required | query input directory |
| `OUT_DIR` | `./results` | output directory, shared with the other pipelines |
| `LIST_DIR` | `./lists` | sample list directory |
| `TMP_DIR` | `$OUT_DIR/tmp` | temporary files |
| `GRCH38` | `$REF_DIR/GRCh38.fa` | reference FASTA |
| `GRCH38_SDF` | `$REF_DIR/GRCh38.sdf` | RTG SDF used by vcfeval |
| `SAMPLE_LIST` | `$LIST_DIR/benchmark_samples.list` | samples evaluated |
| `HIFI_DV_DIR` | `$TRUTH_DIR/hifi_deepvariant` | HiFi DeepVariant calls |
| `LR_SV_DIR` | `$TRUTH_DIR/long_read_sv` | long-read SV calls |
| `GRAPHS` | `kpr hprc_kpr hprc` | graph call sets |
| `LINEAR_SMALL` | `gatk` | GRCh38 small variant call set |
| `SV_SUPPORT_LEVELS` | `ge1 ge2 ge3` | SV truth sets evaluated |
| `TRUTH_MIN_GQ`, `TRUTH_MIN_DP` | `20`, `7` | genotype quality and depth of the small variant truth |
| `HET_VAF_MIN`, `HET_VAF_MAX` | `0.30`, `0.70` | VAF range of heterozygous truth genotypes |
| `HOM_VAF_MIN` | `0.90` | minimum VAF of homozygous alternative truth genotypes |
| `SUPPORT_REFDIST`, `SUPPORT_PCTSIZE`, `SUPPORT_PCTSEQ` | `500`, `0.7`, `0.7` | SVIM-asm versus read-based caller matching |
| `READSUP_MIN_ALT`, `READSUP_MIN_DP` | `5`, `10` | read support of an SV truth call |
| `TV_REFDIST`, `TV_PCTSIZE`, `TV_PCTSEQ`, `TV_PCTOVL` | `500`, `0.8`, `0.8`, `0.8` | Truvari matching of the SV benchmark |
| `TV_SIZEMIN`, `TV_SIZEFILT`, `TV_SIZEMAX` | `51`, `30`, `50000` | Truvari size range |
| `SV_MIN_LEN` | `51` | minimum length of an SV truth call |
| `JOBS` | `8` | samples or call sets processed in parallel |
| `THREADS` | `4` | threads per hap.py run |

---

## 5. What the steps do

The pipeline has two branches, one for small variants and one for SVs. Each
builds truth sets and queries separately and compares them at the benchmark
step:

```
HiFi DeepVariant                     SVIM-asm + pbsv / Sniffles / SVIM
        │                                          │
        ▼                                          ▼
small_truth   GQ, DP, VAF filters    sv_truth_filter       PASS INS/DEL, read support
        │                                          ▼
        │                            sv_truth_support      Truvari, ge1 / ge2 / ge3
        │                                          ▼
        │                            sv_truth_readsupport  ALT ≥ 5 and DP ≥ 10
        │                                          │
small_query   graph and GATK calls   sv_query              graph, Manta and Delly calls
        │                                          │
        ▼                                          ▼
small_bench   hap.py (vcfeval)       sv_bench              Truvari bench
        │                                          │
        └──────────────────┬───────────────────────┘
                           ▼
summary        precision, recall, F1, paired Wilcoxon tests
```

The numbered comments in the script (`1-1)`, `2-1)`, …) correspond to the
items below.

### Step 1. `small_truth` (per sample)

Input: `HIFI_DV_DIR`. Output: `truth/small/SAMPLE/`.

High-confidence genotypes of the HiFi DeepVariant calls on `chr1`–`chr22` are
kept:

| Filter | SNP | Small deletion and insertion |
|---|---|---|
| Records | `bcftools norm -f GRCh38 -m -any`, bi-allelic, PASS | `-v indels` |
| Genotype quality and depth | GQ ≥ 20, DP ≥ 7 | GQ ≥ 20, DP ≥ 7 |
| VAF | heterozygous 0.30–0.70, homozygous alternative ≥ 0.90 | the same |
| VAF source | `bcftools +fill-tags -t FORMAT/VAF`, from AD | DeepVariant `FORMAT/VAF` |

- **1-1) SNPs.** Output: `SAMPLE.SNP.truth.vcf.gz`.
- **1-2) Small deletions and insertions.** Records with `strlen(REF) >
  strlen(ALT)` are deletions and the others insertions. Output:
  `SAMPLE.DEL.truth.vcf.gz`, `SAMPLE.INS.truth.vcf.gz`.

### Step 2. `sv_truth_filter` (per sample and caller)

Input: `LR_SV_DIR`. Output: `truth/sv/SAMPLE/SAMPLE.CALLER.alt.vcf.gz`.

- **2-1) Call filtering.** `sv_truth_helper.py filter` keeps PASS insertions
  and deletions of at least 50 bp on `chr1`–`chr22`, `chrX` and `chrY` with a
  complete non-reference genotype. The read-based callers also need read
  support:

  | Caller | Additional filter |
  |---|---|
  | SVIM-asm | none |
  | pbsv | DP ≥ 7, ALT reads (AD[1]) ≥ 3 |
  | SVIM | DP ≥ 7, ALT reads (AD[1]) ≥ 3, QUAL ≥ 10 |
  | Sniffles | GQ ≥ 20, DR + DV ≥ 7, DV ≥ 3 |

  The SVIM-asm records are given unique IDs (`SAMPLE_ASM_00000001`, …), used
  in step 3.

### Step 3. `sv_truth_support` (per sample)

Input and output: `truth/sv/SAMPLE/`.

- **3-1) Matching.** The SVIM-asm calls are compared with each read-based call
  set:

  ```bash
  truvari bench -b SAMPLE.svim_asm.alt.vcf.gz -c SAMPLE.CALLER.alt.vcf.gz -f GRCh38.fa \
    --refdist 500 --pctsize 0.7 --pctseq 0.7 --sizemin 50 -o bench_CALLER
  ```

- **3-2) Support levels.** `sv_truth_helper.py collect` counts the read-based
  callers whose `tp-base` holds each SVIM-asm SV (`SAMPLE.SV.support.tsv`).
  `levels` writes three truth sets, with the assembly genotypes:

  | Level | SVIM-asm SVs supported by |
  |---|---|
  | `ge1` | at least 1 read-based caller |
  | `ge2` | at least 2 read-based callers |
  | `ge3` | all 3 read-based callers |

### Step 4. `sv_truth_readsupport` (per sample and level)

Input: `truth/sv/SAMPLE/SAMPLE.SV.truth.LEVEL.vcf.gz` and the raw read-based
calls. Output: `truth/sv/SAMPLE/LEVEL/`.

- **4-1) Read support.** `sv_truth_readsupport.py` matches every truth SV
  again to the Sniffles, pbsv and SVIM calls of the sample. A match needs the
  same SV type, start positions within 500 bp, size similarity ≥ 0.5 and, for
  insertions, inserted-sequence similarity ≥ 0.5. For every caller the best
  match is used:

  | Caller | Depth | ALT reads |
  |---|---|---|
  | Sniffles | DR + DV | DV |
  | pbsv | FORMAT/DP | AD[1] |
  | SVIM | FORMAT/DP | AD[1] |

  A truth SV is kept only when one caller has a match with ALT reads ≥ 5 and
  depth ≥ 10. Both values must come from the same caller. SVs without such a
  match are removed, including SVs whose matched calls do not report depth or
  ALT reads. Output: `truth.readsupport.vcf.gz`, the per-SV values in
  `evidence.tsv` and the counts in `filter_counts.tsv`.
- **4-2) Deletions and insertions.** The INFO/SVLEN header is set to
  `Number=A`, records are split with `bcftools norm -m -any` and sorted, and
  `prepare_sv_autosome.py` keeps the deletions or the insertions of at least
  51 bp on `chr1`–`chr22`. Output: `truth.DEL.vcf.gz`, `truth.INS.vcf.gz`.

### Step 5. `small_query` (per sample and call set)

Input: `QUERY_DIR`. Output: `query/small/SAMPLE/`.

Every call set is restricted to `chr1`–`chr22` and to PASS or missing FILTER,
normalised with `bcftools norm -f GRCh38 -m -any`, and limited to bi-allelic
records with a heterozygous or homozygous alternative genotype without a
missing allele (`(GT="het" || GT="AA") && GT!~"\."`).

- **5-1) SNPs.** Output: `CALLSET.SNP.vcf.gz`.
- **5-2) Small deletions and insertions.** Output: `CALLSET.DEL.vcf.gz`,
  `CALLSET.INS.vcf.gz`.

### Step 6. `small_bench` (per sample, call set and type)

Output: `bench/small/SAMPLE/CALLSET.TYPE.*`.

- **6-1) hap.py.**

  ```bash
  hap.py SAMPLE.TYPE.truth.vcf.gz CALLSET.TYPE.vcf.gz -r GRCh38.fa -R autosomes.bed \
    -o CALLSET.TYPE --engine vcfeval --engine-vcfeval-template GRCh38.sdf --no-roc
  ```

  SNPs are restricted to the autosomes with `-R`; small deletions and
  insertions use the same BED as confident regions (`-f`). With `-f`, the
  hap.py quantification fails at the size of the SNP truth sets. Metrics are
  read from the `Filter = ALL` row of `CALLSET.TYPE.summary.csv`.

### Step 7. `sv_query` (per sample)

Input: `QUERY_DIR`. Output: `query/sv/SAMPLE/`.

Every call set keeps PASS records, is split with `bcftools norm -m -any`, and
keeps non-reference genotypes without a missing allele (`GT="alt" &&
GT!~"\."`).

- **7-1) Graphs.** The deletion and insertion files are used as they are.
- **7-2) Manta.** The insertion query holds the insertion and the duplication
  calls together, so that tandem duplications can match truth insertions.
- **7-3) Delly.** Deletions only.

### Step 8. `sv_bench` (per sample, call set, type and level)

Output: `bench/sv/SAMPLE/CALLSET/LEVEL/TYPE/`.

- **8-1) Truvari.** Every query is compared with the `ge1`, `ge2` and `ge3`
  truth sets:

  ```bash
  truvari bench -b truth.TYPE.vcf.gz -c CALLSET.TYPE.vcf.gz -o TYPE --includebed autosomes.bed \
    --refdist 500 --pctsize 0.8 --pctseq 0.8 --pctovl 0.8 \
    --sizemin 51 --sizefilt 30 --sizemax 50000 --passonly --pick single
  ```

  | Option | Meaning |
  |---|---|
  | `--refdist 500` | breakpoints within 500 bp |
  | `--pctsize 0.8` | size similarity ≥ 80% |
  | `--pctseq 0.8` | sequence similarity ≥ 80% |
  | `--pctovl 0.8` | reciprocal overlap ≥ 80% |
  | `--sizemin 51 --sizefilt 30 --sizemax 50000` | truth SVs of 51 bp to 50 kb; query calls of at least 30 bp can match |
  | `--pick single` | one-to-one matching |
  | `--typeignore` | insertions only, so that Manta duplications can match truth insertions |

### Step 9. `summary`

Input: `bench/`. Output: `summary/`.

**9-1)** `summarize_benchmark.py` collects TP, FP, FN, precision, recall and
F1 of every run, computes their mean and SD over the samples, and compares
the graphs with paired two-sided exact Wilcoxon signed-rank tests across the
samples. With seven samples the smallest possible p-value is 0.016, reached
when all differences go in the same direction.

### Running the helper scripts on their own

```bash
python3 sv_truth_readsupport.py --truth SAMPLE.SV.truth.ge1.vcf.gz \
  --sniffles sniffles_SAMPLE.PASS.SVLEN51.chr.vcf.gz \
  --pbsv pbsv_SAMPLE.PASS.SVLEN51.chr.vcf.gz \
  --svim svim_SAMPLE.PASS.SVLEN51.chr.vcf.gz \
  --min-alt 5 --min-dp 10 --outdir ge1
python3 summarize_benchmark.py --outdir results/benchmarking \
  --samples KOR-001 KOR-002 --graphs kpr hprc_kpr hprc
```

---

## 6. Output

```
results/benchmarking/
├── autosomes.bed
├── truth/
│   ├── small/SAMPLE/       SNP, DEL and INS truth sets
│   └── sv/SAMPLE/          filtered calls, Truvari support runs, support table,
│       └── LEVEL/          ge1, ge2 and ge3 truth sets, read-support evidence
├── query/
│   ├── small/SAMPLE/       SNP, DEL and INS queries
│   └── sv/SAMPLE/          DEL and INS queries
├── bench/
│   ├── small/SAMPLE/       hap.py output
│   └── sv/SAMPLE/CALLSET/LEVEL/TYPE/   Truvari output
└── summary/
```

Main tables:

- `summary/kor_small_detail.tsv` — hap.py truth and query totals, TP, FN, FP,
  precision, recall and F1 per sample, call set and type
- `summary/kor_sv_detail.tsv` — the same for Truvari, per sample, call set,
  support level and type
- `summary/kor_small_mean.tsv`, `summary/kor_sv_mean.tsv` — means and SDs over
  the samples
- `summary/kor_wilcoxon.tsv` — paired Wilcoxon tests between the graphs
- `truth/sv/SAMPLE/LEVEL/filter_counts.tsv` — SVs kept by the read-support
  filter

These tables correspond to Tables S28–S30 of the manuscript.

---

## Software

hap.py 0.3.15 · RTG Tools 3.13 · Truvari 5.4.0 · bcftools 1.17 · samtools ·
Python 3 with pysam 0.23.3

All executables must be available on `PATH`.
