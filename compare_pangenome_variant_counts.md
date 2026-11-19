# compare_pangenome_variant_counts.sh

Compares structural variants (SVs) and small variants from the KPR–HPRC pangenome graph at both the locus and allele levels, classifying variants as KPR-specific, shared, or HPRC-specific.

Part of the [KPR-pipeline](README.md) repository.

The input is the Minigraph-Cactus [2] graph built from the KPR and HPRC [3] assemblies, with variants deconstructed against GRCh38 using vg [4]. Two counting approaches are used: locus-based counting treats each variable site as one locus regardless of the number of alternative alleles, whereas allele-based counting counts distinct alternative alleles separately. We referred to the Asian Pangenome Project phase 1 (APGp1) [1], which describes both approaches in detail. At the locus level, we identified 7.8K KPR-specific SVs and 2.06M KPR-specific small variants; at the allele level, we identified 18.1K KPR-specific SV alleles and 2.36M KPR-specific small variant alleles.

---

## 1. Set up

### 1.1 Software

Install the programs under [Software](#software) and make sure they are on
`PATH`. The versions used in the study can be installed into one
mamba/conda environment:

```bash
mamba create -n kpr_variant_counts -c conda-forge -c bioconda \
  python=3.8 bcftools=1.24 htslib=1.24 samtools=1.24 vcflib=1.0.15 \
  vcfbub=0.1.2 vsearch=2.32.0 pysam=0.22.1 rtg-tools=3.13 vg=1.63.1
conda activate kpr_variant_counts
```

PanSVMerger [1] is not packaged and is used from its repository, at the
commit used in the study:

```bash
git clone https://github.com/Asian-Pan-Genome/PanSVMerger.git
git -C PanSVMerger checkout 9d9827f
```

`sv_merge_alleles.py` imports the clustering functions of
`PanSVMerger/scripts/SV.merge_variants.py` unchanged. Set `PANSVMERGER` if the
repository is not cloned next to the pipeline.

RTG Tools asks about usage logging the first time it runs, which can halt a
background job. Answer it once interactively, or set `RTG_TALKBACK=false` and
`RTG_USAGE=false` in the `rtg.cfg` of the installation, before running the
pipeline.

### 1.2 Input

The input is produced by
[`pangenome_construction_and_statistics.sh`](pangenome_construction_and_statistics.md):

```
results/
├── graph_variants/hprc_kpr/reference/hprc_kpr.GRCh38.vcf.gz   vg deconstruct -P GRCh38 -a
└── index/hprc_kpr/hprc_kpr.xg                                 GRCh38 sequence of the graph
```

The VCF is the raw `vg deconstruct` output, without `vcfbub`; bubble selection
is done in `prepare_input`.

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
bash compare_pangenome_variant_counts.sh
```

Run everything in order:

```bash
bash compare_pangenome_variant_counts.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash compare_pangenome_variant_counts.sh prepare_input
JOBS=23 bash compare_pangenome_variant_counts.sh sv_decompose
JOBS=4 MERGE_WORKERS=16 bash compare_pangenome_variant_counts.sh sv_merge
JOBS=24 bash compare_pangenome_variant_counts.sh small_decompose
JOBS=12 bash compare_pangenome_variant_counts.sh small_merge
JOBS=12 bash compare_pangenome_variant_counts.sh count
bash compare_pangenome_variant_counts.sh summary
```

Test on one chromosome first. With a separate `OUT_DIR`, point `INPUT_VCF`
and `GRAPH_XG` to the existing input, since their defaults are under `OUT_DIR`:

```bash
CHROMS=chr22 OUT_DIR=test \
INPUT_VCF=results/graph_variants/hprc_kpr/reference/hprc_kpr.GRCh38.vcf.gz \
GRAPH_XG=results/index/hprc_kpr/hprc_kpr.xg \
bash compare_pangenome_variant_counts.sh all
```

SV steps run per chromosome, `small_decompose` per sample. `sv_decompose` runs
one single-threaded `vcfwave` per chromosome and needs up to 35 GB each;
`sv_merge` uses `MERGE_WORKERS` cores per chromosome. Each job writes its own
log to `$OUT_DIR/logs/variant_comparison/`, and every step writes to
`$OUT_DIR/logs/pipeline.log`. If any job of a step fails, the step exits with
an error once its remaining jobs have finished, and the following steps are not
run.

---

## 3. Steps

| Step | Description |
|---|---|
| `prepare_input` | Select top-level bubbles with `vcfbub`, standardize contig names, remove the CHM13 sample and variants located in centromeres and telomeres, and prepare the SV input and the reference sequence |
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
| `RTG_MAXLEN` | `5000` | records with a longer allele, and inversions, are not passed to rtg |
| `RTG_MEM` | `8g` | rtg heap size |
| `SORT_MEM` | `8G` | `bcftools sort` buffer |
| `PANSVMERGER` | `./PanSVMerger/scripts/SV.merge_variants.py` | PanSVMerger script |
| `JOBS` | `4` | chromosomes or samples processed in parallel |

---

## 5. What the steps do

The pipeline preprocesses the `vg deconstruct` VCF of the graph, processes SVs and small variants separately, and combines their results for locus- and allele-level comparisons between KPR and HPRC.

```text
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
               bcftools norm -m-                        bcftools norm -c s -f GRCh38 -m -any
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

### Prepare input

Bubbles of the `vg deconstruct` VCF are selected with `vcfbub -l 0 -r 10000000`, as in APGp1 [1]: top-level bubbles are kept, and a bubble whose reference allele is longer than 10 Mb is replaced by its child bubbles. Contig names are standardized, and the CHM13 reference haplotype is removed. Variants whose start position (POS) lies in a centromere or telomere are excluded. Sites with no carrier are removed, and the SV input and the GRCh38 reference sequence are prepared.

### SV decomposition (per chromosome)

Graph-derived variants are decomposed with `vcfwave` [5] and split into bi-allelic records with `bcftools norm -m-` [6], and variants with a REF–ALT length difference greater than 50 bp are retained as SVs. Records with overlapping reference spans are consolidated with `vcfcreatemulti` [5] to define SV loci.

### SV allele merging (per chromosome)

Within each SV locus, PanSVMerger [1] and VSEARCH [7] group similar alternative alleles based on sequence identity and length similarity. The resulting allele clusters are used for allele-level comparison, with genotypes and allele counts updated accordingly.

### Small variant decomposition (per sample)

Following the graph-based small variant processing of APGp1 [1], each sample is processed separately using `vcfwave` [5], `bcftools norm` [6], and `rtg vcfdecompose` [8] to obtain consistently represented SNVs and indels. Inversions and records with an allele longer than 5 kb (`RTG_MAXLEN`) are not passed to `rtg vcfdecompose`, so small variants within such records are not counted. Variants with a REF–ALT length difference of 50 bp or less are retained, and duplicate records from the two haplotypes of the same sample are combined.

### Small variant merging (per chromosome)

The per-sample small variant VCFs are combined using `bcftools merge`, followed by recalculation of allele statistics. SNVs and indels are also written to separate files for downstream counting.

### Variant counting criteria

At the locus level, variants with overlapping reference spans are counted as a single site, regardless of whether their alternative alleles match. At the allele level, distinct small variant alleles are identified by their chromosome, position, REF, and ALT, while SV alleles are distinguished using the PanSVMerger clusters. Variants are classified as KPR-specific, shared, or HPRC-specific according to which group carries them.

### Variant counting (per chromosome)

The `count_groups.py` script compares KPR and HPRC genotypes and counts SVs and small variants at both the locus and allele levels. A locus is shared when both groups carry an alternative allele at the same site, even if their alternative alleles differ.

### Summary

The `summarize_variant_comparison.py` script aggregates chromosome-level counts into summary tables for KPR-specific, shared, and HPRC-specific variants. It also calculates the novel percentage, the number of KPR-specific variants relative to the variants already present in HPRC:

```
novel % = KPR-specific / (shared + HPRC-specific) × 100
```

### Running the helper scripts independently

The helper scripts can also be run separately on intermediate results when only counting or summarization is needed; the full pipeline handles these steps automatically.

---

## 6. Output

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
  percentage (KPR-specific / (shared + HPRC-specific) × 100)
- `variant_comparison_by_class.tsv` — the same counts by variant class (SNV,
  insertion, deletion, complex). A complex locus holds alleles of more than one
  class, or multi-base substitutions
- `variant_comparison_shared_loci.tsv` — shared SV loci carrying a
  KPR-specific allele, and the number of such alleles
- `sv/CHR/CHR.sv.merged.vcf.gz` — SV loci after allele merging
- `small_variants/merged/CHR/CHR.small.vcf.gz` — merged small variants

---

## 7. References

1. Wu D, et al. (2026). Hidden genetic diversity in 320 nearly-complete East Asian genome assemblies. *bioRxiv*. https://doi.org/10.64898/2026.07.27.740875
2. Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nature Biotechnology* 42:663–673. https://doi.org/10.1038/s41587-023-01793-w
3. Liao WW, et al. (2023). A draft human pangenome reference. *Nature* 617:312–324. https://doi.org/10.1038/s41586-023-05896-x
4. Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nature Biotechnology* 36:875–879. https://doi.org/10.1038/nbt.4227
5. Garrison E, et al. (2022). A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. *PLOS Computational Biology* 18:e1009123. https://doi.org/10.1371/journal.pcbi.1009123
6. Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008. https://doi.org/10.1093/gigascience/giab008
7. Rognes T, et al. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584. https://doi.org/10.7717/peerj.2584
8. Cleary JG, et al. (2015). Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. *bioRxiv* 023754. https://doi.org/10.1101/023754

---

## Software

| Software | Version | Used for |
|---|---|---|
| vcfbub | 0.1.2 | bubble selection (https://github.com/pangenome/vcfbub) |
| vcflib | 1.0.15 | `vcfwave`, `vcfcreatemulti` |
| PanSVMerger | commit 9d9827f | SV allele clustering (https://github.com/Asian-Pan-Genome/PanSVMerger) |
| VSEARCH | 2.32.0 | sequence clustering called by PanSVMerger |
| RTG Tools | 3.13 | `rtg vcfdecompose` |
| bcftools / htslib | 1.24 | VCF filtering, normalization, merging, `bgzip` |
| samtools | 1.24 | `samtools faidx` |
| vg | 1.63.1 | `vg paths`, reference sequence from the graph |
| Python | 3.8 | helper scripts |
| pysam | 0.22.1 | `sv_merge_alleles.py` |

Bash is required. All executables must be available on `PATH`, and all helper
scripts are in the repository root.
