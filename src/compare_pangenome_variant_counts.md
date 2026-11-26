# 4.compare_pangenome_variant_counts.sh

Comparison of KPR-private, shared and HPRC-private graph variants at the locus and allele levels (Figures 3c–d; Table S21).

Part of the [KPR-pipeline](../README.md) repository.

## 1. Set up

### 1.1 Software

Install the programs under [Software and versions](#software-and-versions) and make sure they are on
`PATH`. The versions used in the study can be installed into one
mamba/conda environment (confirm the solver supports your platform):

```bash
mamba create -n kpr_variant_counts -c conda-forge -c bioconda \
  python=3.8 bcftools=1.24 htslib=1.24 samtools=1.24 vcflib=1.0.15 \
  vcfbub=0.1.2 vsearch=2.32.0 pysam=0.22.1 rtg-tools=3.13 vg=1.63.1
conda activate kpr_variant_counts
```

PanSVMerger [1] is not packaged and is used from its repository, at the
commit used in the study:

```bash
git clone https://github.com/Asian-Pan-Genome/PanSVMerger.git src/PanSVMerger
git -C src/PanSVMerger checkout 9d9827f
```

`sv_merge_alleles.py` imports the clustering functions of
`src/PanSVMerger/scripts/SV.merge_variants.py` unchanged. Set `PANSVMERGER` if the
repository is installed elsewhere. PanSVMerger is an external dependency
and is not bundled with this repository.

RTG Tools asks about usage logging the first time it runs, which can halt a
background job. Answer it once interactively, or set `RTG_TALKBACK=false` and
`RTG_USAGE=false` in the `rtg.cfg` of the installation, before running the
pipeline.

### 1.2 Input

Minigraph–Cactus can generate a graph-derived VCF during pangenome
construction using the `--vcf` option. This option is enabled in the
`construct` step of
[`3.pangenome_construction_and_statistics.sh`](pangenome_construction_and_statistics.md).
A GRCh38-referenced VCF can also be extracted from the completed graph
using `vg deconstruct -P GRCh38 -a`, as implemented in stage 3's
`variants` step. The latter step writes the default input files for
this workflow:

```text
results/
├── graph_variants/hprc_kpr/reference/hprc_kpr.GRCh38.vcf.gz   raw GRCh38-referenced VCF
└── index/hprc_kpr/hprc_kpr.xg                                 indexed pangenome graph
```

The input must be an unfiltered GRCh38-referenced graph VCF, with the
reference paths and nested variants required by this analysis. To use
the VCF generated directly by Minigraph–Cactus with `--vcf`, set
`INPUT_VCF` to its path and check that its reference and sample paths
and variant representation match the expected input. Do not apply
`vcfbub` beforehand; bubble selection is performed by `prepare_input`.

### 1.3 Reference files

In `REF_DIR`, downloaded from UCSC:

| File | Source |
|---|---|
| `hg38.centromeres.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz |
| `hg38.gap.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz |

Centromeres and the telomere entries of the gap table are excluded.

### 1.4 Population table

`LIST_DIR/population.tsv` holds one sample per line, sample ID and population
separated by a tab. The following identifiers are illustrative placeholders,
not sample names from the study. The sample IDs are the sample columns of the input VCF.
Samples of `TARGET_POP` (`KOR`) are compared with all other samples:

```
KOR_SAMPLE_01	KOR
HPRC_SAMPLE_01	CHN
HPRC_SAMPLE_02	AFR
```

### 1.5 Environment

```bash
export REF_DIR=/path/to/reference     # required
export OUT_DIR=/path/to/results       # default: ./results
export LIST_DIR=/path/to/lists        # default: ./lists
export JOBS=4                         # chromosomes or samples in parallel
```


## 2. Run

Print the available steps:

```bash
bash 4.compare_pangenome_variant_counts.sh
```

Run everything in order:

```bash
bash 4.compare_pangenome_variant_counts.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash 4.compare_pangenome_variant_counts.sh prepare_input
JOBS=23 bash 4.compare_pangenome_variant_counts.sh sv_decompose
JOBS=4 MERGE_WORKERS=16 bash 4.compare_pangenome_variant_counts.sh sv_merge
JOBS=24 bash 4.compare_pangenome_variant_counts.sh small_decompose
JOBS=12 bash 4.compare_pangenome_variant_counts.sh small_merge
JOBS=12 bash 4.compare_pangenome_variant_counts.sh count
bash 4.compare_pangenome_variant_counts.sh summary
```

Test on one chromosome first. With a separate `OUT_DIR`, point `INPUT_VCF`
and `GRAPH_XG` to the existing input, since their defaults are under `OUT_DIR`:

```bash
CHROMS=chr22 OUT_DIR=test \
INPUT_VCF=results/graph_variants/hprc_kpr/reference/hprc_kpr.GRCh38.vcf.gz \
GRAPH_XG=results/index/hprc_kpr/hprc_kpr.xg \
bash 4.compare_pangenome_variant_counts.sh all
```

SV steps run per chromosome, `small_decompose` per sample. `sv_decompose` runs
one single-threaded `vcfwave` per chromosome and needs up to 35 GB each;
`sv_merge` uses `MERGE_WORKERS` cores per chromosome. Each job writes its own
log to `$OUT_DIR/logs/variant_comparison/`, and every step writes to
`$OUT_DIR/logs/pipeline.log`. If any job of a step fails, the step exits with
an error once its remaining jobs have finished, and the following steps are not
run.


## 3. Steps

| Step | Description |
|---|---|
| `prepare_input` | Run `vcfbub`, rename/reference-filter contigs and samples, then exclude centromere/telomere-overlapping variants; recalculate AC/AN and prepare the SV input and reference sequence |
| `sv_decompose` | `vcfwave` decomposition, SV alleles, including reassigned inversions and long alleles, SV loci with `vcfcreatemulti` |
| `sv_merge` | PanSVMerger allele clustering within each SV locus |
| `small_decompose` | Per sample: `vcfwave`, `bcftools norm`, `rtg vcfdecompose`, small variants of 50 bp or less |
| `small_merge` | Per chromosome: sample merge, SNV and indel files |
| `count` | KPR-private, shared and HPRC-private variants per chromosome |
| `summary` | Summary tables at the locus and allele level |
| `all` | All of the above, in order |

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
| `SV_MIN_LEN` | `50` | Length-difference threshold; inversions and long alleles are reassigned |
| `MERGE_IDENTITY` | `0.9` | vsearch identity used by PanSVMerger |
| `MERGE_WORKERS` | `16` | parallel vsearch processes per chromosome |
| `MULTI_MEM_GB` | `32` | memory cap of `vcfcreatemulti` |
| `RTG_MAXLEN` | `5000` | records with a longer allele, and inversions, are not passed to rtg |
| `RTG_MEM` | `8g` | rtg heap size |
| `SORT_MEM` | `8G` | `bcftools sort` buffer |
| `PANSVMERGER` | `./src/PanSVMerger/scripts/SV.merge_variants.py` | PanSVMerger script |
| `JOBS` | `4` | chromosomes or samples processed in parallel |

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

- `variant_comparison_summary.tsv` — KPR-private, shared and HPRC-private
  SVs and small variants at the locus and allele levels, with the
  proportion of KPR-private variants among all three classes (%)
- `variant_comparison_by_class.tsv` — the same counts by variant class (SNV,
  insertion, deletion, complex). A complex locus holds alleles of more than one
  class, or multi-base substitutions
- `variant_comparison_shared_loci.tsv` — shared SV loci carrying a
  KPR-private allele, and the number of such alleles
- `sv/CHR/CHR.sv.merged.vcf.gz` — SV loci after allele merging
- `small_variants/merged/CHR/CHR.small.vcf.gz` — merged small variants

## 6. Graph variant extraction and filtering

We extracted variants from HPRC+KPR using vg deconstruct (-P GRCh38 -a). The analysis retained chr1–22 and chrX and excluded the CHM13 reference haplotype. vcfbub v0.1.2 (-l 0 -r 10000000) selected nested bubbles in place of top-level reference alleles longer than 10 Mb. After reference/sample selection, UCSC hg38 centromere and telomere intervals were excluded by record overlap. AC and AN were recalculated and records without observed ALT alleles removed.

## 7. Small variant decomposition

For each sample, vcfwave (-L 100000) decomposed complex alleles and bcftools norm v1.24 (-c s -f GRCh38 -m -any) normalized and split the records. Inversions and alleles longer than 5 kb were excluded from small variant counts and reassigned to the SV branch. RTG Tools v3.13 vcfdecompose (--break-mnps --break-indels) decomposed remaining MNPs and complex indels. After retaining length differences ≤50 bp, dedup_sample.py consolidated duplicate genotypes. bcftools merge (-m none --missing-to-ref) merged the decomposed per-sample VCFs; AC, AN and AF were recalculated. This missing-to-reference setting applies to the decomposed graph comparison, not to Mendelian analysis.

## 8. Structural variant allele clustering

The joint multi-sample VCF was decomposed using vcfwave (-L 100000) and bcftools norm (-m-). The SV branch retained length differences >50 bp and the reassigned inversion/long-allele records; vcfcreatemulti merged overlapping reference intervals. PanSVmerger (commit 9d9827f) and VSEARCH v2.32.0 (--cluster_fast --id 0.9 --strand both) clustered similar alleles, splitting clusters at allele-length gaps greater than 50 bp. The wrapper preserves genotypes and uses one representative sequence per cluster; alleles clustering with REF become reference alleles. Loci with no non-reference cluster are removed. After genotype reassignment AC, AN and AF are updated. Retain the clustering logs and complex variant classes when comparing the SV outputs with the reported INS and DEL counts.

## 9. Locus-level and allele-level comparison

A variant is carried by a population panel when at least one called haplotype contains its ALT allele. Missing genotypes are ignored in panel presence counting. Small-variant loci are the union of overlapping REF intervals; small-variant alleles are identical CHROM, POS, REF and ALT records. SV loci are merged records from vcfcreatemulti; SV alleles are PanSVmerger clusters. Count KPR-private, shared and HPRC-private loci and alleles separately; count KPR-private SV alleles inside shared loci as an additional category.

For the study, the proportion private within KPR is KPR-private / (KPR-private + shared). The proportion added relative to HPRC is KPR-private / (HPRC-private + shared). The generic summary's union-based proportion has denominator KPR-private + shared + HPRC-private. Use the first two ratios for the Results and Table S21.

## References

1. Wu D, et al. (2026). Hidden genetic diversity in 320 nearly-complete East Asian genome assemblies. *bioRxiv*. https://doi.org/10.64898/2026.07.27.740875
2. Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nature Biotechnology* 42:663–673. https://doi.org/10.1038/s41587-023-01793-w
3. Liao WW, et al. (2023). A draft human pangenome reference. *Nature* 617:312–324. https://doi.org/10.1038/s41586-023-05896-x
4. Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nature Biotechnology* 36:875–879. https://doi.org/10.1038/nbt.4227
5. Garrison E, et al. (2022). A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. *PLOS Computational Biology* 18:e1009123. https://doi.org/10.1371/journal.pcbi.1009123
6. Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008. https://doi.org/10.1093/gigascience/giab008
7. Rognes T, et al. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584. https://doi.org/10.7717/peerj.2584
8. Cleary JG, et al. (2015). Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. *bioRxiv* 023754. https://doi.org/10.1101/023754


## Software and versions

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

Bash is required. All executables must be available on `PATH`, and the supporting Python scripts are under `src/`.
