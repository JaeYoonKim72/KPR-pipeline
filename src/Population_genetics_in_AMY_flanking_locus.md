# 8.Population_genetics_in_AMY_flanking_locus.sh

Population genetics of the single-copy distal and proximal flanking regions of the AMY locus (Figure 6; Tables S30–S32).

Part of the [KPR-pipeline](../README.md) repository.

## 1. Set up

The two GRCh38 flanks are chr1:103,456,163–103,561,526 (105,364 bp; bundle 0) and chr1:103,760,698–103,826,698 (66,001 bp; bundle 1a). Exclude the intervening SVR. The study used the GRCh38 StrictMask (20160622) intersected with non-N positions, giving callable lengths of 90,193 and 59,429 bp.

Use an indexed joint diploid biallelic SNP VCF for diversity, LD and FST, a separately prepared phased joint VCF for XP-nSL, and disjoint sample lists. SNP data for KOR were taken before missingness and MAF filtering, phasing and imputation; genotype DP <10 or GQ <20 was set to missing. Six KOR individuals with inbreeding coefficients beyond three scaled MADs from the median were excluded; the retained 302 individuals had no pair above KING-robust kinship 0.0884. Comparison data were the 2,504 unrelated Phase 3 individuals from the 30× GRCh38 1000 Genomes dataset.

| Population | Individuals |
|---|---:|
| KOR | 302 |
| CHN (CHB + CHS) | 208 |
| CDX | 93 |
| JPN | 104 |
| KHV | 99 |
| SAS | 489 |
| EUR | 503 |
| AMR | 347 |
| AFR | 661 |

Within each population, retain biallelic SNPs with call rate ≥90%, one REF–ALT record per position, and recalculate AC/AN using non-missing genotypes. Apply the callable mask. For XP-nSL use the MAF ≥0.05-filtered, SHAPEIT5-phased and IMPUTE5-imputed KOR dataset and phased 1000 Genomes data. Match sites by POS, REF and ALT.

~~~bash
mamba create -n amy_popgen -c conda-forge -c bioconda python bcftools bedtools htslib selscan=3.0.1
conda activate amy_popgen
export VCF=/path/to/joint.qc.biallelic.snps.vcf.gz
export XPNSL_VCF=/path/to/joint.phased.biallelic.snps.vcf.gz
export GENOME_VCF="$VCF"
export CALLABLE_BED=/path/to/GRCh38.StrictMask.nonN.bed.gz
export POPS="KOR=/path/to/KOR.txt CHN=/path/to/CHN.txt CDX=/path/to/CDX.txt JPN=/path/to/JPN.txt KHV=/path/to/KHV.txt SAS=/path/to/SAS.txt EUR=/path/to/EUR.txt AMR=/path/to/AMR.txt AFR=/path/to/AFR.txt"
export OUT_DIR=/path/to/results/amy_flanking
~~~

## 2. Run

~~~bash
bash 8.Population_genetics_in_AMY_flanking_locus.sh all
~~~

A single population runs within-population analyses; with multiple populations the first population is the query for FST and XP-nSL. Run the study's individual and genomic resampling settings as follows:

~~~bash
INDIVIDUAL_CI=1 IND_BOOT_REPS=10000 IND_JACKKNIFE_FRACTION=1 \
GENOMIC_BOOTSTRAP=1 GENOMIC_CHROMS="$(printf 'chr%s,' {1..22} | sed 's/,$//')" \
GENOMIC_BLOCK_MB=5 GENOMIC_BOOT_REPS=10000 GENOMIC_WINDOW_STEP=1000 \
GENOMIC_FST_CALLABLE_MATCH=-1 \
bash 8.Population_genetics_in_AMY_flanking_locus.sh all
~~~

The wrapper calculates focal statistics, raw phased XP-nSL, individual resampling and genomic ranks for diversity and FST. Sections 7 and 9 describe the separate normalization, subsampling, LD background and XP-nSL genomic analyses. The independent nsl command is available for other studies; within-population nSL was not reported in this study.

## 3. Steps

| Step | Primary statistic | Required population count | Main result |
|---|---|---|---|
| `diversity` | π, projected θW, projected Tajima's D | ≥1 | `stats/diversity.tsv` |
| `ld` | Mean within-flank corrected dosage r² | ≥1 | `stats/within_ld.tsv` |
| `across_ld` | Mean distal × proximal corrected dosage r² | ≥1 | `stats/across_ld.tsv` |
| `fst` | Hudson FST | ≥2 | `stats/hudson_fst.tsv` |
| `nsl` | Optional stand-alone unphased nSL (not in all) | ≥1, selscan | `selscan/` |
| `xpnsl` | Phased XP-nSL raw core scores | ≥2, selscan | `selscan/` |
| `individual_ci` | Bootstrap and jackknife for SNP-derived distal/proximal flanking-region statistics | ≥1 | `individual_ci/individual_ci.tsv` |
| `genomic` | 1-kb genomic background windows, block-bootstrap genomic percentile | ≥1; genomewide VCF | `genomic/genomic_percentile_ci.tsv` |
| `all` | Diversity, LD, FST, phased XP-nSL; optional individual CI and genomic ranks | ≥1 | `status.tsv` |

## 4. Configuration

| Variable / option | Default | Description |
|---|---|---|
| `VCF` / `--vcf` | required | Joint diploid biallelic SNP VCF and index |
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
| `NSL_CORE_MAF` / `--core-maf` | `0.05` | Core MAF for stand-alone nSL; XP-nSL uses post-scoring balanced MAF |
| `NSL_MAX_EXTEND` / `--max-extend-nsl` | `50` | Maximum SNP markers for nSL/XP-nSL extension |
| `THREADS` | `4` | selscan threads (no outer parallelism) |
| `INDIVIDUAL_CI` / `--individual-ci` | `0` (OFF) | Enable individual-level resampling CIs |
| `IND_BOOT_REPS` / `--bootstrap-reps` | `100` | Individual bootstrap replicates; paper used 10,000 |
| `IND_JACKKNIFE_FRACTION` / `--jackknife-fraction` | `0.30` | ceil(30% N) deleted-person replicates/group; paper used 1.0 |
| `GENOMIC_BOOTSTRAP` / `--genomic` | `0` (OFF) | Enable genomewide windows and block bootstrap |
| `GENOMIC_CHROMS` / `--chromosomes` | `chr1` | Configured background chromosomes; study chr1–chr22 |
| `GENOMIC_BLOCK_MB` / `--block-mb` | `50` | Physical block size; study 5 Mb |
| `GENOMIC_WINDOW_STEP` / `--window-step` | `1000` | Sliding-window step in bp; paper also used 1 kb |
| `GENOMIC_BOOT_REPS` / `--genomic-reps` | `100` | Genome block draws; paper used 10,000 |
| `GENOMIC_FST_CALLABLE_MATCH` / `--fst-callable-match` | `0.10` | Optional callable-bp matching; study -1 (disabled) |
| `GENOME_VCF` / `--genome-vcf` | `VCF` | Genomewide joint VCF for configured chromosomes |
| `GENOME_CHROM_SIZES` / `--chrom-sizes` | built-in GRCh38 | Optional reference `.fai` to override chromosome lengths |
| `SEED` / `--seed` | `20261117` | Deterministic resampling seed |
| `XPNSL_VCF` / `--xpnsl-vcf` | `VCF` | Separate phased VCF for haplotype scoring |

## 5. Output

~~~text
results/amy_flanking_population_genetics/
├── stats/
│   ├── diversity.tsv             π, θW, Tajima's D
│   ├── within_ld.tsv             within-flank corrected r²
│   ├── across_ld.tsv             distal × proximal corrected r²
│   └── hudson_fst.tsv            first-vs-other Hudson FST
├── individual_ci/                only when INDIVIDUAL_CI=1
│   ├── individual_ci.tsv         percentile/shifted bootstrap and JK intervals
│   └── individual_replicates.tsv person-resampling replicate values
├── genomic/                      only when GENOMIC_BOOTSTRAP=1
│   └── genomic_percentile_ci.tsv upper-tail ranks and block bootstrap 95% interval
├── intermediate/                 region VCFs used for selscan
├── selscan/                      raw XP-nSL and optional nSL outputs (if installed)
├── logs/
└── status.tsv                    DONE, RAW_ONLY, SKIPPED, FAILED
~~~

Correspondence with paper supplementary tables:

| Table | Within/between | Current public-code distal/proximal flanking-region coverage | Additional published metrics |
|---|---|---|---|
| S30 | Within-population diversity | π, θW and exploratory Tajima D | Projected Tajima D and genomewide ranks/CIs for π, θW and Tajima D |
| S31 | Within- and cross-flank LD | Within-LD and across-SVR corrected r² | Distance-matched LD enrichment and genomewide LD rank/interval |
| S32 | Between-population comparison | Hudson FST and raw XP-nSL | Normalized phased XP-nSL, empirical 93:93 subsampling, complete ranks/CIs |

`status.tsv` distinguishes successfully computed values from raw,
skipped, or failed stages. **Raw XP-nSL is not equivalent to the
normalized phased XP-nSL reported in Table S32.** Unphased within-population
nSL is an optional exploratory output and is **not reported** in the
final supplementary Tables S30–S32.

## 6. Within-population statistics

### 6.1 Nucleotide diversity, Watterson's theta and Tajima's D

Retain polymorphic SNPs (0 < AC < AN) without additional MAF filtering. At each site pi_i = 2 AC_i (AN_i − AC_i) / [AN_i (AN_i − 1)]. Sum pi_i and divide by the callable bases L. No projection is applied to pi. Project the site-frequency spectrum to 168 allele copies using hypergeometric sampling, giving E[S_168]. Calculate thetaW = E[S_168] / (a_167 L), where a_167 is the 167th harmonic number. Calculate Tajima's D from the observed pi sum and E[S_168] with the variance coefficients for 168 alleles. The projection corresponds to 84 diploid individuals at the ≥90% call-rate threshold in CDX. Report pi, thetaW and Tajima's D separately for each flank (Table S30).

### 6.2 Within-region LD

Use biallelic SNPs with MAF ≥0.05. For each SNP pair calculate the squared Pearson correlation of unphased genotypes (0, 1 or 2), using individuals called at both sites. Correct for sample size: r²_corrected = [(n − 1) r²_raw − 1] / (n − 2), retaining negative values. Within-region LD is mean corrected r² for pair distances <100 kb distally and <60 kb proximally (Table S31).

### 6.3 Cross-SVR LD

Use one SNP from each flank (possible distances 199,172–370,535 bp). Compute the mean corrected r² over all eligible pairs. Compare with distance-matched genomic means, as described in section 9.2.

## 7. Between-population statistics

### 7.1 Hudson's FST

For KOR versus each comparison population, retain SNPs with matching POS/REF/ALT, the population call-rate filters and callable mask. No additional MAF filter is applied. The study retained 166 distal and 128 proximal SNPs. For frequencies p and q and called allele numbers n and m, N_i = (p−q)² − p(1−p)/(n−1) − q(1−q)/(m−1) and D_i = p(1−q)+q(1−p). Hudson's FST is sum(N_i)/sum(D_i), retaining negative estimates. Omit sites with D_i=0 or insufficient called alleles. For individual resampling, hold the observed SNP panel fixed and recalculate the allele counts and estimator.

### 7.2 Phased XP-nSL

We used selscan v3.0.1 with KOR as query, each comparison population as reference and the two flanks analyzed separately. Run XP-nSL with --keep-low-freq and --max-extend-nsl 50, allowing low-frequency markers in the phased input to contribute to extension. The scoring command uses --maf 0 to avoid the pre-scoring frequency cutoff; select valid scored cores after calculation. For each matched SNP calculate the allele frequency in each population, let p_bal = (p_KOR + p_reference)/2 and MAF_bal = min(p_bal, 1−p_bal). Retain scored cores with MAF_bal ≥0.05. Apply this filter to focal and background score files using their matched phased VCFs and retain the population-pair order.

### 7.3 XP-nSL normalization

For each pair, score chromosomes 2–22 with the same phasing mode, marker extension and balanced core rule. Estimate the pair-specific mean and SD from valid background cores and standardize scores as z_i = (raw_i − mean_bg)/sd_bg. Apply that normalization to the focal raw scores. For selscan norm, use the independently prepared, balanced-core-filtered background score files:

~~~bash
mapfile -t bg < <(find /path/to/pair/chr2_22 -name '*.xpnsl.out' | sort)
selscan norm --xpnsl --log pair.xpnsl.norm.log --files "${bg[@]}"
selscan norm --xpnsl --log-input pair.xpnsl.norm.log --files distal.xpnsl.out proximal.xpnsl.out
~~~

Use the matching background normalization for each ordered pair. Report number of valid scored cores, arithmetic mean signed z, mean |z|, fraction z ≥2 and fraction z ≤−2. Positive signed XP-nSL indicates longer haplotype homozygosity in KOR and negative values indicate longer homozygosity in the comparison population (Table S32; Figure 6g).

### 7.4 Equal-sample-size subsampling

Sample 93 individuals without replacement from each population, independently for 100 replicates. For CDX all 93 individuals enter every replicate and only the KOR subset varies. Recalculate phased XP-nSL, select valid balanced cores, and use the corresponding background normalization for the same sample-size analysis. Summarize the median and 2.5th/97.5th percentiles of replicate mean signed XP-nSL. These are empirical subsampling intervals, not 95% confidence intervals. Keep the sample lists, scoring/normalization files and replicate summaries separately.

## 8. Individual bootstrap and jackknife

### 8.1 Individual bootstrap

Resample individuals with replacement to the observed sample size (10,000 replicates). For FST, resample the two populations independently while retaining the original SNP panel. Calculate the 2.5th and 97.5th percentiles of the replicate estimates. The wrapper also reports mean-shifted intervals as an additional statistic; the study reports the raw percentile interval.

### 8.2 Delete-one jackknife

Exclude each individual once, holding the other population fixed for FST. For two-population FST sum the two group-specific jackknife variance components and calculate SE. Report the observed estimate ±1.96 × SE. The study used complete delete-one jackknife intervals as the primary uncertainty measure for thetaW, Tajima's D and FST because individual bootstrap distributions were shifted relative to observed estimates. The generic fractional jackknife option is a separate reduced-replicate calculation.

## 9. Genome-wide comparison and block bootstrap

### 9.1 Genomic sliding windows

Slide physical windows of 105,364 bp and 66,001 bp along chr1–22 with a 1-kb step, excluding windows overlapping the AMY locus and those with insufficient data. Apply the same genotype, allele and callable criteria to the background and focal regions. For the study the window lengths were matched but the callable-base lengths and SNP numbers were not. Set GENOMIC_FST_CALLABLE_MATCH=-1. Use chr1 background medians for pi, thetaW and Tajima's D, and autosomal background medians for within-region LD and FST. Calculate medians from the eligible window values.

### 9.2 Genome-wide LD enrichment

**Within-flank LD.** For each population, calculate mean sample-size-corrected r² in sliding autosomal windows of the corresponding flank length, with a 1-kb step, the same ≥90% call-rate and MAF ≥0.05 filters, and pair separations <100 kb (distal) or <60 kb (proximal). Exclude the AMY locus. Compare the focal mean to the autosomal window median and calculate the genomic rank. For the within-region background, retain the window means and their physical block assignments.

**Cross-SVR LD.** Use background flank pairs preserving the distal and proximal window lengths and their relative spacing. Bin eligible SNP-pair distances into 1-kb bins. Within each population let w_k be the fraction of eligible AMY pairs in distance bin k, and B_k the mean corrected genomic r² in that bin. The distance-matched background B = sum(w_k B_k). Calculate B separately for chromosome 1 and for the autosomes, using the same sample, SNP, callable and corrected-r² filters. If F is the mean AMY cross-SVR corrected r², enrichment = F/B. For Figure 6b, the error-bar endpoints are the 2.5th and 97.5th percentiles of the individual-bootstrap AMY r² divided by the fixed chr1 background mean B. The background denominator is not resampled for this plotted interval.

The focal estimate is provided in stats/across_ld.tsv; calculate the distance-bin background separately and save a TSV with population and matched_mean_corrected_r2 columns. The following command combines the focal and background outputs with the bootstrap interval from individual_ci/individual_ci.tsv:

~~~bash
python3 - "$OUT_DIR/stats/across_ld.tsv" /path/to/chr1_matched_ld.tsv "$OUT_DIR/individual_ci/individual_ci.tsv" <<'PY'
import csv
import sys

def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

bg = {r["population"]: float(r["matched_mean_corrected_r2"]) for r in rows(sys.argv[2])}
ci = {r["population"]: r for r in rows(sys.argv[3]) if r["statistic"] == "across_corrected_r2"}
print("population\tfocal_r2\tmatched_r2\tfold\tfold_025\tfold_975")
for r in rows(sys.argv[1]):
    pop = r["population"]
    b = bg[pop]
    if b <= 0:
        raise ValueError("Non-positive background: " + pop)
    f = float(r["mean_corrected_r2"])
    c = ci[pop]
    print(pop, f, b, f/b, float(c["bootstrap_raw_025"])/b,
          float(c["bootstrap_raw_975"])/b, sep="\t")
PY
~~~

The genomewide LD window and matched distance-bin computations are separate from amy_flank_genomic_background.py, which calculates diversity and FST window ranks.

### 9.3 Upper-tail genomic rank

For focal estimate T and W eligible background windows, upper-tail rank (%) = 100 × [number of windows greater than T + 0.5 × number tied with T] / W. Rank is based on Hudson's FST, mean |XP-nSL| or the corresponding within-population statistic. Smaller ranks indicate a higher relative position. Calculate autosomal ranks and supplementary chr1 ranks from the corresponding window distributions. The wrapper compares values at the same reported numeric precision.

### 9.4 Chromosome-stratified block bootstrap

Partition each autosome into contiguous 5-Mb physical blocks and assign windows to blocks by their midpoint. In each of 10,000 replicates, resample the same number of physical blocks within each chromosome with replacement. Pool the selected blocks' eligible windows and recalculate the rank; report the 2.5th and 97.5th percentiles, valid replicate count, occupied blocks and total blocks. The length-matched windows overlap and the ranks are used for relative genomic comparison.

### 9.5 Genome-wide XP-nSL comparison

Use the autosomal phased, pair-specific normalized XP-nSL scores and the balanced post-scoring MAF criterion in section 7. Calculate the mean absolute normalized score in each 1-kb-step window matching the corresponding flank length. Retain only windows with at least 10 valid scored core SNPs. Apply the upper-tail rank and 5-Mb block bootstrap described above. The raw scoring step supplies the flank scores; genomewide scoring, pair-specific normalization, window means and 93:93 subsampling are separate analyses of the phased inputs.

## References

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


## Software and versions

| Software | Version | Operation |
|---|---|---|
| bcftools | 1.15+ | VCF preparation |
| PLINK | 2.0 | Inbreeding and relatedness QC |
| SHAPEIT5 | 5.1.1 | KOR phasing |
| IMPUTE5 | study input | KOR imputation |
| selscan | 3.0.1 | Phased XP-nSL and normalization |
| Python | 3.9+ | SNP statistics, resampling and genomic ranks |
| bedtools, htslib | compatible versions | Mask and indexed VCF preparation |
