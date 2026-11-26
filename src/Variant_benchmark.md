# 5.Variant_benchmark.sh

Genotype-aware benchmarking of GIAB and KOR variant calls and trio inheritance inconsistency analysis (Figures 4a–e; Tables S24–S26).

Part of the [KPR-pipeline](../README.md) repository.

## 1. Set up

Use the indexed GRCh38 FASTA, matching RTG SDF, per-sample confidence BEDs and variant call sets in GRCh38 coordinates. Install bcftools, samtools, hap.py, RTG Tools, Truvari v5.5.0 and pysam (Python ≥3.12 for Truvari).

~~~bash
export REF_DIR=/path/to/reference
export TRUTH_DIR=/path/to/truth_input
export QUERY_DIR=/path/to/query_input
export CONF_BED_DIR=/path/to/confidence_beds
export SAMPLE_LIST=/path/to/benchmark_samples.list
export OUT_DIR=/path/to/benchmark_results
samtools faidx "$REF_DIR/GRCh38.fa"
rtg format -o "$REF_DIR/GRCh38.sdf" "$REF_DIR/GRCh38.fa"
~~~

For KOR, use the seven samples KOR-001, KOR-002, KOR-006, KOR-007, KOR-008, KOR-009 and KOR-010. Input layouts:

~~~text
truth_input/
  hifi_deepvariant/SAMPLE.hifi_dv.biSNP.addfilter.VAF0208.vcf.gz
  hifi_deepvariant/SAMPLE.hifi_dv.Small.vcf.gz
  long_read_sv/svim_asm_SAMPLE.PASS.SVLEN51.chr.vcf.gz
  long_read_sv/pbsv_SAMPLE.PASS.SVLEN51.chr.vcf.gz
  long_read_sv/sniffles_SAMPLE.PASS.SVLEN51.chr.vcf.gz
  long_read_sv/svim_SAMPLE.PASS.SVLEN51.chr.vcf.gz
  giab_small/HG001.vcf.gz ... HG007.vcf.gz
confidence_beds/
  SAMPLE.small.bed
  SAMPLE.sv.bed
query_input/
  kpr/kpr.SAMPLE.snp.vcf.gz
  kpr/kpr.SAMPLE.indel.vcf.gz
  kpr/kpr.SAMPLE.sv.INS.vcf.gz
  kpr/kpr.SAMPLE.sv.DEL.vcf.gz
  hprc/...
  hprc_kpr/...
  deepvariant/deepvariant.SAMPLE.snp.vcf.gz
  deepvariant/deepvariant.SAMPLE.indel.vcf.gz
  manta/manta.SAMPLE.sv.INS.vcf.gz
  delly/delly.SAMPLE.sv.DEL.vcf.gz
~~~

The benchmark uses the actual high-confidence BED for each truth set. Prepare the Manta INS/DUP input as the insertion call set used in the benchmark, with normalized insertion representation and genotype information. Optional Manta DEL inputs can be compared separately.

## 2. Run

~~~bash
export TRUTH_MODE=kor
bash 5.Variant_benchmark.sh check
bash 5.Variant_benchmark.sh all
~~~

For GIAB v4.2.1 small-variant truth files, use a sample list containing HG001–HG007. Set GIAB_SMALL_DIR to the directory containing one truth VCF per sample. Supply HG002 T2T-HG002-Q100 v1.1 SV truth for the HG002 SV benchmark and optionally the GRCh38 LRS query INS/DEL VCFs.

~~~bash
export TRUTH_MODE=giab
export SAMPLE_LIST=/path/to/GIAB7.list
export GIAB_SMALL_DIR=/path/to/giab_small
export GIAB_SV_VCF=/path/to/HG002.Q100.GRCh38.vcf.gz
export GIAB_LRS_SV_DIR=/path/to/giab_lrs_sv
export OUT_DIR=/path/to/giab_results
bash 5.Variant_benchmark.sh all
~~~

GIAB_LRS_SV_DIR contains HG002.sv.INS.vcf.gz and HG002.sv.DEL.vcf.gz. Omit it for SRS-only evaluation. GIAB_SMALL_VCF can instead specify a single-sample truth file. Run KOR and GIAB in separate output directories.

For the trio workflow, configure the child, parental and confidence-region inputs and run the mendelian step independently (section 10).

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
| `ALLOW_AUTOSOME_AS_CONFIDENCE` | 0 | Whole-autosome comparison option; study uses truth confidence BEDs |
| `GBZ_EVIDENCE_DIR` | Unset | Independently prepared Truvari results for personalized-GBZ allele representation |
| `JOBS` / `THREADS` | 8 / 4 | Parallel jobs and hap.py threads |
| `MEND_TRIO` / `MEND_ARM` | AJ / graph | Trio and reference arm |
| `MEND_MIN_GQ` / `MEND_MIN_DP` | 20 / 10 | Mendelian genotype filter; linear SV DP unavailable |
| `GIAB_SMALL_DIR` | `$TRUTH_DIR/giab_small` | Sample-specific v4.2.1 VCFs |
| `GIAB_SMALL_VCF` | unset | Single-sample truth VCF override |
| `GIAB_SV_VCF` | unset | HG002 Q100 SV truth VCF |
| `GIAB_LRS_SV_DIR` | unset | HG002 GRCh38 LRS query INS/DEL files |

## 5. Output

The summary filenames retain the `kor_` prefix in both truth modes. Use the individual and mean tables for S24 and S26 and the trio summary/adjudication tables for S25.

## 6. Construction of KOR truth sets

PacBio LRS DeepVariant v1.6 calls were restricted to autosomal PASS SNPs and 1–50-bp indels with GQ ≥20, DP ≥7, heterozygous VAF 0.30–0.70 and homozygous ALT VAF ≥0.90. The SV truth procedure used assembly SVIM-asm and three read-based callers, the filtering criteria in stage 1, and Truvari support (500 bp, ≥70% size and sequence similarity). Truth SVs required ≥5 ALT reads and DP ≥10 in at least one matching read caller. We generated ≥1, ≥2 and three-caller-supported truth sets, with ≥2 as the primary comparison.

## 7. GIAB benchmarking

GIAB v4.2.1 small variants were evaluated for HG001–HG007. HG002 SVs were evaluated against T2T-HG002-Q100 v1.1 in GRCh38 coordinates, including GRCh38-based SRS and LRS queries and three graph-based SRS query panels. The HG002 SV truth denominator after PASS, spanning-deletion and 51–50,000-bp size filtering was 32,503 (20,815 INS and 11,688 DEL). The official high-confidence BEDs were applied to the comparisons.

## 8. Variant benchmarking and metrics

Normalize multiallelic calls with bcftools norm -m -any, retain autosomal PASS non-reference genotypes, and classify SNP, small INS/DEL (1–50 bp) and SV INS/DEL (51–50,000 bp). hap.py used RTG vcfeval and the corresponding high-confidence BED. Truvari bench v5.5.0 used the following comparison parameters:

~~~bash
truvari bench -b truth.vcf.gz -c query.vcf.gz -o bench \
  --includebed confidence.bed --refdist 500 --pctsize 0.8 \
  --pctseq 0.8 --pctovl 0.8 --pick single --sizemin 51 \
  --sizefilt 30 --sizemax 50000 --passonly
~~~

TPs require allele and genotype concordance. For small variants, TP = TRUTH.TP; presence = (TRUTH.TP + BK=am)/TRUTH.TOTAL; genotype concordance = TRUTH.TP/(TRUTH.TP + BK=am). Precision, recall and F1 use METRIC.Precision, METRIC.Recall and METRIC.F1_Score from hap.py summary.csv.

For SVs, presence = TP-base/base cnt; genotype concordance = gt_concordance. Recall = TP-base_TP-gt/base cnt; precision = TP-comp_TP-gt/comp cnt. F1 is the harmonic mean. GT-discordant truth-side and query-side matches contribute to FNs and FPs respectively (TP-base_FP-gt and TP-comp_FP-gt). Summarize KOR seven-sample means at the individual level.

## 9. False-negative classification

Small-variant FN annotations BK=. indicate detection failures; BK=lm and BK=am indicate allele and genotype matching failures. For SVs, inspect unmatched truth variants in fn.vcf.gz and query calls within ±500 bp. For each personalized GBZ used in mapping, deconstruct the relevant graph in GRCh38 coordinates and compare its represented alleles to the truth with the same Truvari matching parameters. A truth allele absent from that graph, or a truth call without any query call within the positional tolerance, is a detection failure. Remaining unmatched calls with a corresponding query call are matching failures. Supply the separately generated graph-representation evidence to GBZ_EVIDENCE_DIR as specified by classify_benchmark_fn.py; without it the output records the positional evidence, not a complete graph-representation adjudication.

## 10. Mendelian genotype inheritance inconsistency

The AJ trio contains HG002 (child), HG003 (father), HG004 (mother); the HC trio contains HG005 (child), HG006 (father), HG007 (mother). Use the leave-one-sample-out graphs described in stage 6 for the relevant child. DeepVariant trio gVCFs require GQ ≥20 and DP ≥10 in every member. SVs require PASS, 51–50,000 bp and GQ ≥20; graph SVs additionally require DP ≥10, and Delly precise calls are used for the linear DEL comparison. Merge genotype VCFs with bcftools merge, preserving genuinely called 0/0 and leaving missing genotypes missing. RTG mendelian annotates incompatible sites with INFO/MCV.

~~~bash
export MEND_TRIO=AJ
export MEND_ARM=graph
export MEND_REFERENCE=kpr
export MEND_VCF_DIR=/path/to/trio_genotype_vcfs
export MEND_TRUTH_SMALL=/path/to/HG002.small.truth.vcf.gz
export MEND_CONF_SMALL=/path/to/HG002.small.bed
export MEND_TRUTH_SV=/path/to/HG002.sv.truth.vcf.gz
export MEND_CONF_SV=/path/to/HG002.sv.bed
bash 5.Variant_benchmark.sh mendelian
~~~

For linear comparisons set MEND_ARM and the Manta INS/Delly DEL input directories according to trio_mendelian_inconsistency.sh. The observed rate is the number of incompatible autosomal child non-reference genotypes divided by all testable diploid trio genotypes at child non-reference sites. Truth-set adjudication uses the child allele and genotype within confidence regions: concordant calls are true inconsistent, absent or discordant child calls false inconsistent. The false rate divides false by adjudicated observed inconsistencies. Figure 4c reports the AJ trio, excluding the GRCh38 SV bars; the HC trio has small-variant truth but no corresponding SV truth set (Table S25).

## 11. Statistical comparison

For the seven KOR samples, compare paired sample-level SV recall values for the three pangenome references. Use two-sided exact paired Wilcoxon signed-rank tests and Holm-adjust the three pairwise P values. The combined SV recall is calculated from summed genotype-concordant truth counts and truth denominators for INS and DEL within each sample before the paired comparisons.

## References

1. Wagner J, Olson ND, Harris L, et al. (2022). Benchmarking challenging small variants with linked and long reads. *Cell Genomics*. https://doi.org/10.1016/j.xgen.2022.100128
2. English AC, Menon VK, Gibbs RA, Metcalf GA, Sedlazeck FJ. (2022). Truvari: refined structural variant comparison preserves allelic diversity. *Genome Biology*. https://doi.org/10.1186/s13059-022-02840-6
3. Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*. https://doi.org/10.1093/gigascience/giab008
4. Cleary JG, Braithwaite R, Gaastra K, et al. (2015). Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. *bioRxiv*. https://doi.org/10.1101/023754
5. Wu D, Yang C, Chen Q, et al. (2026). Hidden genetic diversity in 320 nearly-complete East Asian genome assemblies. *bioRxiv*. https://doi.org/10.64898/2026.07.27.740875
6. Suzuki Y, Owa C, Kobayashi H, et al. (2026). Accessing medically relevant complex regions with a pangenome graph of 20 near-complete Japanese haplotypes. *Nature Communications*. https://doi.org/10.1038/s41467-026-73461-x


## Software and versions

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
