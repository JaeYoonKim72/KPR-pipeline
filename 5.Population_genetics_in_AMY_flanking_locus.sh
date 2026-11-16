#!/usr/bin/env bash
# AMY flanks (GRCh38), Bolognini et al. 2024 bundle 0 / bundle 1a.
# One population -> within; multiple populations -> within + first-vs-rest FST/XP-nSL.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VCF="${VCF:-}"
OUT_DIR="${OUT_DIR:-./results/amy_flanking_population_genetics}"
DISTAL="${DISTAL:-chr1:103456163-103561526}"
PROXIMAL="${PROXIMAL:-chr1:103760698-103826698}"
CALLABLE_BED="${CALLABLE_BED:-}"
PROJECTION_N="${PROJECTION_N:-168}"
CALL_RATE="${CALL_RATE:-0.90}"
THREADS="${THREADS:-4}"
POPS="${POPS:-}"
# Fast/optional and publication-scale uncertainty analyses.
# Defaults are safe exploratory settings; see the docs for manuscript settings.
INDIVIDUAL_CI="${INDIVIDUAL_CI:-0}"
IND_BOOT_REPS="${IND_BOOT_REPS:-100}"
IND_JACKKNIFE_FRACTION="${IND_JACKKNIFE_FRACTION:-0.30}"
GENOMIC_BOOTSTRAP="${GENOMIC_BOOTSTRAP:-0}"
GENOMIC_CHROMS="${GENOMIC_CHROMS:-chr1}"
GENOMIC_BLOCK_MB="${GENOMIC_BLOCK_MB:-50}"
GENOMIC_WINDOW_STEP="${GENOMIC_WINDOW_STEP:-1000}"
GENOMIC_BOOT_REPS="${GENOMIC_BOOT_REPS:-100}"
GENOMIC_FST_CALLABLE_MATCH="${GENOMIC_FST_CALLABLE_MATCH:-0.10}"
GENOME_VCF="${GENOME_VCF:-$VCF}"
GENOME_CHROM_SIZES="${GENOME_CHROM_SIZES:-}"
SEED="${SEED:-20261117}"
LD_MAF="${LD_MAF:-0.05}"
LD_DISTAL_MAX="${LD_DISTAL_MAX:-100000}"
LD_PROXIMAL_MAX="${LD_PROXIMAL_MAX:-60000}"
NSL_CORE_MAF="${NSL_CORE_MAF:-0.05}"
NSL_MAX_EXTEND="${NSL_MAX_EXTEND:-50}"
STEP="${1:-help}"; shift || true
usage() {
cat <<'EOF'
Usage: bash 5.Population_genetics_in_AMY_flanking_locus.sh all [options]
Steps: all, diversity, ld, across_ld, fst, nsl, xpnsl, individual_ci, genomic
Required: VCF=/path/to/joint.phased.biallelic.snps.vcf.gz
          POPS="KOR=/path/KOR.samples CHN=/path/CHN.samples"
          OR repeated --pop NAME=sample_file
Options: --vcf FILE --pop NAME=FILE --out DIR --mask BED
         --distal chr1:START-END --proximal chr1:START-END
         --projection 168 --call-rate 0.90
         --individual-ci --bootstrap-reps 100 --jackknife-fraction 0.30
         --genomic --chromosomes chr1 --block-mb 50 --genomic-reps 100
         --window-step 1000 --genome-vcf FILE --chrom-sizes FILE
         --fst-callable-match 0.10 --ld-maf 0.05
         --ld-distal-max 100000 --ld-proximal-max 60000
         --core-maf 0.05 --max-extend-nsl 50 --seed 20261117
A single --pop computes within only, two or more also computes comparisons
against the first population. Default coordinates are GRCh38 Bolognini 2024.
Requires Python3 and bcftools. selscan=3.0.1 is optional for nSL/XP-nSL.
Uncertainty computations default OFF. If enabled, individual CI uses 100
bootstrap reps and approximate 30%-subset delete-one JK; genomic CI uses
chr1, physical 50-Mb blocks, 1-kb sliding windows and 100 block reps.
For manuscript settings: IND_BOOT_REPS=10000 IND_JACKKNIFE_FRACTION=1
GENOMIC_CHROMS=chr1,...,chr22 GENOMIC_BLOCK_MB=5 GENOMIC_BOOT_REPS=10000.
Genomic percentile currently covers pi, thetaW, Tajima D and Hudson FST;
other statistics require independently generated normalized genomewide windows.
EOF
}
case "$STEP" in
  help|-h|--help) usage; exit 0 ;;
  all|diversity|ld|across_ld|fst|nsl|xpnsl|individual_ci|genomic) ;;
  *) usage; exit 2 ;;
esac
POPARGS=()
while (( $# )); do
 case "$1" in
 --vcf) VCF="$2"; shift 2;;
 --pop) POPARGS+=("$2"); shift 2;;
 --out) OUT_DIR="$2"; shift 2;;
 --mask) CALLABLE_BED="$2"; shift 2;;
 --distal) DISTAL="$2"; shift 2;;
 --proximal) PROXIMAL="$2"; shift 2;;
 --projection) PROJECTION_N="$2"; shift 2;;
 --call-rate) CALL_RATE="$2"; shift 2;;
 --individual-ci) INDIVIDUAL_CI=1; shift;;
 --genomic) GENOMIC_BOOTSTRAP=1; shift;;
 --bootstrap-reps) IND_BOOT_REPS="$2"; shift 2;;
 --jackknife-fraction) IND_JACKKNIFE_FRACTION="$2"; shift 2;;
 --chromosomes) GENOMIC_CHROMS="$2"; shift 2;;
 --block-mb) GENOMIC_BLOCK_MB="$2"; shift 2;;
 --genomic-reps) GENOMIC_BOOT_REPS="$2"; shift 2;;
 --window-step) GENOMIC_WINDOW_STEP="$2"; shift 2;;
 --genome-vcf) GENOME_VCF="$2"; shift 2;;
 --chrom-sizes) GENOME_CHROM_SIZES="$2"; shift 2;;
 --fst-callable-match) GENOMIC_FST_CALLABLE_MATCH="$2"; shift 2;;
 --ld-maf) LD_MAF="$2"; shift 2;;
 --ld-distal-max) LD_DISTAL_MAX="$2"; shift 2;;
 --ld-proximal-max) LD_PROXIMAL_MAX="$2"; shift 2;;
 --core-maf) NSL_CORE_MAF="$2"; shift 2;;
 --max-extend-nsl) NSL_MAX_EXTEND="$2"; shift 2;;
 --seed) SEED="$2"; shift 2;;
 *) echo "[ERROR] unknown argument: $1" >&2; exit 2;;
 esac
done
# --vcf is parsed after the environment defaults; inherit it for genomic runs.
[[ -n "$GENOME_VCF" ]] || GENOME_VCF="$VCF"
if [[ -n "$POPS" ]]; then
 read -r -a FROMENV <<< "$POPS"
 POPARGS=("${FROMENV[@]}" "${POPARGS[@]}")
fi
[[ -s "$VCF" ]] || { echo "[ERROR] VCF is required" >&2; exit 2; }
(( ${#POPARGS[@]} )) || { echo "[ERROR] supply --pop or POPS" >&2; exit 2; }
for binary in python3 bcftools; do
 command -v "$binary" >/dev/null || { echo "[ERROR] missing $binary" >&2; exit 2; }
done
mkdir -p "$OUT_DIR"/{stats,selscan,intermediate,logs,individual_ci,genomic}
printf 'step\tpopulation\tregion\tstatus\tnote\n' > "$OUT_DIR/status.tsv"
status() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" >> "$OUT_DIR/status.tsv"; }
NAMES=(); LISTS=(); PYPOPS=()
for x in "${POPARGS[@]}"; do
 [[ "$x" == *=* ]] || { echo "[ERROR] expected NAME=FILE: $x" >&2; exit 2; }
 name="${x%%=*}"; list="${x#*=}"
 [[ "$name" =~ ^[A-Za-z0-9_-]+$ && -s "$list" ]] || { echo "[ERROR] invalid population: $x" >&2; exit 2; }
 NAMES+=("$name"); LISTS+=("$list"); PYPOPS+=(--pop "$x")
done
MASKOPT=()
if [[ -n "$CALLABLE_BED" ]]; then
 [[ -s "$CALLABLE_BED" ]] || { echo "[ERROR] mask not found" >&2; exit 2; }
 MASKOPT=(--callable-bed "$CALLABLE_BED")
else
 echo "[WARNING] No StrictMask BED: physical length, not callable bp, used for diversity denominator." >&2
fi
# Use indexed region extraction to avoid streaming whole-autosome VCFs for
# focal estimators and individual resampling. Never reuse stale extracted VCFs.
prepare_focal_vcf() {
 if [[ -z "${FOCAL_VCF:-}" ]]; then
   FOCAL_VCF="$OUT_DIR/intermediate/focal.vcf.gz"
   bcftools view -r "$DISTAL,$PROXIMAL" -m2 -M2 -v snps \
     -Oz -o "$FOCAL_VCF" "$VCF"
   bcftools index -f "$FOCAL_VCF"
 fi
}
calc() {
 prepare_focal_vcf
 python3 "$SCRIPT_DIR/amy_flank_popgen.py" --vcf "$FOCAL_VCF" --out "$OUT_DIR/stats" \
 --distal "$DISTAL" --proximal "$PROXIMAL" --projection "$PROJECTION_N" \
 --callrate "$CALL_RATE" --ld-maf "$LD_MAF" --distal-ld-max "$LD_DISTAL_MAX" \
 --proximal-ld-max "$LD_PROXIMAL_MAX" --only "$1" "${MASKOPT[@]}" "${PYPOPS[@]}"
 status "$1" all both DONE focal_point_estimate
}
SELSCAN_OK=0
check_selscan() {
 local detected=""
 if command -v selscan >/dev/null 2>&1; then
   detected="$(selscan --version 2>&1 || true)"
   [[ "$detected" == *3.0.1* ]] && SELSCAN_OK=1
 fi
 if (( ! SELSCAN_OK )); then
  echo "[WARNING] selscan v3.0.1 not found (detected: ${detected:-none})." >&2
  echo "[WARNING] Install with mamba/conda in an activated environment:" >&2
  echo "  mamba create -n amy_popgen -c conda-forge -c bioconda selscan=3.0.1" >&2
  echo "  conda activate amy_popgen" >&2
  echo "[WARNING] nSL/XP-nSL skipped; other statistics continue." >&2
 fi
}
EXTRACTED=0
extract_inputs() {
 (( EXTRACTED )) && return 0
 for reg in distal proximal; do
  if [[ "$reg" == distal ]]; then coords="$DISTAL"; else coords="$PROXIMAL"; fi
  for i in "${!NAMES[@]}"; do
   name="${NAMES[$i]}"
   dest="$OUT_DIR/intermediate/$name.$reg.vcf.gz"
   if [[ ! -s "$dest" ]]; then
    if [[ -n "$CALLABLE_BED" ]]; then
      bcftools view -r "$coords" -T "$CALLABLE_BED" -S "${LISTS[$i]}" -m2 -M2 -v snps \
        -Oz -o "$dest" "$VCF"
    else
      bcftools view -r "$coords" -S "${LISTS[$i]}" -m2 -M2 -v snps \
        -Oz -o "$dest" "$VCF"
    fi
    bcftools index -f "$dest"
   fi
  done
 done
 EXTRACTED=1
}
run_nsl() {
 if (( ! SELSCAN_OK )); then
   for name in "${NAMES[@]}"; do
    for reg in distal proximal; do status nsl "$name" "$reg" SKIPPED_MISSING_SELSCAN need_v3.0.1; done
   done
   return
 fi
 extract_inputs
 for name in "${NAMES[@]}"; do
  for reg in distal proximal; do
   prefix="$OUT_DIR/selscan/$name.$reg.nsl"
   if selscan --nsl --unphased --vcf "$OUT_DIR/intermediate/$name.$reg.vcf.gz" \
     --maf "$NSL_CORE_MAF" --keep-low-freq --max-extend-nsl "$NSL_MAX_EXTEND" --threads "$THREADS" \
     --out "$prefix" >"$OUT_DIR/logs/$name.$reg.nsl.log" 2>&1; then
     status nsl "$name" "$reg" RAW_ONLY needs_autosomal_AF_bin_normalization
   else
     status nsl "$name" "$reg" FAILED see_log
   fi
  done
 done
}
run_xpnsl() {
 if (( ${#NAMES[@]} < 2 )); then
   status xpnsl NA both SKIPPED_ONE_POP requires_two_populations
   return
 fi
 if (( ! SELSCAN_OK )); then
   for ((i=1;i<${#NAMES[@]};i++)); do
    for reg in distal proximal; do
      status xpnsl "${NAMES[0]}_vs_${NAMES[$i]}" "$reg" SKIPPED_MISSING_SELSCAN need_v3.0.1
    done
   done
   return
 fi
 extract_inputs
 for ((i=1;i<${#NAMES[@]};i++)); do
  pair="${NAMES[0]}_vs_${NAMES[$i]}"
  for reg in distal proximal; do
   work="$OUT_DIR/intermediate/$pair.$reg"
   python3 "$SCRIPT_DIR/amy_flank_population_stats.py" merge \
     --query "$OUT_DIR/intermediate/${NAMES[0]}.$reg.vcf.gz" \
     --reference-vcf "$OUT_DIR/intermediate/${NAMES[$i]}.$reg.vcf.gz" \
     --query-label "${NAMES[0]}" --reference-label "${NAMES[$i]}" \
     --output "$work.merged.vcf" --query-list "$work.query.ids" --reference-list "$work.ref.ids" \
     --query-matched "$work.query.vcf" --reference-matched "$work.ref.vcf" \
     >"$OUT_DIR/logs/$pair.$reg.match.log"
   prefix="$OUT_DIR/selscan/$pair.$reg.xpnsl"
   if selscan --xpnsl --vcf "$work.query.vcf" --vcf-ref "$work.ref.vcf" \
      --maf "$NSL_CORE_MAF" --keep-low-freq --max-extend-nsl "$NSL_MAX_EXTEND" --threads "$THREADS" \
      --out "$prefix" >"$OUT_DIR/logs/$pair.$reg.xpnsl.log" 2>&1; then
     status xpnsl "$pair" "$reg" RAW_ONLY needs_genomewide_normalization
   else
     status xpnsl "$pair" "$reg" FAILED see_log
   fi
  done
 done
}
# Optional CI for SNP-level pi/thetaW/TajimaD/within-LD/across-LD/FST.
# selscan-based nSL/XP-nSL CI needs separate reruns + genomewide normalization.
run_individual_ci() {
 prepare_focal_vcf
 python3 "$SCRIPT_DIR/amy_flank_individual_ci.py" \
  --vcf "$FOCAL_VCF" --out "$OUT_DIR/individual_ci" \
  --distal "$DISTAL" --proximal "$PROXIMAL" \
  --projection "$PROJECTION_N" --callrate "$CALL_RATE" \
  --bootstrap-reps "$IND_BOOT_REPS" \
  --jackknife-fraction "$IND_JACKKNIFE_FRACTION" \
  --ld-maf "$LD_MAF" --distal-ld-max "$LD_DISTAL_MAX" \
  --proximal-ld-max "$LD_PROXIMAL_MAX" --seed "$SEED" \
  "${MASKOPT[@]}" "${PYPOPS[@]}"
 status individual_ci all both DONE "SNP_statistics_only;nSL_XPnSL_not_resampled"
}

run_genomic() {
 [[ -s "$GENOME_VCF" ]] || { echo "[ERROR] GENOME_VCF missing: $GENOME_VCF" >&2; return 1; }
 # A stand-alone genomic run creates focal point values if not yet present.
 if [[ ! -s "$OUT_DIR/stats/diversity.tsv" ]]; then calc diversity; fi
 if (( ${#NAMES[@]} > 1 )) && [[ ! -s "$OUT_DIR/stats/hudson_fst.tsv" ]]; then calc fst; fi
 local additional=()
 [[ -n "$GENOME_CHROM_SIZES" ]] && additional=(--chrom-sizes "$GENOME_CHROM_SIZES")
 python3 "$SCRIPT_DIR/amy_flank_genomic_background.py" \
  --vcf "$GENOME_VCF" --focal-stats "$OUT_DIR/stats" \
  --out "$OUT_DIR/genomic" --distal "$DISTAL" --proximal "$PROXIMAL" \
  --chromosomes "$GENOMIC_CHROMS" --block-mb "$GENOMIC_BLOCK_MB" \
  --window-step "$GENOMIC_WINDOW_STEP" --bootstrap-reps "$GENOMIC_BOOT_REPS" \
  --fst-callable-match "$GENOMIC_FST_CALLABLE_MATCH" \
  --projection "$PROJECTION_N" --callrate "$CALL_RATE" --seed "$SEED" \
  "${MASKOPT[@]}" "${additional[@]}" "${PYPOPS[@]}"
 status genomic all both DONE "pi_thetaW_TajimaD_HudsonFST_only"
}
case "$STEP" in
 all)
   calc all
   check_selscan
   run_nsl
   run_xpnsl
   if [[ "$INDIVIDUAL_CI" == 1 ]]; then run_individual_ci; fi
   if [[ "$GENOMIC_BOOTSTRAP" == 1 ]]; then run_genomic; fi
   ;;
 diversity|ld|across_ld|fst) calc "$STEP";;
 nsl) check_selscan; run_nsl;;
 xpnsl) check_selscan; run_xpnsl;;
 individual_ci) run_individual_ci;;
 genomic) run_genomic;;
esac
echo "[INFO] Finished: $OUT_DIR/status.tsv"
