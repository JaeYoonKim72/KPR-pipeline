#!/usr/bin/env bash
#
# Constructing a Korean pangenome reveals distinct haplotypes at the amylase
# locus
#
# pangenome_variant_benchmarking.sh
#
# Usage: bash pangenome_variant_benchmarking.sh <step>
#        (see docs for the step list and for the required inputs)

set -euo pipefail

###############################################################################
# Configuration
###############################################################################

JOBS=${JOBS:-8}
THREADS=${THREADS:-4}

REF_DIR=${REF_DIR:?set REF_DIR to the directory containing the reference files}
TRUTH_DIR=${TRUTH_DIR:?set TRUTH_DIR to the directory containing the truth-set inputs}
QUERY_DIR=${QUERY_DIR:?set QUERY_DIR to the directory containing the call sets to evaluate}
OUT_DIR=${OUT_DIR:-$(pwd)/results}
LIST_DIR=${LIST_DIR:-$(pwd)/lists}
TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}
SCRIPT_DIR=${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

GRCH38=${GRCH38:-${REF_DIR}/GRCh38.fa}
GRCH38_SDF=${GRCH38_SDF:-${REF_DIR}/GRCh38.sdf}

SAMPLE_LIST=${SAMPLE_LIST:-${LIST_DIR}/benchmark_samples.list}
HIFI_DV_DIR=${HIFI_DV_DIR:-${TRUTH_DIR}/hifi_deepvariant}
LR_SV_DIR=${LR_SV_DIR:-${TRUTH_DIR}/long_read_sv}

GRAPHS=${GRAPHS:-"kpr hprc_kpr hprc"}
LINEAR_SMALL=${LINEAR_SMALL:-gatk}
SV_SUPPORT_LEVELS=${SV_SUPPORT_LEVELS:-"ge1 ge2 ge3"}

TRUTH_MIN_GQ=${TRUTH_MIN_GQ:-20}
TRUTH_MIN_DP=${TRUTH_MIN_DP:-7}
HET_VAF_MIN=${HET_VAF_MIN:-0.30}
HET_VAF_MAX=${HET_VAF_MAX:-0.70}
HOM_VAF_MIN=${HOM_VAF_MIN:-0.90}

SUPPORT_REFDIST=${SUPPORT_REFDIST:-500}
SUPPORT_PCTSIZE=${SUPPORT_PCTSIZE:-0.7}
SUPPORT_PCTSEQ=${SUPPORT_PCTSEQ:-0.7}
READSUP_MIN_ALT=${READSUP_MIN_ALT:-5}
READSUP_MIN_DP=${READSUP_MIN_DP:-10}

TV_REFDIST=${TV_REFDIST:-500}
TV_PCTSIZE=${TV_PCTSIZE:-0.8}
TV_PCTSEQ=${TV_PCTSEQ:-0.8}
TV_PCTOVL=${TV_PCTOVL:-0.8}
TV_SIZEMIN=${TV_SIZEMIN:-51}
TV_SIZEFILT=${TV_SIZEFILT:-30}
TV_SIZEMAX=${TV_SIZEMAX:-50000}
SV_MIN_LEN=${SV_MIN_LEN:-51}

SMALL_GT_EXPR='((GT="het" || GT="AA") && GT!~"\\.")'
SV_GT_EXPR='GT="alt" && GT!~"\\."'

BENCH_OUT=${OUT_DIR}/benchmarking
TRUTH_OUT=${BENCH_OUT}/truth
QUERY_OUT=${BENCH_OUT}/query
RESULT_OUT=${BENCH_OUT}/bench
BED=${BENCH_OUT}/autosomes.bed
LOG_OUT=${OUT_DIR}/logs
JOB_LOG=${LOG_OUT}/benchmarking

mkdir -p "${TRUTH_OUT}" "${QUERY_OUT}" "${RESULT_OUT}" "${LOG_OUT}" "${JOB_LOG}" \
         "${TMP_DIR}" "${LIST_DIR}"

export HGREF=${GRCH38}

log() { echo -e "[$(date '+%F %T')] $*" | tee -a "${LOG_OUT}/pipeline.log"; }

samples() { grep -v '^[[:space:]]*$' "${SAMPLE_LIST}"; }

index_vcf() { bcftools index -f -t "$1"; }

run_jobs() {
  local FN=$1; shift
  local PIDS=() X P FAIL=0
  for X in "$@"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do sleep 5; done
    ( ${FN} ${X} ) > "${JOB_LOG}/${FN}.${X// /.}.log" 2>&1 &
    PIDS+=($!)
  done
  for P in "${PIDS[@]}"; do wait "${P}" || FAIL=1; done
  (( FAIL == 0 )) || { log "  ${FN} failed, see ${JOB_LOG}/${FN}.*.log"; exit 1; }
}

autosome_bed() {
  [[ -s "${GRCH38}.fai" ]] || samtools faidx "${GRCH38}"
  [[ -s "${BED}" ]] || awk -v OFS='\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $1, 0, $2}' \
    "${GRCH38}.fai" > "${BED}"
  [[ $(wc -l < "${BED}") -eq 22 ]] || { log "  ${BED} must hold chr1-chr22"; exit 1; }
}

###############################################################################
# Step 1. Small variant truth sets
###############################################################################

small_truth_sample() {
  local S=$1
  local OUT="${TRUTH_OUT}/small/${S}"
  local SNP_IN="${HIFI_DV_DIR}/${S}.hifi_dv.biSNP.addfilter.VAF0208.vcf.gz"
  local SMALL_IN="${HIFI_DV_DIR}/${S}.hifi_dv.Small.vcf.gz"
  local GQ="FMT/GQ>=${TRUTH_MIN_GQ} && FMT/DP>=${TRUTH_MIN_DP}"
  local VAF="(GT=\"het\" && FMT/VAF>=${HET_VAF_MIN} && FMT/VAF<=${HET_VAF_MAX}) || (GT=\"AA\" && FMT/VAF>=${HOM_VAF_MIN})"
  mkdir -p "${OUT}"

  # 1-1) SNPs
  bcftools view -R "${BED}" "${SNP_IN}" -Ou \
    | bcftools norm -f "${GRCH38}" -m -any -Ou \
    | bcftools view -m2 -M2 -v snps -f PASS -Ou \
    | bcftools +fill-tags -Ou -- -t FORMAT/VAF \
    | bcftools view -i "${GQ} && (${VAF})" -Oz -o "${OUT}/${S}.SNP.truth.vcf.gz"
  index_vcf "${OUT}/${S}.SNP.truth.vcf.gz"

  # 1-2) small deletions and insertions
  bcftools view -R "${BED}" -v indels -i "${GQ} && (${VAF}) && strlen(REF)>strlen(ALT)" \
    "${SMALL_IN}" -Oz -o "${OUT}/${S}.DEL.truth.vcf.gz"
  bcftools view -R "${BED}" -v indels -i "${GQ} && (${VAF}) && strlen(ALT)>strlen(REF)" \
    "${SMALL_IN}" -Oz -o "${OUT}/${S}.INS.truth.vcf.gz"
  index_vcf "${OUT}/${S}.DEL.truth.vcf.gz"
  index_vcf "${OUT}/${S}.INS.truth.vcf.gz"
}

small_truth() {
  log "step 1: small variant truth sets"
  autosome_bed
  run_jobs small_truth_sample $(samples)
  log "step 1 finished"
}

###############################################################################
# Step 2. Long-read SV call filtering
###############################################################################

sv_filter_caller() {
  local S=$1 C=$2
  local OUT="${TRUTH_OUT}/sv/${S}"
  mkdir -p "${OUT}"

  # 2-1) PASS INS/DEL, complete non-reference genotype, read support
  python3 "${SCRIPT_DIR}/sv_truth_helper.py" filter \
    --input "${LR_SV_DIR}/${C}_${S}.PASS.SVLEN51.chr.vcf.gz" \
    --output "${OUT}/${S}.${C}.alt.vcf.gz" \
    --caller "${C}" --sample "${S}"
  index_vcf "${OUT}/${S}.${C}.alt.vcf.gz"
}

sv_truth_filter() {
  log "step 2: long-read SV call filtering"
  local ITEMS=() S C
  for S in $(samples); do
    for C in svim_asm pbsv sniffles svim; do ITEMS+=("${S} ${C}"); done
  done
  run_jobs sv_filter_caller "${ITEMS[@]}"
  log "step 2 finished"
}

###############################################################################
# Step 3. Assembly SVs supported by read-based callers
###############################################################################

sv_support_caller() {
  local S=$1 C=$2
  local D="${TRUTH_OUT}/sv/${S}"
  rm -rf "${D}/bench_${C}"

  # 3-1) SVIM-asm SVs matched by each read-based caller
  truvari bench -b "${D}/${S}.svim_asm.alt.vcf.gz" -c "${D}/${S}.${C}.alt.vcf.gz" \
    -f "${GRCH38}" -o "${D}/bench_${C}" \
    --refdist "${SUPPORT_REFDIST}" --pctsize "${SUPPORT_PCTSIZE}" --pctseq "${SUPPORT_PCTSEQ}" \
    --sizemin 50
}

sv_support_levels() {
  local S=$1
  local D="${TRUTH_OUT}/sv/${S}"

  # 3-2) number of supporting callers, truth sets ge1, ge2 and ge3
  python3 "${SCRIPT_DIR}/sv_truth_helper.py" collect --sample "${S}" --dir "${D}"
  python3 "${SCRIPT_DIR}/sv_truth_helper.py" levels --sample "${S}" --dir "${D}"
}

sv_truth_support() {
  log "step 3: assembly SVs supported by read-based callers"
  local ITEMS=() S C
  for S in $(samples); do
    for C in pbsv sniffles svim; do ITEMS+=("${S} ${C}"); done
  done
  run_jobs sv_support_caller "${ITEMS[@]}"
  run_jobs sv_support_levels $(samples)
  log "step 3 finished"
}

###############################################################################
# Step 4. Read-support filter of the SV truth sets
###############################################################################

prepare_sv() {
  local IN=$1 OUT=$2 TYPE=$3
  local T
  T=$(mktemp -d "${TMP_DIR}/prep.XXXXXX")
  python3 "${SCRIPT_DIR}/sv_truth_helper.py" fixheader --input "${IN}" --output "${T}/h.vcf"
  bcftools norm -m -any "${T}/h.vcf" -Ov -o "${T}/split.vcf"
  bcftools sort -T "${T}/sort" "${T}/split.vcf" -Ov -o "${T}/sorted.vcf"
  python3 "${SCRIPT_DIR}/prepare_sv_autosome.py" --input "${T}/sorted.vcf" \
    --output "${T}/typed.vcf" --type "${TYPE}" --min-size "${SV_MIN_LEN}"
  bcftools view "${T}/typed.vcf" -Oz -o "${OUT}"
  index_vcf "${OUT}"
  rm -rf "${T}"
}

sv_readsupport_level() {
  local S=$1 L=$2
  local D="${TRUTH_OUT}/sv/${S}"
  local OUT="${D}/${L}"
  local C T F
  mkdir -p "${OUT}"
  for C in sniffles pbsv svim; do
    F="${LR_SV_DIR}/${C}_${S}.PASS.SVLEN51.chr.vcf.gz"
    [[ -s "${F}.tbi" || -s "${F}.csi" ]] || index_vcf "${F}"
  done

  # 4-1) ALT reads >= READSUP_MIN_ALT and depth >= READSUP_MIN_DP in the same caller
  python3 "${SCRIPT_DIR}/sv_truth_readsupport.py" \
    --truth "${D}/${S}.SV.truth.${L}.vcf.gz" \
    --sniffles "${LR_SV_DIR}/sniffles_${S}.PASS.SVLEN51.chr.vcf.gz" \
    --pbsv "${LR_SV_DIR}/pbsv_${S}.PASS.SVLEN51.chr.vcf.gz" \
    --svim "${LR_SV_DIR}/svim_${S}.PASS.SVLEN51.chr.vcf.gz" \
    --min-alt "${READSUP_MIN_ALT}" --min-dp "${READSUP_MIN_DP}" \
    --outdir "${OUT}"

  # 4-2) autosomal deletions and insertions
  for T in DEL INS; do
    prepare_sv "${OUT}/truth.readsupport.vcf.gz" "${OUT}/truth.${T}.vcf.gz" "${T}"
  done
}

sv_truth_readsupport() {
  log "step 4: read-support filter of the SV truth sets"
  local ITEMS=() S L
  for S in $(samples); do
    for L in ${SV_SUPPORT_LEVELS}; do ITEMS+=("${S} ${L}"); done
  done
  run_jobs sv_readsupport_level "${ITEMS[@]}"
  log "step 4 finished"
}

###############################################################################
# Step 5. Small variant queries
###############################################################################

small_query_set() {
  local S=$1 G=$2
  local IN="${QUERY_DIR}/${G}/${G}.${S}"
  local OUT="${QUERY_OUT}/small/${S}"
  local T X
  T=$(mktemp -d "${TMP_DIR}/sq.XXXXXX")
  mkdir -p "${OUT}"

  # 5-1) SNPs
  bcftools view -T "${BED}" -f PASS,. "${IN}.snp.vcf.gz" -Ou \
    | bcftools norm -f "${GRCH38}" -m -any -Ou \
    | bcftools view -m2 -M2 -v snps -i "${SMALL_GT_EXPR}" -Ou \
    | bcftools sort -T "${T}" -Oz -o "${OUT}/${G}.SNP.vcf.gz"

  # 5-2) small deletions and insertions
  bcftools view -T "${BED}" -f PASS,. "${IN}.indel.vcf.gz" -Ou \
    | bcftools norm -f "${GRCH38}" -m -any -Ou \
    | bcftools sort -T "${T}" -Ob -o "${T}/indel.bcf"
  bcftools view -m2 -M2 -v indels -i "(STRLEN(REF)>STRLEN(ALT)) && ${SMALL_GT_EXPR}" \
    "${T}/indel.bcf" -Oz -o "${OUT}/${G}.DEL.vcf.gz"
  bcftools view -m2 -M2 -v indels -i "(STRLEN(REF)<STRLEN(ALT)) && ${SMALL_GT_EXPR}" \
    "${T}/indel.bcf" -Oz -o "${OUT}/${G}.INS.vcf.gz"

  for X in SNP DEL INS; do index_vcf "${OUT}/${G}.${X}.vcf.gz"; done
  rm -rf "${T}"
}

small_query() {
  log "step 5: small variant queries"
  autosome_bed
  local ITEMS=() S G
  for S in $(samples); do
    for G in ${GRAPHS} ${LINEAR_SMALL}; do ITEMS+=("${S} ${G}"); done
  done
  run_jobs small_query_set "${ITEMS[@]}"
  log "step 5 finished"
}

###############################################################################
# Step 6. Small variant benchmark
###############################################################################

small_bench_set() {
  local S=$1 G=$2 TYPE=$3
  local OUT="${RESULT_OUT}/small/${S}"
  local REGION=(-f "${BED}")
  [[ ${TYPE} == SNP ]] && REGION=(-R "${BED}")
  mkdir -p "${OUT}"

  # 6-1) hap.py with the vcfeval engine
  hap.py "${TRUTH_OUT}/small/${S}/${S}.${TYPE}.truth.vcf.gz" \
    "${QUERY_OUT}/small/${S}/${G}.${TYPE}.vcf.gz" \
    -r "${GRCH38}" "${REGION[@]}" -o "${OUT}/${G}.${TYPE}" --threads "${THREADS}" \
    --engine vcfeval --engine-vcfeval-template "${GRCH38_SDF}" --no-roc
}

small_bench() {
  log "step 6: small variant benchmark"
  local ITEMS=() S G T
  for S in $(samples); do
    for G in ${GRAPHS} ${LINEAR_SMALL}; do
      for T in SNP DEL INS; do ITEMS+=("${S} ${G} ${T}"); done
    done
  done
  run_jobs small_bench_set "${ITEMS[@]}"
  log "step 6 finished"
}

###############################################################################
# Step 7. SV queries
###############################################################################

sv_query_file() {
  local OUT=$1 TYPE_EXPR=$2; shift 2
  local PARTS=() I=0 F P
  echo SAMPLE > "${OUT}.name"
  for F in "$@"; do
    P="${OUT}.part${I}.vcf.gz"; I=$((I + 1))
    bcftools view -f PASS "${F}" -Ou \
      | bcftools norm -m -any -Ou \
      | bcftools view -i "(${SV_GT_EXPR}) && (${TYPE_EXPR})" -Ou \
      | bcftools reheader -s "${OUT}.name" - \
      | bcftools sort -Oz -o "${P}"
    index_vcf "${P}"
    PARTS+=("${P}")
  done
  if (( ${#PARTS[@]} == 1 )); then
    mv "${PARTS[0]}" "${OUT}"
  else
    bcftools concat -a "${PARTS[@]}" -Ou | bcftools sort -Oz -o "${OUT}"
  fi
  index_vcf "${OUT}"
  rm -f "${OUT}".part* "${OUT}.name"
}

sv_query_sample() {
  local S=$1
  local OUT="${QUERY_OUT}/sv/${S}"
  local G T
  mkdir -p "${OUT}"

  # 7-1) graph-based SVs
  for G in ${GRAPHS}; do
    for T in DEL INS; do
      sv_query_file "${OUT}/${G}.${T}.vcf.gz" 'POS>0' "${QUERY_DIR}/${G}/${G}.${S}.sv.${T}.vcf.gz"
    done
  done

  # 7-2) Manta deletions, and insertions with duplications
  sv_query_file "${OUT}/manta.DEL.vcf.gz" 'INFO/SVTYPE="DEL"' \
    "${QUERY_DIR}/manta/manta.${S}.sv.DEL.vcf.gz"
  sv_query_file "${OUT}/manta.INS.vcf.gz" 'INFO/SVTYPE="INS" || INFO/SVTYPE="DUP"' \
    "${QUERY_DIR}/manta/manta.${S}.sv.INS.vcf.gz" "${QUERY_DIR}/manta/manta.${S}.sv.DUP.vcf.gz"

  # 7-3) Delly deletions
  sv_query_file "${OUT}/delly.DEL.vcf.gz" 'INFO/SVTYPE="DEL"' \
    "${QUERY_DIR}/delly/delly.${S}.sv.DEL.vcf.gz"
}

sv_query() {
  log "step 7: SV queries"
  run_jobs sv_query_sample $(samples)
  log "step 7 finished"
}

###############################################################################
# Step 8. SV benchmark
###############################################################################

sv_bench_set() {
  local S=$1 G=$2 TYPE=$3 L=$4
  local OUT="${RESULT_OUT}/sv/${S}/${G}/${L}/${TYPE}"
  local EXTRA=()
  [[ ${TYPE} == INS ]] && EXTRA=(--typeignore)
  mkdir -p "$(dirname "${OUT}")"
  rm -rf "${OUT}"

  # 8-1) Truvari bench against the read-supported truth set
  truvari bench -b "${TRUTH_OUT}/sv/${S}/${L}/truth.${TYPE}.vcf.gz" \
    -c "${QUERY_OUT}/sv/${S}/${G}.${TYPE}.vcf.gz" \
    -o "${OUT}" --includebed "${BED}" \
    --refdist "${TV_REFDIST}" --pctsize "${TV_PCTSIZE}" --pctseq "${TV_PCTSEQ}" --pctovl "${TV_PCTOVL}" \
    --sizemin "${TV_SIZEMIN}" --sizefilt "${TV_SIZEFILT}" --sizemax "${TV_SIZEMAX}" \
    --passonly --pick single "${EXTRA[@]}"
}

sv_bench() {
  log "step 8: SV benchmark"
  autosome_bed
  local ITEMS=() S L G T
  for S in $(samples); do
    for L in ${SV_SUPPORT_LEVELS}; do
      for G in ${GRAPHS} manta; do
        for T in DEL INS; do ITEMS+=("${S} ${G} ${T} ${L}"); done
      done
      ITEMS+=("${S} delly DEL ${L}")
    done
  done
  run_jobs sv_bench_set "${ITEMS[@]}"
  log "step 8 finished"
}

###############################################################################
# Step 9. Summary tables
###############################################################################

summary() {
  log "step 9: summary tables"

  # 9-1) per-sample tables, means and paired Wilcoxon tests between the graphs
  python3 "${SCRIPT_DIR}/summarize_benchmark.py" --outdir "${BENCH_OUT}" \
    --samples $(samples) --graphs ${GRAPHS} \
    | tee -a "${LOG_OUT}/pipeline.log"

  log "step 9 finished"
}

###############################################################################
# Entry point
###############################################################################

usage() {
  cat <<EOF
Usage: bash $(basename "$0") <step>

  small_truth            HiFi DeepVariant small variant truth sets (SNP, DEL, INS)
  sv_truth_filter        filtering of the SVIM-asm, pbsv, Sniffles and SVIM calls
  sv_truth_support       SVIM-asm SVs supported by 1, 2 or 3 read-based callers
  sv_truth_readsupport   ALT reads >= ${READSUP_MIN_ALT} and depth >= ${READSUP_MIN_DP}, DEL and INS truth sets
  small_query            small variant queries (SNP, DEL, INS)
  small_bench            hap.py benchmark of the small variants
  sv_query               SV queries (DEL, INS)
  sv_bench               Truvari benchmark of the SVs
  summary                summary tables and paired Wilcoxon tests
  all                    run every step in order
EOF
}

main() {
  case "${1:-}" in
    small_truth)          small_truth ;;
    sv_truth_filter)      sv_truth_filter ;;
    sv_truth_support)     sv_truth_support ;;
    sv_truth_readsupport) sv_truth_readsupport ;;
    small_query)          small_query ;;
    small_bench)          small_bench ;;
    sv_query)             sv_query ;;
    sv_bench)             sv_bench ;;
    summary)              summary ;;
    all)
      small_truth
      sv_truth_filter
      sv_truth_support
      sv_truth_readsupport
      small_query
      small_bench
      sv_query
      sv_bench
      summary
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
