#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == help || "${1:-}" == --help || "${1:-}" == -h ]]; then
  printf '%s\n' 'Usage: MEND_REF=... MEND_SDF=... MEND_VCF_DIR=... MEND_TRIO=AJ|HC MEND_TRUTH_SMALL=... MEND_CONF_SMALL=... bash src/trio_mendelian_inconsistency.sh'
  exit 0
fi

MEND_REF=${MEND_REF:?Set MEND_REF to GRCh38 FASTA}
MEND_SDF=${MEND_SDF:?Set MEND_SDF to matching RTG SDF}
MEND_VCF_DIR=${MEND_VCF_DIR:?Set MEND_VCF_DIR to input genotype VCF directory}
MEND_SMALL_DIR=${MEND_SMALL_DIR:-${MEND_VCF_DIR}}
MEND_INS_DIR=${MEND_INS_DIR:-${MEND_VCF_DIR}}
MEND_DEL_DIR=${MEND_DEL_DIR:-${MEND_VCF_DIR}}
MEND_TRIO=${MEND_TRIO:-AJ}
MEND_ARM=${MEND_ARM:-graph}
MEND_REFERENCE=${MEND_REFERENCE:-${MEND_ARM}}
MEND_OUTDIR=${MEND_OUTDIR:-results/benchmarking/mendelian}
MEND_THREADS=${MEND_THREADS:-4}
MEND_TRUTH_SMALL=${MEND_TRUTH_SMALL:?Set MEND_TRUTH_SMALL to child small-variant truth VCF}
MEND_CONF_SMALL=${MEND_CONF_SMALL:?Set MEND_CONF_SMALL to corresponding high-confidence BED}
MEND_TRUTH_SV=${MEND_TRUTH_SV:-}
MEND_CONF_SV=${MEND_CONF_SV:-}
MEND_MIN_GQ=${MEND_MIN_GQ:-20}
MEND_MIN_DP=${MEND_MIN_DP:-10}
MEND_MIN_SV=51
MEND_MAX_SV=50000
MEND_SDF=${MEND_SDF}
MEND_ROOT="${MEND_OUTDIR}/${MEND_REFERENCE}.${MEND_TRIO}"

case "${MEND_TRIO}" in
  AJ) CHILD=HG002; FATHER=HG003; MOTHER=HG004 ;;
  HC) CHILD=HG005; FATHER=HG006; MOTHER=HG007 ;;
  *) printf 'Unknown trio: %s\n' "${MEND_TRIO}" >&2; exit 1 ;;
esac
[[ "${MEND_ARM}" == graph || "${MEND_ARM}" == linear ]] || { printf 'MEND_ARM must be graph or linear\n' >&2; exit 1; }
if [[ "${MEND_TRIO}" == HC && -n "${MEND_TRUTH_SV}" ]]; then
  printf 'HC trio has no supplied validated GIAB SV truth; omit MEND_TRUTH_SV\n' >&2
  exit 1
fi
for exe in bcftools rtg python3 hap.py truvari bgzip; do
  command -v "${exe}" >/dev/null || { printf 'Missing executable: %s\n' "${exe}" >&2; exit 1; }
done
mkdir -p "${MEND_ROOT}"
[[ -s "${MEND_REF}" && -d "${MEND_SDF}" && -s "${MEND_CONF_SMALL}" ]] || { printf 'Missing reference, SDF or small truth-confidence BED\n' >&2; exit 1; }
[[ -s "${MEND_TRUTH_SMALL}" ]] || { printf 'Missing child small truth VCF\n' >&2; exit 1; }
if [[ -n "${MEND_TRUTH_SV}" ]]; then
  [[ -s "${MEND_TRUTH_SV}" && -s "${MEND_CONF_SV}" ]] || { printf 'Missing SV truth or confidence BED\n' >&2; exit 1; }
fi

autosomes=$(seq 1 22 | sed 's/^/chr/' | paste -sd, -)
log() { printf '[%s] %s\n' "$(date +%F_%T)" "$*" >&2; }
index() { bcftools index -f -t "$1"; }

normalize() {
  local cls=$1 sample=$2 src=$3 out="${MEND_ROOT}/${cls}.${sample}.vcf.gz"
  [[ -s "${src}" ]] || { printf 'Missing VCF: %s\n' "${src}" >&2; return 1; }
  local tmp="${out}.tmp.vcf.gz"
  bcftools norm -m -any -f "${MEND_REF}" -c w "${src}" -Ou \
    | bcftools annotate --set-id '%CHROM_%POS_%REF_%ALT' -Ou \
    | bcftools view -r "${autosomes}" -f PASS -Ou \
    | bcftools norm -d exact -Oz -o "${tmp}"
  bcftools reheader -s <(printf '%s\n' "${sample}") "${tmp}" -o "${out}"
  rm -f "${tmp}"
  index "${out}"
}

filter_trio() {
  local cls=$1 raw="${MEND_ROOT}/${cls}.merged.vcf.gz" out="${MEND_ROOT}/${cls}.filtered.vcf.gz"
  local criterion
  if [[ "${cls}" == SNP || "${cls}" == Indel || "${MEND_ARM}" == graph ]]; then
    criterion="FMT/GQ<${MEND_MIN_GQ} || FMT/DP<${MEND_MIN_DP} || FMT/GQ=\".\" || FMT/DP=\".\""
  else
    criterion="FMT/GQ<${MEND_MIN_GQ} || FMT/GQ=\".\""
  fi
  bcftools filter -S . -e "${criterion}" "${raw}" -Ou \
    | bcftools view -Oz -o "${out}"
  index "${out}"
}

run_class() {
  local cls=$1 srcdir term raw ann query region truth source
  case "${cls}" in
    SNP|Indel) srcdir="${MEND_SMALL_DIR}"; region="${MEND_CONF_SMALL}"; truth="${MEND_TRUTH_SMALL}" ;;
    INS) srcdir="${MEND_INS_DIR}"; region="${MEND_CONF_SV}"; truth="${MEND_TRUTH_SV}" ;;
    DEL) srcdir="${MEND_DEL_DIR}"; region="${MEND_CONF_SV}"; truth="${MEND_TRUTH_SV}" ;;
  esac

  for sample in "${CHILD}" "${FATHER}" "${MOTHER}"; do
    normalize "${cls}" "${sample}" "${srcdir}/${sample}.vcf.gz"
  done

  raw="${MEND_ROOT}/${cls}.merged.vcf.gz"
  bcftools merge -Oz -o "${raw}" \
    "${MEND_ROOT}/${cls}.${CHILD}.vcf.gz" \
    "${MEND_ROOT}/${cls}.${FATHER}.vcf.gz" \
    "${MEND_ROOT}/${cls}.${MOTHER}.vcf.gz"
  index "${raw}"
  filter_trio "${cls}"
  local typ="${MEND_ROOT}/${cls}.typed.vcf.gz"
  case "${cls}" in
    SNP) term='TYPE="snp"' ;;
    Indel) term='TYPE="indel" && ABS(ILEN)>=1 && ABS(ILEN)<=50' ;;
    INS) term='MAX(ILEN)>=51 && MAX(ILEN)<=50000' ;;
    DEL) term='MIN(ILEN)<=-51 && MIN(ILEN)>=-50000' ;;
  esac
  if [[ "${MEND_ARM}" == linear && ( "${cls}" == INS || "${cls}" == DEL ) ]]; then
    case "${cls}" in
      INS) term='INFO/SVTYPE="INS" && ABS(INFO/SVLEN)>=51 && ABS(INFO/SVLEN)<=50000' ;;
      DEL) term='INFO/SVTYPE="DEL" && (INFO/END-POS)>=51 && (INFO/END-POS)<=50000 && INFO/PRECISE=1' ;;
    esac
  fi
  bcftools view -i "${term}" "${MEND_ROOT}/${cls}.filtered.vcf.gz" -Oz -o "${typ}"
  index "${typ}"

  cat > "${MEND_ROOT}/${cls}.ped" <<PED
#fam	ind	pat	mat	sex	phen
${MEND_TRIO}	${CHILD}	${FATHER}	${MOTHER}	1	0
${MEND_TRIO}	${FATHER}	0	0	1	0
${MEND_TRIO}	${MOTHER}	0	0	2	0
PED

  ann="${MEND_ROOT}/${cls}.annot.vcf.gz"
  rtg mendelian --input "${typ}" --output "${ann}" \
    --pedigree "${MEND_ROOT}/${cls}.ped" --template "${MEND_SDF}" \
    > "${MEND_ROOT}/${cls}.rtg_report.txt" 2>&1
  index "${ann}"
  bcftools view -i 'INFO/MCV!="."' "${ann}" -Oz -o "${MEND_ROOT}/${cls}.mcv.vcf.gz"
  index "${MEND_ROOT}/${cls}.mcv.vcf.gz"
  query="${MEND_ROOT}/${cls}.child_mcv.vcf.gz"
  bcftools view -s "${CHILD}" -i 'GT="alt"' "${MEND_ROOT}/${cls}.mcv.vcf.gz" -Oz -o "${query}"
  index "${query}"
  if [[ -z "${truth}" ]]; then
    return 0
  fi
  if [[ "${cls}" == SNP || "${cls}" == Indel ]]; then
    hap.py "${truth}" "${query}" -r "${MEND_REF}" -f "${region}" \
      --engine vcfeval --engine-vcfeval-template "${MEND_SDF}" \
      -o "${MEND_ROOT}/${cls}.truth" --threads "${MEND_THREADS}" --no-roc -V
  else
    local out="${MEND_ROOT}/${cls}.truvari"
    rm -rf "${out}"
    truvari bench -b "${truth}" -c "${query}" -o "${out}" \
      --includebed "${region}" --refdist 500 --pctsize 0.80 \
      --pctseq 0.80 --pctovl 0.80 --sizemin 51 --sizefilt 30 \
      --sizemax 50000 --passonly --pick single
  fi
}

for cls in SNP Indel INS DEL; do
  log "Processing ${MEND_TRIO} ${MEND_ARM} ${cls}"
  run_class "${cls}"
done

python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/summarize_mendelian.py" \
  --directory "${MEND_ROOT}" --reference "${MEND_REFERENCE}" --trio "${MEND_TRIO}" \
  --child "${CHILD}" --father "${FATHER}" --mother "${MOTHER}" \
  --small-bed "${MEND_CONF_SMALL}" --sv-bed "${MEND_CONF_SV}"
log "Finished: ${MEND_ROOT}/mendelian_summary.tsv"
