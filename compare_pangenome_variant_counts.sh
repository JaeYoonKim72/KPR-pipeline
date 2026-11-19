#!/usr/bin/env bash
#
# Constructing a Korean pangenome reveals distinct haplotypes at the amylase
# locus
#
# compare_pangenome_variant_counts.sh
#
# Usage: bash compare_pangenome_variant_counts.sh <step>
#        (see docs for the step list and for the required inputs)

set -euo pipefail

###############################################################################
# Configuration
###############################################################################

JOBS=${JOBS:-4}                       # chromosomes or samples processed in parallel

REF_DIR=${REF_DIR:-}                  # required by prepare_input
OUT_DIR=${OUT_DIR:-$(pwd)/results}
LIST_DIR=${LIST_DIR:-$(pwd)/lists}
TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}
SCRIPT_DIR=${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

# Graph whose variants are compared. The VCF produced by
# pangenome_construction_and_statistics.sh (step variants) is used as input:
# vg deconstruct -P GRCh38 -a, all nesting levels, no vcfbub.
GRAPH=${GRAPH:-hprc_kpr}
INPUT_VCF=${INPUT_VCF:-${OUT_DIR}/graph_variants/${GRAPH}/reference/${GRAPH}.GRCh38.vcf.gz}
GRAPH_XG=${GRAPH_XG:-${OUT_DIR}/index/${GRAPH}/${GRAPH}.xg}

# Sample groups. Two columns, sample<TAB>population; samples of TARGET_POP are
# compared with all other samples.
POP_TABLE=${POP_TABLE:-${LIST_DIR}/population.tsv}
TARGET_POP=${TARGET_POP:-KOR}
TARGET_LABEL=${TARGET_LABEL:-KPR}
OTHER_LABEL=${OTHER_LABEL:-HPRC}

# Reference haplotype stored as a sample column of the VCF, removed before the
# analysis, and the prefix of the reference contig names.
EXCLUDE_SAMPLES=${EXCLUDE_SAMPLES:-CHM13v2}
REF_SAMPLE=${REF_SAMPLE:-GRCh38}
CHROM_PREFIX=${CHROM_PREFIX:-${REF_SAMPLE}#0#}
CHROMS=${CHROMS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"}

# Bubble selection, as in APGp1: top-level bubbles, bubbles with a reference
# allele above 10 Mb replaced by their child bubbles.
VCFBUB_LEVEL=${VCFBUB_LEVEL:-0}
VCFBUB_MAX_REF=${VCFBUB_MAX_REF:-10000000}

# Alleles longer than this are not realigned by vcfwave.
WAVE_MAXLEN=${WAVE_MAXLEN:-100000}

# Size threshold: |len(ALT) - len(REF)| > SV_MIN_LEN is a structural variant,
# anything shorter a small variant.
SV_MIN_LEN=${SV_MIN_LEN:-50}

# PanSVMerger allele clustering: vsearch identity, parallel vsearch processes
# per chromosome, and the memory cap of vcfcreatemulti in GB.
MERGE_IDENTITY=${MERGE_IDENTITY:-0.9}
MERGE_WORKERS=${MERGE_WORKERS:-16}
MULTI_MEM_GB=${MULTI_MEM_GB:-32}

# Records with an allele longer than this, and inversions, are not passed to
# rtg vcfdecompose (they are structural variants; rtg splits an inversion into
# thousands of artificial SNVs and its memory grows with the allele length).
RTG_MAXLEN=${RTG_MAXLEN:-5000}
RTG_MEM=${RTG_MEM:-8g}
SORT_MEM=${SORT_MEM:-8G}

# --- reference files ---------------------------------------------------------
# UCSC hg38 centromere table and gap table (type = telomere), excluded regions
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
CENTROMERES=${CENTROMERES:-${REF_DIR}/hg38.centromeres.txt.gz}
GAPS=${GAPS:-${REF_DIR}/hg38.gap.txt.gz}

# PanSVMerger, commit 9d9827f
#   https://github.com/Asian-Pan-Genome/PanSVMerger
PANSVMERGER=${PANSVMERGER:-${SCRIPT_DIR}/PanSVMerger/scripts/SV.merge_variants.py}

CMP_OUT=${OUT_DIR}/variant_comparison
IN_OUT=${CMP_OUT}/input
SV_OUT=${CMP_OUT}/sv
SMALL_OUT=${CMP_OUT}/small_variants
COUNT_OUT=${CMP_OUT}/counts
LOG_OUT=${OUT_DIR}/logs
JOB_LOG=${LOG_OUT}/variant_comparison

setup_dirs() {
  mkdir -p "${IN_OUT}" "${SV_OUT}" "${SMALL_OUT}" "${COUNT_OUT}" "${LOG_OUT}" \
           "${JOB_LOG}" "${TMP_DIR}" "${LIST_DIR}"
}

log() { echo -e "[$(date '+%F %T')] $*" | tee -a "${LOG_OUT}/pipeline.log"; }

# run_jobs <function> <item>... : run "<function> <item>" for every item, JOBS
# at a time, one log file per item; stop when any of them fails.
run_jobs() {
  local FN=$1; shift
  local PIDS=() X P FAIL=0
  for X in "$@"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do sleep 5; done
    ( "${FN}" "${X}" ) > "${JOB_LOG}/${FN}.${X}.log" 2>&1 &
    PIDS+=($!)
  done
  for P in "${PIDS[@]}"; do wait "${P}" || FAIL=1; done
  (( FAIL == 0 )) || { log "  ${FN} failed, see ${JOB_LOG}/${FN}.*.log"; exit 1; }
}

###############################################################################
# Step 1. Input preparation
###############################################################################

prepare_input() {
  : "${REF_DIR:?set REF_DIR to the directory containing the reference files}"
  log "step 1: input preparation"

  # 1-1) centromeres and telomeres, excluded from the analysis
  { zcat "${CENTROMERES}" | awk -v OFS='\t' '{print $2, $3, $4}'
    zcat "${GAPS}" | awk -v OFS='\t' '$8 == "telomere" {print $2, $3, $4}'
  } | sort -k1,1 -k2,2n > "${IN_OUT}/exclude.cen_telo.bed"

  # 1-2) contig rename table, GRCh38#0#chr1 -> chr1
  bcftools view -h "${INPUT_VCF}" | grep '^##contig=<ID=' \
    | sed 's/^##contig=<ID=\([^,>]*\).*/\1/' \
    | awk -v p="${CHROM_PREFIX}" '{n = $1; if (index(n, p) == 1) n = substr(n, length(p) + 1); print $1 "\t" n}' \
    > "${IN_OUT}/rename_chrs.txt"

  # 1-3) common preprocessing of SVs and small variants: bubble selection with
  #      vcfbub, contig renaming, reference haplotype removed, centromeres and
  #      telomeres removed, sites without a carrier removed
  vcfbub -l "${VCFBUB_LEVEL}" -r "${VCFBUB_MAX_REF}" -i "${INPUT_VCF}" \
    | bcftools annotate --rename-chrs "${IN_OUT}/rename_chrs.txt" -Ou \
    | bcftools view -s "^${EXCLUDE_SAMPLES}" --force-samples -t "$(echo ${CHROMS} | tr ' ' ',')" -Ou \
    | bcftools view -T "^${IN_OUT}/exclude.cen_telo.bed" -Ou \
    | bcftools +fill-tags -Ou -- -t AC,AN \
    | bcftools view -c 1 -Oz -o "${IN_OUT}/all.vcf.gz"
  bcftools index -f -t "${IN_OUT}/all.vcf.gz"
  bcftools query -l "${IN_OUT}/all.vcf.gz" > "${IN_OUT}/samples.txt"

  # 1-4) SV input: sites with an allele longer than SV_MIN_LEN. Shorter sites
  #      cannot yield an SV, so this only shortens the vcfwave run.
  bcftools view -i "strlen(REF)>${SV_MIN_LEN} || strlen(ALT)>${SV_MIN_LEN}" \
    "${IN_OUT}/all.vcf.gz" -Oz -o "${IN_OUT}/sv_input.vcf.gz"
  bcftools index -f -t "${IN_OUT}/sv_input.vcf.gz"

  # 1-5) linear reference taken from the graph, used by bcftools norm
  vg paths -x "${GRAPH_XG}" -S "${REF_SAMPLE}" -F \
    | sed "s/^>${CHROM_PREFIX}/>/" > "${IN_OUT}/${REF_SAMPLE}.fa"
  samtools faidx "${IN_OUT}/${REF_SAMPLE}.fa"

  log "step 1 finished"
}

###############################################################################
# Step 2. SV decomposition and loci
###############################################################################

sv_decompose_chr() {
  local CHR=$1
  local OUT="${SV_OUT}/${CHR}"
  mkdir -p "${OUT}" "${TMP_DIR}/sv.${CHR}"

  # 2-1) realignment and decomposition with vcfwave, split into bi-allelic
  #      records
  bcftools view -r "${CHR}" "${IN_OUT}/sv_input.vcf.gz" -Ov \
    | vcfwave -L "${WAVE_MAXLEN}" --quiet \
    | bcftools norm -m- -Ou \
    | bcftools view -e 'ALT="*"' -Ou \
    | bcftools sort -m "${SORT_MEM}" -T "${TMP_DIR}/sv.${CHR}/sort" -Ou \
    | bcftools +fill-tags -Ou -- -t AC,AN \
    | bcftools view -c 1 -Oz -o "${OUT}/${CHR}.bi.vcf.gz"
  bcftools index -f -t "${OUT}/${CHR}.bi.vcf.gz"

  # 2-2) SV alleles, |len(ALT) - len(REF)| > SV_MIN_LEN
  bcftools view -i "abs(strlen(REF)-strlen(ALT))>${SV_MIN_LEN}" \
    "${OUT}/${CHR}.bi.vcf.gz" -Oz -o "${OUT}/${CHR}.sv.bi.vcf.gz"
  bcftools index -f -t "${OUT}/${CHR}.sv.bi.vcf.gz"

  # 2-3) SV loci: overlapping SV alleles merged into one multi-allelic record
  bcftools view "${OUT}/${CHR}.sv.bi.vcf.gz" -Ov \
    | ( ulimit -v $(( MULTI_MEM_GB * 1024 * 1024 )); vcfcreatemulti --quiet ) \
    | bgzip -c > "${OUT}/${CHR}.sv.multi.vcf.gz"
  bcftools index -f -t "${OUT}/${CHR}.sv.multi.vcf.gz"

  rm -rf "${TMP_DIR}/sv.${CHR}"
}

sv_decompose() {
  log "step 2: SV decomposition and loci"
  run_jobs sv_decompose_chr ${CHROMS}
  log "step 2 finished"
}

###############################################################################
# Step 3. SV allele merging
###############################################################################

sv_merge_chr() {
  local CHR=$1
  local OUT="${SV_OUT}/${CHR}"
  mkdir -p "${TMP_DIR}/merge.${CHR}"

  # 3-1) similar alleles of each locus clustered with PanSVMerger (vsearch
  #      identity >= 0.9, alleles differing by more than 50 bp kept apart);
  #      genotypes, AC and AN are recomputed from the clusters
  python3 "${SCRIPT_DIR}/sv_merge_alleles.py" \
    --in "${OUT}/${CHR}.sv.multi.vcf.gz" \
    --out "${OUT}/${CHR}.sv.merged.vcf.gz" \
    --pansvmerger "${PANSVMERGER}" \
    --identity "${MERGE_IDENTITY}" \
    --workers "${MERGE_WORKERS}" \
    --tmpdir "${TMP_DIR}/merge.${CHR}"

  rm -rf "${TMP_DIR}/merge.${CHR}"
}

sv_merge() {
  log "step 3: SV allele merging"
  run_jobs sv_merge_chr ${CHROMS}
  log "step 3 finished"
}

###############################################################################
# Step 4. Small variant decomposition, per sample
###############################################################################

small_decompose_sample() {
  local S=$1
  local OUT="${SMALL_OUT}/samples/${S}"
  local T="${TMP_DIR}/small.${S}"
  mkdir -p "${OUT}" "${T}"

  # 4-1) sites carried by the sample, realigned with vcfwave, normalised to the
  #      reference and split into bi-allelic records
  bcftools view -a -I -s "${S}" -c 1:nref "${IN_OUT}/all.vcf.gz" -Ov \
    | vcfwave -L "${WAVE_MAXLEN}" --quiet \
    | bcftools norm -c s -f "${IN_OUT}/${REF_SAMPLE}.fa" -m -any -Oz -o "${T}/${S}.wave.vcf.gz"

  # 4-2) MNPs and complex indels broken into SNVs and simple indels; the graph
  #      path annotation (INFO/AT) is dropped first
  bcftools view -i "ALT!=\"*\" && INFO/INV!=1 && strlen(REF)<=${RTG_MAXLEN} && strlen(ALT)<=${RTG_MAXLEN}" \
    "${T}/${S}.wave.vcf.gz" -Ou \
    | bcftools annotate -x INFO/AT -Oz -o "${T}/${S}.pre.vcf.gz"
  rtg RTG_MEM="${RTG_MEM}" vcfdecompose -i "${T}/${S}.pre.vcf.gz" -o "${T}/${S}.decomp.vcf.gz" \
    --break-mnps --break-indels

  # 4-3) small variants, |len(ALT) - len(REF)| <= SV_MIN_LEN; records of the
  #      same allele on the two haplotypes merged into one record
  bcftools view -i "ALT!=\"*\" && abs(strlen(REF)-strlen(ALT))<=${SV_MIN_LEN}" \
    "${T}/${S}.decomp.vcf.gz" -Ou \
    | bcftools sort -m "${SORT_MEM}" -T "${T}/sort" -Ov \
    | python3 "${SCRIPT_DIR}/dedup_sample.py" --report "${OUT}/${S}.dedup.txt" \
    | bcftools view -Oz -o "${OUT}/${S}.small.vcf.gz"
  bcftools index -f -t "${OUT}/${S}.small.vcf.gz"

  rm -rf "${T}"
}

small_decompose() {
  log "step 4: small variant decomposition"
  run_jobs small_decompose_sample $(cat "${IN_OUT}/samples.txt")
  log "step 4 finished"
}

###############################################################################
# Step 5. Small variant merging, per chromosome
###############################################################################

small_merge_chr() {
  local CHR=$1
  local OUT="${SMALL_OUT}/merged/${CHR}"
  mkdir -p "${OUT}"

  awk -v d="${SMALL_OUT}/samples" '{print d "/" $1 "/" $1 ".small.vcf.gz"}' \
    "${IN_OUT}/samples.txt" > "${OUT}/sample_vcfs.txt"

  # 5-1) samples merged; a sample without a record carries the reference allele
  #      there, since each sample file holds only its non-reference records
  bcftools merge -m none --missing-to-ref -r "${CHR}" -l "${OUT}/sample_vcfs.txt" -Ou \
    | bcftools +fill-tags -Ou -- -t AC,AN,AF \
    | bcftools view -c 1 -Oz -o "${OUT}/${CHR}.small.vcf.gz"
  bcftools index -f -t "${OUT}/${CHR}.small.vcf.gz"

  # 5-2) SNVs and indels
  bcftools view -v snps "${OUT}/${CHR}.small.vcf.gz" -Oz -o "${OUT}/${CHR}.snv.vcf.gz"
  bcftools view -v indels "${OUT}/${CHR}.small.vcf.gz" -Oz -o "${OUT}/${CHR}.indel.vcf.gz"
  bcftools index -f -t "${OUT}/${CHR}.snv.vcf.gz"
  bcftools index -f -t "${OUT}/${CHR}.indel.vcf.gz"
}

small_merge() {
  log "step 5: small variant merging"
  run_jobs small_merge_chr ${CHROMS}
  log "step 5 finished"
}

###############################################################################
# Step 6. Group comparison counts
###############################################################################

# Genotype table of a VCF, the input of count_groups.py.
gt_table() { bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$1"; }

count_chr() {
  local CHR=$1
  local SV="${SV_OUT}/${CHR}/${CHR}.sv.merged.vcf.gz"
  local SM="${SMALL_OUT}/merged/${CHR}"
  local ARGS=(--pop "${POP_TABLE}" --target "${TARGET_POP}")

  # 6-1) target-specific, shared and other-specific variants
  #        SV     locus  : vcfcreatemulti locus
  #               allele : PanSVMerger allele cluster
  #        small  locus  : variants with overlapping reference spans
  #               allele : identical CHROM, POS, REF and ALT
  {
    gt_table "${SV}"                   | python3 "${SCRIPT_DIR}/count_groups.py" --mode locus   --varset sv    "${ARGS[@]}"
    gt_table "${SM}/${CHR}.snv.vcf.gz"   | python3 "${SCRIPT_DIR}/count_groups.py" --mode exact   --varset small "${ARGS[@]}"
    gt_table "${SM}/${CHR}.indel.vcf.gz" | python3 "${SCRIPT_DIR}/count_groups.py" --mode exact   --varset small "${ARGS[@]}"
    gt_table "${SM}/${CHR}.small.vcf.gz" | python3 "${SCRIPT_DIR}/count_groups.py" --mode overlap --varset small "${ARGS[@]}"
  } | awk -v c="${CHR}" 'BEGIN {OFS = "\t"} {print c, $0}' > "${COUNT_OUT}/${CHR}.tsv"
}

count() {
  log "step 6: group comparison counts"
  run_jobs count_chr ${CHROMS}
  log "step 6 finished"
}

###############################################################################
# Step 7. Summary tables
###############################################################################

summary() {
  log "step 7: summary tables"

  # 7-1) counts summed over the chromosomes, novel % = target-specific /
  #      (shared + other-specific)
  python3 "${SCRIPT_DIR}/summarize_variant_comparison.py" \
    --counts "${COUNT_OUT}" --out "${CMP_OUT}" \
    --target-label "${TARGET_LABEL}" --other-label "${OTHER_LABEL}" \
    | tee -a "${LOG_OUT}/pipeline.log"

  log "step 7 finished"
}

###############################################################################
# Entry point
###############################################################################

usage() {
  cat <<EOF
Usage: bash $(basename "$0") <step>

  prepare_input     vcfbub, contig renaming, reference haplotype, centromere
                    and telomere removal, reference FASTA
  sv_decompose      vcfwave, SV alleles (> ${SV_MIN_LEN} bp), loci with vcfcreatemulti
  sv_merge          PanSVMerger allele clustering within each locus
  small_decompose   per sample: vcfwave, bcftools norm, rtg vcfdecompose
  small_merge       per chromosome: sample merge, SNV and indel files
  count             ${TARGET_LABEL}-specific, shared and ${OTHER_LABEL}-specific variants
  summary           summary tables at the locus and allele level
  all               run every step in order
EOF
}

main() {
  case "${1:-}" in
    prepare_input|sv_decompose|sv_merge|small_decompose|small_merge|count|summary|all)
      setup_dirs ;;
  esac
  case "${1:-}" in
    prepare_input)   prepare_input ;;
    sv_decompose)    sv_decompose ;;
    sv_merge)        sv_merge ;;
    small_decompose) small_decompose ;;
    small_merge)     small_merge ;;
    count)           count ;;
    summary)         summary ;;
    all)
      prepare_input
      sv_decompose
      sv_merge
      small_decompose
      small_merge
      count
      summary
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
