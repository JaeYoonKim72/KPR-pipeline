#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-8}

DATA_DIR=${DATA_DIR:?set DATA_DIR to the directory containing the short-read FASTQ files}
REF_DIR=${REF_DIR:?set REF_DIR to the directory containing the reference files}
GRAPH_DIR=${GRAPH_DIR:?set GRAPH_DIR to the directory containing the pangenome graphs}
OUT_DIR=${OUT_DIR:-$(pwd)/results}
LIST_DIR=${LIST_DIR:-$(pwd)/lists}
TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}

# Pangenome references to evaluate. For each name the following files are
# expected under ${GRAPH_DIR}/<name>/ :
#   <name>.gbz      graph used by vg giraffe, vg pack and vg call
#   <name>.xg       graph used by vg surject and vg paths
#   <name>.hapl     haplotype index used for haplotype sampling
#   <name>.snarls   snarls used by vg call, generated here if absent
GRAPHS=${GRAPHS:-"kpr hprc hprc_kpr"}

# Name of the linear reference embedded in the graph. Alignments are projected
# onto its paths and variants are reported in its coordinates.
REF_SAMPLE=${REF_SAMPLE:-GRCh38}

# k-mer counting parameters used to build the per-sample KFF index.
KMC_K=${KMC_K:-29}
KMC_MEM=${KMC_MEM:-128}

# Size threshold separating indels from structural variants.
SV_MIN_LEN=${SV_MIN_LEN:-50}

# Genotype filters applied to the classified call sets.
MIN_DP=${MIN_DP:-7}
MIN_GQ=${MIN_GQ:-10}
MIN_QUAL=${MIN_QUAL:-30}

# DeepVariant v1.6 is run through the Parabricks container.
DV_IMAGE=${DV_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1}

# GRCh38 primary assembly, used for the linear reference baseline
#   https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
GRCH38=${GRCH38:-${REF_DIR}/GRCh38.fa}

BQSR_KNOWN_SITES=${BQSR_KNOWN_SITES:-}
BQSR_JAVA_MEM=${BQSR_JAVA_MEM:-16g}

# One sample identifier per line. For every identifier the reads are expected
# at ${DATA_DIR}/<sample>_1.fastq.gz and ${DATA_DIR}/<sample>_2.fastq.gz
SAMPLE_LIST=${SAMPLE_LIST:-${LIST_DIR}/samples.list}

READ_OUT=${OUT_DIR}/reads
INDEX_OUT=${OUT_DIR}/kmer_index
REF_OUT=${OUT_DIR}/reference
GAF_OUT=${OUT_DIR}/gaf
BAM_OUT=${OUT_DIR}/bam
SNV_OUT=${OUT_DIR}/small_variants
SV_OUT=${OUT_DIR}/structural_variants
CLASS_OUT=${OUT_DIR}/variant_classes
LINEAR_OUT=${OUT_DIR}/grch38
LOG_OUT=${OUT_DIR}/logs

mkdir -p "${READ_OUT}" "${INDEX_OUT}" "${REF_OUT}" "${GAF_OUT}" "${BAM_OUT}" \
         "${SNV_OUT}" "${SV_OUT}" "${CLASS_OUT}" "${LINEAR_OUT}" "${LOG_OUT}" \
         "${TMP_DIR}" "${LIST_DIR}"

log() { echo -e "[$(date '+%F %T')] $*" | tee -a "${LOG_OUT}/pipeline.log"; }

verify_bqsr_inputs() {
  command -v gatk >/dev/null 2>&1 || { echo "[ERROR] gatk is required for BQSR" >&2; return 1; }
  command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools is required for BQSR" >&2; return 1; }
  [[ -n "${BQSR_KNOWN_SITES//[[:space:]]/}" ]] || {
    echo "[ERROR] BQSR_KNOWN_SITES must specify indexed GRCh38 known-sites VCF(s)" >&2
    return 1
  }
  local FILE
  for FILE in $BQSR_KNOWN_SITES; do
    [[ -s "$FILE" ]] || { echo "[ERROR] known-sites VCF not found: $FILE" >&2; return 1; }
    [[ -s "$FILE.tbi" || -s "$FILE.csi" ]] || {
      echo "[ERROR] known-sites VCF requires an index: $FILE" >&2
      return 1
    }
  done
}

prepare_bqsr_known_args() {
  local GRAPH="$1" FILE TARGET TARGET_DIR MAP INV CHROMS INDEX=0
  BQSR_KNOWN_ARGS=()
  if [[ -n "$GRAPH" ]]; then
    MAP="${REF_OUT}/${GRAPH}/${REF_SAMPLE}.rename_chrs.txt"
    [[ -s "$MAP" ]] || { echo "[ERROR] missing graph contig map: $MAP" >&2; return 1; }
    TARGET_DIR="${REF_OUT}/${GRAPH}/bqsr_known_sites"
    mkdir -p "$TARGET_DIR"
    INV="${TARGET_DIR}/plain_to_graph_chrs.tsv"
    if [[ ! -s "$INV" || "$MAP" -nt "$INV" ]]; then
      awk -F '\t' 'BEGIN {OFS="\t"} NF>=2 {print $2,$1}' "$MAP" > "$INV"
    fi
    [[ -s "$INV" ]] || { echo "[ERROR] empty contig mapping: $MAP" >&2; return 1; }
  fi

  for FILE in $BQSR_KNOWN_SITES; do
    TARGET="$FILE"
    if [[ -n "$GRAPH" ]]; then
      TARGET="${TARGET_DIR}/known_sites_${INDEX}.vcf.gz"
      CHROMS=$(awk 'NR==FNR {ok[$1]=1; next} ($1 in ok) {printf "%s%s",sep,$1;sep=","}' <(bcftools index -s "$FILE") "$INV")
      [[ -n "$CHROMS" ]] || {
        echo "[ERROR] known-sites VCF and graph reference have no shared GRCh38 contigs: $FILE" >&2
        return 1
      }
      if [[ ! -s "$TARGET" || "$FILE" -nt "$TARGET" || "$INV" -nt "$TARGET" ]]; then
        bcftools view -r "$CHROMS" "$FILE" -Ou \
          | bcftools annotate --rename-chrs "$INV" -Oz -o "$TARGET"
        bcftools index -f -t "$TARGET"
      fi
      [[ -s "$TARGET.tbi" || -s "$TARGET.csi" ]] || bcftools index -f -t "$TARGET"
    fi
    BQSR_KNOWN_ARGS+=(--known-sites "$TARGET")
    INDEX=$((INDEX + 1))
  done
}

run_bqsr() {
  local PREFIX="$1" REF="$2"
  [[ -s "$REF" && -s "${PREFIX}.rmdup.bam" ]] || {
    echo "[ERROR] missing BQSR reference or BAM: $REF / ${PREFIX}.rmdup.bam" >&2
    return 1
  }
  gatk --java-options "-Xmx${BQSR_JAVA_MEM}" BaseRecalibrator \
    -R "$REF" -I "${PREFIX}.rmdup.bam" \
    "${BQSR_KNOWN_ARGS[@]}" -O "${PREFIX}.BQSR.table"
  gatk --java-options "-Xmx${BQSR_JAVA_MEM}" ApplyBQSR \
    -R "$REF" -I "${PREFIX}.rmdup.bam" \
    --bqsr-recal-file "${PREFIX}.BQSR.table" -O "${PREFIX}.BQSR.bam"
  samtools index -@ "${THREADS}" "${PREFIX}.BQSR.bam"
}

# Step 1. Reference preparation
prepare_reference() {
  log "step 1: reference preparation"

  for G in ${GRAPHS}; do
    XG="${GRAPH_DIR}/${G}/${G}.xg"
    GBZ="${GRAPH_DIR}/${G}/${G}.gbz"
    [[ -s "${XG}" ]] || { log "  ${XG} not found, skipped"; continue; }

    OUT="${REF_OUT}/${G}"
    mkdir -p "${OUT}"

    # 1-1) paths of the linear reference embedded in the graph, used as the
    #      projection target of vg surject
    vg paths -L -x "${XG}" | grep "^${REF_SAMPLE}" | sort -V \
      > "${OUT}/${REF_SAMPLE}.paths.txt"

    # 1-2) sequence of those paths, used as the reference for variant calling
    vg paths -x "${XG}" -S "${REF_SAMPLE}" -F > "${OUT}/${REF_SAMPLE}.fa"
    samtools faidx "${OUT}/${REF_SAMPLE}.fa"
    samtools dict "${OUT}/${REF_SAMPLE}.fa" -o "${OUT}/${REF_SAMPLE}.dict"

    # 1-3) mapping from the path names of the graph to plain chromosome names,
    #      used to rename the contigs of the called VCF files
    awk -F'#' '{print $0"\t"$NF}' "${OUT}/${REF_SAMPLE}.paths.txt" \
      > "${OUT}/${REF_SAMPLE}.rename_chrs.txt"

    # 1-4) snarls, required by vg call
    if [[ ! -s "${GRAPH_DIR}/${G}/${G}.snarls" && -s "${GBZ}" ]]; then
      vg snarls -t "${THREADS}" "${GBZ}" > "${GRAPH_DIR}/${G}/${G}.snarls"
    fi
  done

  log "step 1 finished"
}

# Step 2. Read preparation
read_prep() {
  log "step 2: read preparation"

  while read -r S; do
    [[ -z "${S}" ]] && continue

    # 2-1) remove the read-group tag left in the read names by the sequencing
    #      centre, which vg giraffe does not accept
    for R in 1 2; do
      zcat "${DATA_DIR}/${S}_${R}.fastq.gz" \
        | sed 's/\sRG:Z:1//' \
        | bgzip -c > "${READ_OUT}/${S}_${R}.cleaned.fastq.gz"
    done

    # 2-2) list of the two read files, used by kmc
    {
      echo "${READ_OUT}/${S}_1.cleaned.fastq.gz"
      echo "${READ_OUT}/${S}_2.cleaned.fastq.gz"
    } > "${READ_OUT}/${S}.reads.txt"
  done < "${SAMPLE_LIST}"

  log "step 2 finished"
}

# Step 3. k-mer index of each sample
kmer_index() {
  log "step 3: k-mer counting"

  # kmc writes a KFF file that vg giraffe uses to sample the haplotypes of the
  # graph relevant to the sample, which speeds up mapping.
  while read -r S; do
    [[ -z "${S}" ]] && continue
    mkdir -p "${TMP_DIR}/kmc_${S}"

    kmc -k"${KMC_K}" -m"${KMC_MEM}" -okff -t"${THREADS}" -hp \
      "@${READ_OUT}/${S}.reads.txt" \
      "${INDEX_OUT}/${S}" \
      "${TMP_DIR}/kmc_${S}" \
      > "${LOG_OUT}/${S}.kmc.log" 2>&1

    rm -rf "${TMP_DIR}/kmc_${S}"
  done < "${SAMPLE_LIST}"

  log "step 3 finished"
}

# Step 4. Read mapping to the pangenome graph
giraffe() {
  log "step 4: vg giraffe mapping"

  # vg v1.63.1
  for G in ${GRAPHS}; do
    GBZ="${GRAPH_DIR}/${G}/${G}.gbz"
    HAPL="${GRAPH_DIR}/${G}/${G}.hapl"
    [[ -s "${GBZ}" ]] || { log "  ${GBZ} not found, skipped"; continue; }

    mkdir -p "${GAF_OUT}/${G}"

    while read -r S; do
      [[ -z "${S}" ]] && continue

      vg giraffe -Z "${GBZ}" \
        -f "${READ_OUT}/${S}_1.cleaned.fastq.gz" \
        -f "${READ_OUT}/${S}_2.cleaned.fastq.gz" \
        -o gaf -p -t "${THREADS}" \
        --sample "${S}" --progress \
        --kff-name "${INDEX_OUT}/${S}.kff" \
        --haplotype-name "${HAPL}" \
        2> "${LOG_OUT}/${G}.${S}.giraffe.err" \
        | bgzip -c > "${GAF_OUT}/${G}/${G}.${S}.gaf.gz"
    done < "${SAMPLE_LIST}"
  done

  log "step 4 finished"
}

# Step 5. Projection of the alignments onto the linear reference
surject() {
  log "step 5: vg surject"

  for G in ${GRAPHS}; do
    XG="${GRAPH_DIR}/${G}/${G}.xg"
    PATHS="${REF_OUT}/${G}/${REF_SAMPLE}.paths.txt"
    [[ -s "${XG}" && -s "${PATHS}" ]] || { log "  inputs for ${G} not found, skipped"; continue; }

    mkdir -p "${BAM_OUT}/${G}"

    while read -r S; do
      [[ -z "${S}" ]] && continue

      vg surject -x "${XG}" \
        -G "${GAF_OUT}/${G}/${G}.${S}.gaf.gz" \
        --interleaved \
        -F "${PATHS}" \
        -b -N "${S}" \
        -R "ID:${S} LB:lib1 SM:${S} PL:illumina PU:unit1" \
        1> "${BAM_OUT}/${G}/${G}.${S}.raw.bam" \
        2> "${LOG_OUT}/${G}.${S}.surject.err"
    done < "${SAMPLE_LIST}"
  done

  log "step 5 finished"
}

# Step 6. BAM sorting, duplicate removal and base-quality recalibration
bam_cleanup() {
  log "step 6: BAM cleanup and GATK BQSR"
  verify_bqsr_inputs

  for G in $GRAPHS; do
    [[ -d "${BAM_OUT}/$G" ]] || { log "  ${BAM_OUT}/$G not found, skipped"; continue; }
    REF="${REF_OUT}/$G/${REF_SAMPLE}.fa"
    [[ -s "$REF" ]] || { log "  $REF not found, skipped"; continue; }
    prepare_bqsr_known_args "$G"

    while read -r S; do
      [[ -z "$S" ]] && continue
      B="${BAM_OUT}/$G/$G.$S"

      samtools sort -O BAM --threads "${THREADS}" -T "${TMP_DIR}/$G.$S.tmp" \
        -o "${B}.sorted.bam" "${B}.raw.bam"
      samtools rmdup "${B}.sorted.bam" "${B}.rmdup.bam" \
        2> "${LOG_OUT}/$G.$S.rmdup.err"
      samtools index -@ "${THREADS}" "${B}.rmdup.bam"

      run_bqsr "$B" "$REF"
      rm -f "${B}.sorted.bam"
    done < "${SAMPLE_LIST}"
  done

  log "step 6 finished"
}

# Step 7. Small variant calling
small_variants() {
  log "step 7: small variant calling"

  for G in ${GRAPHS}; do
    REF="${REF_OUT}/${G}/${REF_SAMPLE}.fa"
    RENAME="${REF_OUT}/${G}/${REF_SAMPLE}.rename_chrs.txt"
    [[ -s "${REF}" ]] || { log "  ${REF} not found, skipped"; continue; }

    OUT="${SNV_OUT}/${G}"
    mkdir -p "${OUT}"

    while read -r S; do
      [[ -z "${S}" ]] && continue
      V="${OUT}/${G}.${S}"

      # 7-1) DeepVariant v1.6, run through the Parabricks container. The CPU
      #      build can be substituted here with the same reference and BAM.
      docker run --rm --gpus all \
        --volume "$(pwd)":/workdir --volume "$(pwd)":/outputdir \
        --workdir /workdir \
        "${DV_IMAGE}" \
        pbrun deepvariant \
        --ref "${REF}" \
        --in-bam "${BAM_OUT}/${G}/${G}.${S}.BQSR.bam" \
        --out-variants "${V}.dp.vcf.gz" \
        --logfile "${LOG_OUT}/${G}.${S}.deepvariant.log" \
        --sort-by-haplotypes --parse-sam-aux-fields --mode shortread \
        >& "${LOG_OUT}/${G}.${S}.deepvariant.err"

      # 7-2) keep the PASS sites and rename the contigs from the path names of
      #      the graph to plain chromosome names
      bcftools view -f PASS "${V}.dp.vcf.gz" \
        | bcftools annotate --rename-chrs "${RENAME}" \
          -Oz -o "${V}.PASS.vcf.gz"
      bcftools index -f -t "${V}.PASS.vcf.gz"

      # 7-3) drop the homozygous reference genotypes
      bcftools view -i 'GT!="0/0" & GT!="0|0"' "${V}.PASS.vcf.gz" \
        -Oz -o "${V}.PASS.NoZero.vcf.gz"
      bcftools index -f -t "${V}.PASS.NoZero.vcf.gz"
    done < "${SAMPLE_LIST}"
  done

  log "step 7 finished"
}

# Step 8. Structural variant calling
sv_variants() {
  log "step 8: structural variant calling"

  for G in ${GRAPHS}; do
    GBZ="${GRAPH_DIR}/${G}/${G}.gbz"
    SNARLS="${GRAPH_DIR}/${G}/${G}.snarls"
    FAI="${REF_OUT}/${G}/${REF_SAMPLE}.fa.fai"
    [[ -s "${GBZ}" ]] || { log "  ${GBZ} not found, skipped"; continue; }

    OUT="${SV_OUT}/${G}"
    mkdir -p "${OUT}"

    while read -r S; do
      [[ -z "${S}" ]] && continue
      V="${OUT}/${G}.${S}"

      # 8-1) read support of every node and edge of the graph
      vg pack -x "${GBZ}" -Q 5 -t "${THREADS}" \
        -a "${GAF_OUT}/${G}/${G}.${S}.gaf.gz" \
        -o "${V}.pack" \
        2> "${LOG_OUT}/${G}.${S}.pack.err"

      # 8-2) genotype the bubbles and the alternative paths of the graph
      vg call "${GBZ}" \
        -r "${SNARLS}" \
        -k "${V}.pack" \
        -s "${S}" -S "${REF_SAMPLE}" \
        -a -z -t "${THREADS}" \
        2> "${LOG_OUT}/${G}.${S}.call.err" \
        | bgzip -c > "${V}.sv.vcf.gz"
      bcftools index -f "${V}.sv.vcf.gz"

      # 8-3) keep the PASS sites and align the header with the reference index
      bcftools view -i 'FILTER="PASS"' "${V}.sv.vcf.gz" \
        -Oz -o "${TMP_DIR}/${G}.${S}.sv.tmp.vcf.gz"
      bcftools reheader --fai "${FAI}" \
        -o "${V}.sv.PASS.vcf.gz" "${TMP_DIR}/${G}.${S}.sv.tmp.vcf.gz"
      bcftools index -f -t "${V}.sv.PASS.vcf.gz"
      rm -f "${TMP_DIR}/${G}.${S}.sv.tmp.vcf.gz"

      # 8-4) drop the homozygous reference genotypes, and then the missing ones
      bcftools view -e 'GT="0/0"' "${V}.sv.PASS.vcf.gz" \
        -Oz -o "${V}.sv.PASS.NoZero.vcf.gz"
      bcftools index -f -t "${V}.sv.PASS.NoZero.vcf.gz"

      bcftools view -e 'GT="0/0" | GT="./."' "${V}.sv.PASS.NoZero.vcf.gz" \
        -Oz -o "${V}.sv.PASS.NoZeroNoMiss.vcf.gz"
      bcftools index -f -t "${V}.sv.PASS.NoZeroNoMiss.vcf.gz"
    done < "${SAMPLE_LIST}"
  done

  log "step 8 finished"
}

# Step 9. Variant classification and genotype filtering
# Split one call set into SNPs, bi-allelic SNPs, indels and structural
# variants, and the structural variants into insertions and deletions.
#   $1 = input vcf.gz   $2 = output prefix
classify_vcf() {
  local VCF=$1 OUT=$2

  bcftools view -v snps "${VCF}" \
    | bcftools norm -m +any \
    | bcftools view -V indels -i 'strlen(REF)==1 & strlen(ALT)==1' \
    | bcftools norm -d all - -Oz -o "${OUT}.snp.vcf.gz"

  bcftools view -v snps "${VCF}" \
    | bcftools norm -m +any \
    | bcftools view -V indels -M 2 -m 2 -i 'strlen(REF)==1 & strlen(ALT)==1' \
    | bcftools norm -d all - -Oz -o "${OUT}.bisnp.vcf.gz"

  bcftools view -v indels "${VCF}" \
    | bcftools view -i "strlen(REF) != strlen(ALT) & abs(strlen(REF) - strlen(ALT)) <= ${SV_MIN_LEN}" \
    | bcftools norm -m +any \
    | bcftools view -V snps \
    | bcftools norm -d all - -Oz -o "${OUT}.indel.vcf.gz"

  bcftools view -v indels "${VCF}" \
    | bcftools view -i "strlen(REF) != strlen(ALT) & abs(strlen(REF) - strlen(ALT)) > ${SV_MIN_LEN}" \
    | bcftools norm -m +any \
    | bcftools view -V snps \
    | bcftools norm -d all - -Oz -o "${OUT}.sv.vcf.gz"

  bcftools view -i 'strlen(ALT) > strlen(REF)' "${OUT}.sv.vcf.gz" \
    -Oz -o "${OUT}.sv.INS.vcf.gz"
  bcftools view -i 'strlen(ALT) < strlen(REF)' "${OUT}.sv.vcf.gz" \
    -Oz -o "${OUT}.sv.DEL.vcf.gz"
}

classify() {
  log "step 9: variant classification"

  for G in ${GRAPHS}; do
    OUT="${CLASS_OUT}/${G}"
    mkdir -p "${OUT}"/{class,selected}

    # 9-1) classify the small variant and the structural variant call sets
    for VCF in "${SNV_OUT}/${G}"/*.PASS.vcf.gz "${SNV_OUT}/${G}"/*.PASS.NoZero.vcf.gz \
               "${SV_OUT}/${G}"/*.sv.PASS.vcf.gz "${SV_OUT}/${G}"/*.sv.PASS.NoZero.vcf.gz; do
      [[ -s "${VCF}" ]] || continue
      B=$(basename "${VCF}" .vcf.gz)
      classify_vcf "${VCF}" "${OUT}/class/${B}"
    done

    # 9-2) keep the genotypes supported by enough depth and quality, and count
    #      the remaining records
    {
      echo -e "file\trecords"
      for VCF in "${OUT}"/class/*.vcf.gz; do
        B=$(basename "${VCF}" .vcf.gz)
        bcftools view -i "FORMAT/DP>=${MIN_DP} & FORMAT/GQ>=${MIN_GQ} & QUAL>=${MIN_QUAL}" \
          "${VCF}" -Oz -o "${OUT}/selected/${B}.sel.vcf.gz"
        bcftools index -f "${OUT}/selected/${B}.sel.vcf.gz"
        echo -e "${B}.sel\t$(bcftools index -n "${OUT}/selected/${B}.sel.vcf.gz")"
      done
    } > "${CLASS_OUT}/${G}.variant_counts.tsv"
  done

  log "step 9 finished"
}

# Step 10. Linear reference baseline
grch38_baseline() {
  log "step 10: GRCh38 baseline"

  mkdir -p "${LINEAR_OUT}"/{bam,small_variants,sv}
  verify_bqsr_inputs
  [[ -s "${GRCH38}" ]] || { echo "[ERROR] missing GRCh38 FASTA: ${GRCH38}" >&2; return 1; }
  prepare_bqsr_known_args ""

  while read -r S; do
    [[ -z "${S}" ]] && continue

    # 10-1) alignment to GRCh38
    bwa-mem2 mem -t "${THREADS}" \
      -R "@RG\tID:${S}\tSM:${S}\tLB:lib1\tPL:illumina\tPU:unit1" \
      "${GRCH38}" \
      "${READ_OUT}/${S}_1.cleaned.fastq.gz" \
      "${READ_OUT}/${S}_2.cleaned.fastq.gz" \
      | samtools sort -O BAM --threads "${THREADS}" -T "${TMP_DIR}/${S}.grch38.tmp" \
        -o "${LINEAR_OUT}/bam/${S}.sorted.bam"

    samtools rmdup "${LINEAR_OUT}/bam/${S}.sorted.bam" \
                   "${LINEAR_OUT}/bam/${S}.rmdup.bam" \
      2> "${LOG_OUT}/${S}.grch38.rmdup.err"
    samtools index -@ "${THREADS}" "${LINEAR_OUT}/bam/${S}.rmdup.bam"
    run_bqsr "${LINEAR_OUT}/bam/${S}" "${GRCH38}"
    rm -f "${LINEAR_OUT}/bam/${S}.sorted.bam"

    # 10-2) small variants, DeepVariant v1.6
    V="${LINEAR_OUT}/small_variants/${S}"
    docker run --rm --gpus all \
      --volume "$(pwd)":/workdir --volume "$(pwd)":/outputdir \
      --workdir /workdir \
      "${DV_IMAGE}" \
      pbrun deepvariant \
      --ref "${GRCH38}" \
      --in-bam "${LINEAR_OUT}/bam/${S}.BQSR.bam" \
      --out-variants "${V}.dp.vcf.gz" \
      --logfile "${LOG_OUT}/${S}.grch38.deepvariant.log" \
      --sort-by-haplotypes --parse-sam-aux-fields --mode shortread \
      >& "${LOG_OUT}/${S}.grch38.deepvariant.err"

    bcftools view -f PASS "${V}.dp.vcf.gz" -Oz -o "${V}.PASS.vcf.gz"
    bcftools index -f -t "${V}.PASS.vcf.gz"

    # 10-3) structural variants from the short reads, Manta and Delly, with the
    #       same filters as the assembly pipeline
    configManta.py --bam "${LINEAR_OUT}/bam/${S}.BQSR.bam" \
                   --referenceFasta "${GRCH38}" \
                   --runDir "${LINEAR_OUT}/sv/manta_${S}"
    "${LINEAR_OUT}/sv/manta_${S}/runWorkflow.py" -j "${THREADS}"

    bcftools view -i 'QUAL >= 30 && (INFO/SVTYPE=="INS" || INFO/SVTYPE=="DUP")' \
      -Oz -o "${LINEAR_OUT}/sv/${S}.manta.ins_dup.vcf.gz" \
      "${LINEAR_OUT}/sv/manta_${S}/results/variants/diploidSV.vcf.gz"

    delly call -g "${GRCH38}" -o "${LINEAR_OUT}/sv/${S}.delly.bcf" \
               "${LINEAR_OUT}/bam/${S}.BQSR.bam"
    bcftools view -f PASS -i 'INFO/SVTYPE=="DEL"' \
      -Oz -o "${LINEAR_OUT}/sv/${S}.delly.del.vcf.gz" \
      "${LINEAR_OUT}/sv/${S}.delly.bcf"
  done < "${SAMPLE_LIST}"

  log "step 10 finished"
}

usage() {
  cat <<EOF
Usage: bash $(basename "$0") <step>

  prepare_reference   reference paths, FASTA, dictionary and snarls per graph
  read_prep           read name clean-up and read lists
  kmer_index          kmc k-mer index used for haplotype sampling
  giraffe             read mapping to the pangenome graph
  surject             projection of the alignments onto the linear reference
  bam_cleanup         sort, remove duplicates, recalibrate base qualities (GATK BQSR)
  small_variants      DeepVariant calling and PASS filtering
  sv_variants         vg pack and vg call genotyping and PASS filtering
  classify            classification into SNP, indel and SV, and counting
  grch38_baseline     alignment, duplicate removal, BQSR, DeepVariant, Manta and Delly
  all                 run every step in order

  Set BQSR_KNOWN_SITES to indexed GRCh38 known-sites VCF(s).
  For graph-projected BAMs, the VCF contig names are mapped to graph paths.
  BQSR_JAVA_MEM controls the GATK Java heap (default: 16g).
EOF
}

main() {
  case "${1:-}" in
    prepare_reference) prepare_reference ;;
    read_prep)         read_prep ;;
    kmer_index)        kmer_index ;;
    giraffe)           giraffe ;;
    surject)           surject ;;
    bam_cleanup)       bam_cleanup ;;
    small_variants)    small_variants ;;
    sv_variants)       sv_variants ;;
    classify)          classify ;;
    grch38_baseline)   grch38_baseline ;;
    all)
      prepare_reference
      read_prep
      kmer_index
      giraffe
      surject
      bam_cleanup
      small_variants
      sv_variants
      classify
      grch38_baseline
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
