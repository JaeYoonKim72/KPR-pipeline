#!/usr/bin/env bash
#
# Constructing a Korean pangenome reveals distinct haplotypes at the amylase
# locus
#
# genome_assembly_and_statistics.sh
#
# Usage: bash genome_assembly_and_statistics.sh <step>
#        (see README.md for the step list and for the required inputs)

set -euo pipefail

###############################################################################
# Configuration
###############################################################################

THREADS=${THREADS:-32}
GENOME_SIZE=${GENOME_SIZE:-3200000000}     # 3.2 Gb, used for raw depth only

DATA_DIR=${DATA_DIR:?set DATA_DIR to the directory containing the input FASTQ files}
OUT_DIR=${OUT_DIR:-$(pwd)/results}
REF_DIR=${REF_DIR:?set REF_DIR to the directory containing the reference files}
LIST_DIR=${LIST_DIR:-$(pwd)/lists}

# --- reference files ---------------------------------------------------------
# GRCh38 primary assembly (soft-masked)
#   https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
GRCH38=${GRCH38:-${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa}

# Ensembl Regulatory Build, used by QUAST-LG as the feature track
#   https://ftp.ensembl.org/pub/release-108/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz
GRCH38_GFF=${GRCH38_GFF:-${REF_DIR}/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.gff.gz}

# T2T-CHM13 v2.0
#   https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
CHM13=${CHM13:-${REF_DIR}/chm13v2.0.fa}

# GENCODE v38 annotation
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
GENCODE_GTF=${GENCODE_GTF:-${REF_DIR}/gencode.v38.annotation.gtf}

# Augmented CHM13 satellite library used in the second RepeatMasker round
#   https://github.com/marbl/CHM13  (T2T-CHM13 repeat annotation resources)
SAT_LIB=${SAT_LIB:-${REF_DIR}/chm13_satellite_library.fa}

# Kraken2 database for contaminant screening of the LCL-derived assemblies
#   https://benlangmead.github.io/aws-indexes/k2
KRAKEN_DB=${KRAKEN_DB:-${REF_DIR}/kraken2_db}

# --- sample lists (one sample ID per line) -----------------------------------
HIFI_LIST=${HIFI_LIST:-${LIST_DIR}/hifi.list}         # PacBio HiFi samples
ONT_LIST=${ONT_LIST:-${LIST_DIR}/ont.list}            # ONT PromethION samples
HIC_LIST=${HIC_LIST:-${LIST_DIR}/hic.list}            # samples with Hi-C data
QC_PASS_LIST=${QC_PASS_LIST:-${LIST_DIR}/qc_pass.list}  # written by asm_stats

# --- assembly quality thresholds ---------------------------------------------
MIN_N50=${MIN_N50:-20000000}      # contig N50 >= 20 Mb
MAX_CONTIG=${MAX_CONTIG:-1500}    # number of contigs <= 1,500

QC_OUT=${OUT_DIR}/read_qc
ASM_OUT=${OUT_DIR}/assembly
POLISH_OUT=${OUT_DIR}/polishing
STAT_OUT=${OUT_DIR}/statistics
ANNOT_OUT=${OUT_DIR}/annotation
SV_OUT=${OUT_DIR}/sv
LOG_OUT=${OUT_DIR}/logs

mkdir -p "${QC_OUT}" "${ASM_OUT}" "${POLISH_OUT}" "${STAT_OUT}" \
         "${ANNOT_OUT}" "${SV_OUT}" "${LOG_OUT}" "${LIST_DIR}"

log() { echo -e "[$(date '+%F %T')] $*" | tee -a "${LOG_OUT}/pipeline.log"; }

###############################################################################
# Step 1. Long-read and Hi-C read quality control
###############################################################################

calc_fastq_metrics() {
  local FQ=$1
  local READS BASES

  if [[ ${FQ} == *.gz ]]; then
    READS=$(zcat "${FQ}" | awk 'END{print NR/4}')
    BASES=$(zcat "${FQ}" | awk 'NR%4==2{t+=length($0)} END{print t}')
  else
    READS=$(awk 'END{print NR/4}' "${FQ}")
    BASES=$(awk 'NR%4==2{t+=length($0)} END{print t}' "${FQ}")
  fi

  awk -v f="$(basename "${FQ}")" -v r="${READS}" -v b="${BASES}" -v g="${GENOME_SIZE}" \
    'BEGIN{printf "%s\t%d\t%d\t%.2f\n", f, r, b, b/g}'
}

read_qc() {
  log "step 1: long-read and Hi-C read QC"
  mkdir -p "${QC_OUT}"/{fastqc,trimmed,stats}

  echo -e "file\traw_reads\traw_bases\traw_depth" > "${QC_OUT}/stats/read_metrics.tsv"

  # 1-1) PacBio HiFi: FastQC v0.12.1, HiFiAdapterFilt v3.0
  while read -r S; do
    [[ -z "${S}" ]] && continue

    fastqc -t "${THREADS}" -o "${QC_OUT}/fastqc" "${DATA_DIR}/${S}.hifi.fastq.gz"
    calc_fastq_metrics "${DATA_DIR}/${S}.hifi.fastq.gz" >> "${QC_OUT}/stats/read_metrics.tsv"

    hifiadapterfilt.sh -p "${DATA_DIR}/${S}.hifi" -t "${THREADS}" -o "${QC_OUT}/trimmed"
    calc_fastq_metrics "${QC_OUT}/trimmed/${S}.hifi.filt.fastq.gz" >> "${QC_OUT}/stats/read_metrics.tsv"
  done < "${HIFI_LIST}"

  # 1-2) ONT PromethION: Porechop v0.2.4, Trimmomatic v0.39
  while read -r S; do
    [[ -z "${S}" ]] && continue

    fastqc -t "${THREADS}" -o "${QC_OUT}/fastqc" "${DATA_DIR}/${S}.ont.fastq.gz"
    calc_fastq_metrics "${DATA_DIR}/${S}.ont.fastq.gz" >> "${QC_OUT}/stats/read_metrics.tsv"

    porechop -i "${DATA_DIR}/${S}.ont.fastq.gz" \
             -o "${QC_OUT}/trimmed/${S}.ont.porechop.fastq.gz" \
             --threads "${THREADS}"

    trimmomatic SE -threads "${THREADS}" \
      "${QC_OUT}/trimmed/${S}.ont.porechop.fastq.gz" \
      "${QC_OUT}/trimmed/${S}.ont.filt.fastq.gz" \
      SLIDINGWINDOW:100:10 AVGQUAL:20 MINLEN:80

    calc_fastq_metrics "${QC_OUT}/trimmed/${S}.ont.filt.fastq.gz" >> "${QC_OUT}/stats/read_metrics.tsv"
  done < "${ONT_LIST}"

  # 1-3) Hi-C (Dovetail Omni-C, Illumina NovaSeq X, 2 x 101 bp)
  while read -r S; do
    [[ -z "${S}" ]] && continue

    fastqc -t "${THREADS}" -o "${QC_OUT}/fastqc" \
      "${DATA_DIR}/${S}.hic_1.fastq.gz" "${DATA_DIR}/${S}.hic_2.fastq.gz"
    calc_fastq_metrics "${DATA_DIR}/${S}.hic_1.fastq.gz" >> "${QC_OUT}/stats/read_metrics.tsv"
    calc_fastq_metrics "${DATA_DIR}/${S}.hic_2.fastq.gz" >> "${QC_OUT}/stats/read_metrics.tsv"
  done < "${HIC_LIST}"

  log "step 1 finished"
}

###############################################################################
# Step 2. De novo genome assembly
###############################################################################

assembly() {
  log "step 2: de novo assembly"
  mkdir -p "${ASM_OUT}"/{hifiasm,hifiasm_hic,flye,hifiasm_ont,fasta,screen}

  # 2-1) HiFi assemblies, HiFiasm v0.24
  while read -r S; do
    [[ -z "${S}" ]] && continue
    mkdir -p "${ASM_OUT}/hifiasm/${S}"

    hifiasm -o "${ASM_OUT}/hifiasm/${S}/${S}" \
            --dual-scaf --telo-m CCCTAA \
            -t "${THREADS}" \
            "${QC_OUT}/trimmed/${S}.hifi.filt.fastq.gz"

    for H in 1 2; do
      gfatools gfa2fa "${ASM_OUT}/hifiasm/${S}/${S}.bp.hap${H}.p_ctg.gfa" \
        > "${ASM_OUT}/fasta/${S}_${H}.raw.fa"
    done
  done < "${HIFI_LIST}"

  # 2-2) haplotype-resolved assemblies, HiFiasm v0.24 in Hi-C phasing mode
  #      ultra-long ONT reads are supplied through --ul where available
  while read -r S; do
    [[ -z "${S}" ]] && continue
    mkdir -p "${ASM_OUT}/hifiasm_hic/${S}"

    UL_OPT=""
    if [[ -s "${QC_OUT}/trimmed/${S}.ont.filt.fastq.gz" ]]; then
      UL_OPT="--ul ${QC_OUT}/trimmed/${S}.ont.filt.fastq.gz"
    fi

    # shellcheck disable=SC2086
    hifiasm -o "${ASM_OUT}/hifiasm_hic/${S}/${S}" \
            --dual-scaf --telo-m CCCTAA \
            -t "${THREADS}" \
            ${UL_OPT} \
            --h1 "${DATA_DIR}/${S}.hic_1.fastq.gz" \
            --h2 "${DATA_DIR}/${S}.hic_2.fastq.gz" \
            "${QC_OUT}/trimmed/${S}.hifi.filt.fastq.gz"

    for H in 1 2; do
      gfatools gfa2fa "${ASM_OUT}/hifiasm_hic/${S}/${S}.hic.hap${H}.p_ctg.gfa" \
        > "${ASM_OUT}/fasta/${S}_${H}.raw.fa"
    done
  done < "${HIC_LIST}"

  # 2-3) PromethION assemblies, Flye v2.9.5 (R10.4.1 data basecalled with Guppy)
  #      plus HiFiasm v0.25 in ONT mode as an assembler control
  while read -r S; do
    [[ -z "${S}" ]] && continue

    flye --nano-hq "${QC_OUT}/trimmed/${S}.ont.filt.fastq.gz" \
         --out-dir "${ASM_OUT}/flye/${S}" \
         --genome-size 3g --threads "${THREADS}"
    cp "${ASM_OUT}/flye/${S}/assembly.fasta" "${ASM_OUT}/fasta/${S}_flye.raw.fa"

    hifiasm --ont -o "${ASM_OUT}/hifiasm_ont/${S}" \
            -t "${THREADS}" "${QC_OUT}/trimmed/${S}.ont.filt.fastq.gz"
  done < "${ONT_LIST}"

  # 2-4) removal of non-human contigs from the LCL-derived assemblies
  for FA in "${ASM_OUT}"/fasta/*.raw.fa; do
    B=$(basename "${FA}" .raw.fa)

    kraken2 --db "${KRAKEN_DB}" --threads "${THREADS}" \
            --output "${ASM_OUT}/screen/${B}.kraken.out" \
            --report "${ASM_OUT}/screen/${B}.kraken.report" "${FA}"

    # contigs classified as anything other than Homo sapiens (taxid 9606)
    awk '$1=="C" && $3!=9606 {print $2}' "${ASM_OUT}/screen/${B}.kraken.out" \
      > "${ASM_OUT}/screen/${B}.contaminant.ids"

    seqkit grep -v -f "${ASM_OUT}/screen/${B}.contaminant.ids" "${FA}" \
      > "${ASM_OUT}/fasta/${B}.fa"
    seqkit stat -a "${ASM_OUT}/fasta/${B}.fa" \
      >> "${ASM_OUT}/fasta/assembly_seqkit_stat.txt"
  done

  log "step 2 finished"
}

###############################################################################
# Step 3. Long-read polishing
###############################################################################

polishing() {
  log "step 3: Inspector polishing"
  mkdir -p "${POLISH_OUT}"/{eval_raw,corrected,eval_final,stat}

  polish_one() {                      # $1 = sample, $2 = haplotype, $3 = reads
    local S=$1 H=$2 READS=$3
    local ID="${S}_${H}"

    # Inspector v1.3, evaluation of the raw assembly
    inspector.py -c "${ASM_OUT}/fasta/${ID}.fa" \
                 -r "${READS}" \
                 -o "${POLISH_OUT}/eval_raw/${ID}" \
                 --datatype hifi -t "${THREADS}"

    # correction of small-scale errors only
    inspector-correct.py -i "${POLISH_OUT}/eval_raw/${ID}" \
                         --datatype pacbio-hifi \
                         -o "${POLISH_OUT}/corrected/${ID}" \
                         --skip_structural -t "${THREADS}"

    cp "${POLISH_OUT}/corrected/${ID}/contig_corrected.fa" \
       "${POLISH_OUT}/corrected/${ID}.polished.fa"

    # re-evaluation of the corrected assembly
    inspector.py -c "${POLISH_OUT}/corrected/${ID}.polished.fa" \
                 -r "${READS}" \
                 -o "${POLISH_OUT}/eval_final/${ID}" \
                 --datatype hifi -t "${THREADS}"

    cp "${POLISH_OUT}/eval_final/${ID}/summary_statistics" \
       "${POLISH_OUT}/stat/${ID}_summary_statistics"
  }

  while read -r S; do
    [[ -z "${S}" ]] && continue
    for H in 1 2; do
      polish_one "${S}" "${H}" "${QC_OUT}/trimmed/${S}.hifi.filt.fastq.gz"
    done
  done < <(cat "${HIFI_LIST}" "${HIC_LIST}" | sort -u)

  log "step 3 finished"
}

###############################################################################
# Step 4. Assembly statistics
###############################################################################

asm_stats() {
  log "step 4: assembly statistics"
  mkdir -p "${STAT_OUT}"/{quast_grch38,quast_chm13,ngx,inspector,misjoin,flagger}

  # 4-1) Inspector summary table
  {
    echo -e "assembly\tQV\tmapping_rate\tdepth\tbase_error\texpansion\tcollapse"
    for F in "${POLISH_OUT}"/stat/*_summary_statistics; do
      ID=$(basename "${F}" _summary_statistics)
      QV=$(grep -m1 "QV"       "${F}" | awk '{print $2}')
      MAP=$(grep -m1 "Mapping" "${F}" | awk '{print $4}')
      DEP=$(grep -m1 "Depth"   "${F}" | awk '{print $2}')
      BAS=$(grep -m1 "Base"    "${F}" | awk '{print $3}')
      EXP=$(grep "Small-scale" "${F}" | head -n 2 | tail -n 1 | awk '{print $3}')
      COL=$(grep "Small-scale" "${F}" | head -n 3 | tail -n 1 | awk '{print $3}')
      echo -e "${ID}\t${QV}\t${MAP}\t${DEP}\t${BAS}\t${EXP}\t${COL}"
    done
  } > "${STAT_OUT}/inspector/inspector_summary.tsv"

  # 4-2) QUAST-LG v5.3 against GRCh38 and CHM13
  for FA in "${POLISH_OUT}"/corrected/*.polished.fa; do
    ID=$(basename "${FA}" .polished.fa)

    quast-lg.py "${FA}" \
      -r "${GRCH38}" -g "${GRCH38_GFF}" \
      -o "${STAT_OUT}/quast_grch38/${ID}" \
      --report-all-metrics -e --fragmented --large \
      --est-ref-size 3100000000 --no-icarus -t "${THREADS}"

    quast-lg.py "${FA}" \
      -r "${CHM13}" \
      -o "${STAT_OUT}/quast_chm13/${ID}" \
      --report-all-metrics -e --fragmented --large \
      --est-ref-size 3100000000 --no-icarus -t "${THREADS}"
  done

  # 4-3) NGx curves extracted from the QUAST HTML report (values converted to Mb)
  for D in "${STAT_OUT}"/quast_chm13/*/; do
    ID=$(basename "${D}")
    [[ -s "${D}/report.html" ]] || continue

    tr '\n' ' ' < "${D}/report.html" \
      | grep -oP "(?<=<div id=['\"]coord-ngx-json['\"]>).*?(?=</div>)" \
      | jq '.coord_x' | tr -d '[]\n' | tr ',' '\n' | sed '/^$/d' > "${D}/x.txt"

    tr '\n' ' ' < "${D}/report.html" \
      | grep -oP "(?<=<div id=['\"]coord-ngx-json['\"]>).*?(?=</div>)" \
      | jq '.coord_y' | tr -d '[]\n' | tr ',' '\n' | sed '/^$/d' \
      | awk '{printf "%.6f\n", $1/1000000}' > "${D}/y.txt"

    paste "${D}/x.txt" "${D}/y.txt" > "${STAT_OUT}/ngx/${ID}_ngx.txt"
  done

  # 4-4) assemblies meeting N50 >= 20 Mb and <= 1,500 contigs
  {
    echo -e "assembly\tnum_contigs\tN50\tstatus"
    for D in "${STAT_OUT}"/quast_chm13/*/; do
      ID=$(basename "${D}")
      [[ -s "${D}/report.tsv" ]] || continue
      NC=$(awk -F'\t' '$1=="# contigs"{print $2; exit}' "${D}/report.tsv")
      N50=$(awk -F'\t' '$1=="N50"{print $2; exit}'      "${D}/report.tsv")
      ST=$(awk -v n="${NC}" -v s="${N50}" -v mn="${MIN_N50}" -v mc="${MAX_CONTIG}" \
             'BEGIN{print (s>=mn && n<=mc) ? "PASS" : "FAIL"}')
      echo -e "${ID}\t${NC}\t${N50}\t${ST}"
    done
  } > "${STAT_OUT}/assembly_qc_summary.tsv"

  awk -F'\t' '$4=="PASS"{split($1,a,"_"); print a[1]}' \
    "${STAT_OUT}/assembly_qc_summary.tsv" | sort -u > "${QC_PASS_LIST}"
  log "  assemblies passing quality control: $(wc -l < "${QC_PASS_LIST}") samples"

  # 4-5) interchromosomal misjoins, Minigraph v0.21
  #      contigs aligning to more than one chromosome over at least 100 kb
  for FA in "${POLISH_OUT}"/corrected/*.polished.fa; do
    ID=$(basename "${FA}" .polished.fa)

    minigraph -xasm -K1.9g --show-unmap=yes -t "${THREADS}" \
      "${CHM13}" "${FA}" > "${STAT_OUT}/misjoin/${ID}.paf"

    awk -F'\t' '$6!="*" && ($9-$8)>=100000 {print $1"\t"$6}' \
      "${STAT_OUT}/misjoin/${ID}.paf" \
      | sort -u | cut -f1 | uniq -c \
      | awk -v id="${ID}" '$1>1 {print id"\t"$2"\t"$1}' \
      >> "${STAT_OUT}/misjoin/interchromosomal_misjoins.tsv"
  done

  # 4-6) misassembly flagging, HMM-Flagger v1.0.0
  #      long reads are mapped back to the diploid assembly and each block is
  #      labelled as erroneous, duplicated, haploid or collapsed
  while read -r S; do
    [[ -z "${S}" ]] && continue
    mkdir -p "${STAT_OUT}/flagger/${S}"

    cat "${POLISH_OUT}/corrected/${S}_1.polished.fa" \
        "${POLISH_OUT}/corrected/${S}_2.polished.fa" \
        > "${STAT_OUT}/flagger/${S}/${S}.dip.fa"
    samtools faidx "${STAT_OUT}/flagger/${S}/${S}.dip.fa"

    minimap2 -ax map-hifi -t "${THREADS}" \
      "${STAT_OUT}/flagger/${S}/${S}.dip.fa" \
      "${QC_OUT}/trimmed/${S}.hifi.filt.fastq.gz" \
      | samtools sort -@ 8 -o "${STAT_OUT}/flagger/${S}/${S}.dip.bam" -
    samtools index -@ 8 "${STAT_OUT}/flagger/${S}/${S}.dip.bam"

    bam2cov --bam "${STAT_OUT}/flagger/${S}/${S}.dip.bam" \
            --output "${STAT_OUT}/flagger/${S}/${S}.cov.gz" \
            --threads "${THREADS}"

    hmm_flagger --input "${STAT_OUT}/flagger/${S}/${S}.cov.gz" \
                --outputDir "${STAT_OUT}/flagger/${S}/out" \
                --labelNames Err,Dup,Hap,Col \
                --threads "${THREADS}"
  done < "${QC_PASS_LIST}"

  log "step 4 finished"
}

###############################################################################
# Step 5. Repeat and gene annotation
###############################################################################

annotation() {
  log "step 5: repeat and gene annotation"
  mkdir -p "${ANNOT_OUT}"/{repeatmasker_human,repeatmasker_satellite,vntr,lowcomplexity,liftoff}

  for FA in "${POLISH_OUT}"/corrected/*.polished.fa; do
    ID=$(basename "${FA}" .polished.fa)

    # 5-1) RepeatMasker v4.1.6, human repeat library
    RepeatMasker -species human -pa "${THREADS}" -e ncbi -a -xsmall \
      -dir "${ANNOT_OUT}/repeatmasker_human" "${FA}"

    # 5-2) RepeatMasker v4.1.6, augmented CHM13 satellite library
    RepeatMasker -nolow -s -xsmall -e ncbi -pa "${THREADS}" \
      -lib "${SAT_LIB}" \
      -dir "${ANNOT_OUT}/repeatmasker_satellite" \
      "${ANNOT_OUT}/repeatmasker_human/$(basename "${FA}").masked"

    MASKED="${ANNOT_OUT}/repeatmasker_satellite/$(basename "${FA}").masked.masked"

    # 5-3) VNTRs with motifs of at least 7 bp
    etrf -m 7 "${MASKED}" > "${ANNOT_OUT}/vntr/${ID}.vntr.bed"

    # 5-4) low-complexity regions
    sdust "${MASKED}" > "${ANNOT_OUT}/lowcomplexity/${ID}.sdust.bed"

    # 5-5) gene annotation and gene copy number, Liftoff v1.6.3 with GENCODE v38
    liftoff -g "${GENCODE_GTF}" \
            -o "${ANNOT_OUT}/liftoff/${ID}.gencode38.gff3" \
            -u "${ANNOT_OUT}/liftoff/${ID}.unmapped.txt" \
            -copies -sc 0.95 -p "${THREADS}" \
            "${FA}" "${GRCH38}"

    # 5-6) protein-coding genes present in more than one copy
    awk -F'\t' '$3=="gene"' "${ANNOT_OUT}/liftoff/${ID}.gencode38.gff3" \
      | grep 'gene_type=protein_coding' \
      | sed -n 's/.*gene_name=\([^;]*\).*/\1/p' \
      | sort | uniq -c \
      | awk -v id="${ID}" '$1>1 {print id"\t"$2"\t"$1}' \
      >> "${ANNOT_OUT}/liftoff/duplicated_protein_coding_genes.tsv"
  done

  log "step 5 finished"
}

###############################################################################
# Step 6. Structural variant detection
###############################################################################

sv_calling() {
  log "step 6: structural variant detection"
  mkdir -p "${SV_OUT}"/{bam,pbsv,sniffles,svim,svimasm,merged}

  while read -r S; do
    [[ -z "${S}" ]] && continue

    # 6-1) read alignment, pbmm2 v1.17
    pbmm2 align "${GRCH38}" "${QC_OUT}/trimmed/${S}.hifi.filt.fastq.gz" \
                "${SV_OUT}/bam/${S}.GRCh38.bam" \
                --preset HIFI --sort --rg "@RG\tID:${S}\tSM:${S}" \
                -j "${THREADS}"
    samtools index -@ 8 "${SV_OUT}/bam/${S}.GRCh38.bam"

    # 6-2) pbsv v2.9.0
    pbsv discover "${SV_OUT}/bam/${S}.GRCh38.bam" "${SV_OUT}/pbsv/${S}.svsig.gz"
    pbsv call -j "${THREADS}" "${GRCH38}" \
              "${SV_OUT}/pbsv/${S}.svsig.gz" "${SV_OUT}/pbsv/${S}.pbsv.vcf"

    # 6-3) Sniffles v2.2, at least 10 supporting reads
    sniffles --input "${SV_OUT}/bam/${S}.GRCh38.bam" \
             --reference "${GRCH38}" \
             --vcf "${SV_OUT}/sniffles/${S}.sniffles.vcf" \
             --minsupport 10 --threads "${THREADS}"

    # 6-4) SVIM v2.0.0, at least 10 supporting reads
    svim alignment "${SV_OUT}/svim/${S}" \
         "${SV_OUT}/bam/${S}.GRCh38.bam" "${GRCH38}" \
         --minimum_depth 10
    cp "${SV_OUT}/svim/${S}/variants.vcf" "${SV_OUT}/svim/${S}.svim.vcf"

    # 6-5) SVIM-asm v1.0.3, assembly based calls from both haplotypes
    for H in 1 2; do
      minimap2 -a -x asm5 --cs -r2k -t "${THREADS}" \
        "${GRCH38}" "${POLISH_OUT}/corrected/${S}_${H}.polished.fa" \
        | samtools sort -@ 8 -o "${SV_OUT}/svimasm/${S}_${H}.bam" -
      samtools index -@ 8 "${SV_OUT}/svimasm/${S}_${H}.bam"
    done
    svim-asm diploid "${SV_OUT}/svimasm/${S}" \
             "${SV_OUT}/svimasm/${S}_1.bam" "${SV_OUT}/svimasm/${S}_2.bam" \
             "${GRCH38}"
    cp "${SV_OUT}/svimasm/${S}/variants.vcf" "${SV_OUT}/svimasm/${S}.svimasm.vcf"

    # 6-6) per-sample merge, 100 bp breakpoint distance, SVs >= 51 bp
    ls "${SV_OUT}/pbsv/${S}.pbsv.vcf" \
       "${SV_OUT}/sniffles/${S}.sniffles.vcf" \
       "${SV_OUT}/svim/${S}.svim.vcf" \
       "${SV_OUT}/svimasm/${S}.svimasm.vcf" > "${SV_OUT}/merged/${S}.filelist"

    SURVIVOR merge "${SV_OUT}/merged/${S}.filelist" 100 1 1 0 0 51 \
      "${SV_OUT}/merged/${S}.merged.vcf"
  done < "${QC_PASS_LIST}"

  # 6-7) cohort merge, SVs shared by every sample
  N_SAMPLE=$(wc -l < "${QC_PASS_LIST}")
  ls "${SV_OUT}"/merged/*.merged.vcf > "${SV_OUT}/merged/cohort.filelist"

  SURVIVOR merge "${SV_OUT}/merged/cohort.filelist" 100 "${N_SAMPLE}" 1 0 0 51 \
    "${SV_OUT}/merged/shared_SV.vcf"

  log "step 6 finished"
}

###############################################################################
# Entry point
###############################################################################

usage() {
  cat <<EOF
Usage: bash $(basename "$0") <step>

  read_qc       long-read and Hi-C read quality control
  assembly      HiFiasm, HiFiasm Hi-C mode and Flye assemblies
  polishing     Inspector evaluation and correction
  asm_stats     Inspector, QUAST-LG, NGx, misjoin and HMM-Flagger statistics
  annotation    RepeatMasker, etrf, sdust and Liftoff
  sv_calling    pbsv, Sniffles, SVIM, SVIM-asm and SURVIVOR merging
  all           run every step in order
EOF
}

main() {
  case "${1:-}" in
    read_qc)    read_qc ;;
    assembly)   assembly ;;
    polishing)  polishing ;;
    asm_stats)  asm_stats ;;
    annotation) annotation ;;
    sv_calling) sv_calling ;;
    all)
      read_qc
      assembly
      polishing
      asm_stats
      annotation
      sv_calling
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
