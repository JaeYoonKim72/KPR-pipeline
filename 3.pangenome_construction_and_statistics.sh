#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-16}
CACTUS_MEM=${CACTUS_MEM:-800Gi}       # --consMemory
CACTUS_DISK=${CACTUS_DISK:-1.3Ti}     # --maxDisk

REF_DIR=${REF_DIR:?set REF_DIR to the directory containing the reference files}
SEQFILE_DIR=${SEQFILE_DIR:?set SEQFILE_DIR to the directory containing the seqfiles}
ASM_DIR=${ASM_DIR:-}                  # haplotype assemblies, needed by allele_mapping
OUT_DIR=${OUT_DIR:-$(pwd)/results}
TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}
SCRIPT_DIR=${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
SRC_DIR="${SCRIPT_DIR}/src"

GRAPHS=${GRAPHS:-"kpr hprc_kpr"}

REFERENCE=${REFERENCE:-GRCh38}

FILTER=${FILTER:-9}

REF_PATHS=${REF_PATHS:-"GRCh38 CHM13"}

SV_MIN_LEN=${SV_MIN_LEN:-50}          # small variants versus structural variants
BUBBLE_LARGE=${BUBBLE_LARGE:-10000}   # large bubbles retained for annotation

ALLELE_MIN=${ALLELE_MIN:-5}           # minimum number of assemblies carrying an allele
ALLELE_MAX=${ALLELE_MAX:-29}          # bubbles carried by nearly every assembly

GRCH38=${GRCH38:-${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa}

CHM13=${CHM13:-${REF_DIR}/chm13v2.0.fa}

HPRC_GFA=${HPRC_GFA:-${REF_DIR}/hprc-v1.1-mc-grch38.gfa}

SAT_LIB=${SAT_LIB:-${REF_DIR}/chm13_satellite_library.fa}

GENE_BED=${GENE_BED:-${REF_DIR}/gencode.GRCh38.genes.name.bed}
EXON_BED=${EXON_BED:-${REF_DIR}/gencode.GRCh38.exons.name.bed}

MEDICAL_GENES=${MEDICAL_GENES:-${REF_DIR}/medically_relevant_genes.txt}

GRAPH_OUT=${OUT_DIR}/graph
INDEX_OUT=${OUT_DIR}/index
STAT_OUT=${OUT_DIR}/graph_statistics
VAR_OUT=${OUT_DIR}/graph_variants
GROWTH_OUT=${OUT_DIR}/growth
CALL_OUT=${OUT_DIR}/allele_mapping
BUBBLE_OUT=${OUT_DIR}/bubbles
SEQ_OUT=${OUT_DIR}/nonreference_sequence
LOG_OUT=${OUT_DIR}/logs

mkdir -p "${GRAPH_OUT}" "${INDEX_OUT}" "${STAT_OUT}" "${VAR_OUT}" \
         "${GROWTH_OUT}" "${CALL_OUT}" "${BUBBLE_OUT}" "${SEQ_OUT}" \
         "${LOG_OUT}" "${TMP_DIR}"

log() { echo -e "[$(date '+%F %T')] $*" | tee -a "${LOG_OUT}/pipeline.log"; }

# Step 1. Pangenome graph construction
construct() {
  log "step 1: pangenome graph construction"

  for G in ${GRAPHS}; do
    SEQFILE="${SEQFILE_DIR}/${G}.seqfile"
    [[ -s "${SEQFILE}" ]] || { log "  ${SEQFILE} not found, skipped"; continue; }

    mkdir -p "${GRAPH_OUT}/${G}"

    # shellcheck disable=SC2086
    cactus-pangenome "${TMP_DIR}/js_${G}" "${SEQFILE}" \
      --outName "${G}" --outDir "${GRAPH_OUT}/${G}" \
      --reference ${REFERENCE} \
      --filter "${FILTER}" \
      --batchSystem single_machine \
      --mgCores "${THREADS}" --mapCores "${THREADS}" \
      --consCores "${THREADS}" --indexCores "${THREADS}" \
      --consMemory "${CACTUS_MEM}" --maxDisk "${CACTUS_DISK}" \
      --maxLocalJobs 1000 \
      --giraffe clip filter \
      --chrom-vg clip filter --chrom-og clip filter \
      --gbz clip filter full \
      --gfa clip full \
      --vcf --viz --odgi \
      --workDir "${TMP_DIR}" \
      --logFile "${LOG_OUT}/${G}.cactus.log"
  done

  log "step 1 finished"
}

# Step 2. Graph indexing and format conversion
index_graph() {
  log "step 2: graph indexing"

  for G in ${GRAPHS}; do
    GBZ="${GRAPH_OUT}/${G}/${G}.gbz"
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GBZ}" ]] || { log "  ${GBZ} not found, skipped"; continue; }

    OUT="${INDEX_OUT}/${G}"
    mkdir -p "${OUT}"

    vg convert -t "${THREADS}" "${GBZ}" > "${OUT}/${G}.vg"
    vg index -x "${OUT}/${G}.xg" "${OUT}/${G}.vg" -b "${TMP_DIR}" -t "${THREADS}" -p

    vg gbwt -r "${OUT}/${G}.ri" -Z "${GBZ}"
    vg haplotypes -v 2 -H "${OUT}/${G}.hapl" \
                  -d "${GRAPH_OUT}/${G}/${G}.dist" \
                  -r "${OUT}/${G}.ri" "${GBZ}"

    vg gbwt -x "${OUT}/${G}.xg" -o "${OUT}/${G}.gbwt" -E

    vg paths -L -x "${OUT}/${G}.xg" > "${OUT}/${G}.paths.txt"
    cut -d'#' -f1   "${OUT}/${G}.paths.txt" | sort -u > "${OUT}/${G}.paths.samples"
    cut -d'#' -f1,2 "${OUT}/${G}.paths.txt" | sort -u > "${OUT}/${G}.paths.haplotypes"

    vg convert -W -t "${THREADS}" -g "${GFA}" -f > "${OUT}/${G}.v1.gfa"
    odgi build -g "${OUT}/${G}.v1.gfa" -o "${OUT}/${G}.v1.og" -s -O -t "${THREADS}"
    odgi sort -i "${OUT}/${G}.v1.og" -o "${OUT}/${G}.v1.sort.og" -O

    gfatools gfa2fa -l 60 "${GFA}" > "${OUT}/${G}.nodes.fa"
  done

  log "step 2 finished"
}

# Step 3. Graph statistics and pangenome size
graph_stats() {
  log "step 3: graph statistics"

  for G in ${GRAPHS}; do
    OG="${INDEX_OUT}/${G}/${G}.v1.sort.og"
    [[ -s "${OG}" ]] || { log "  ${OG} not found, skipped"; continue; }

    odgi stats -i "${OG}" -S > "${STAT_OUT}/${G}.odgi_stats.txt"

    odgi heaps -i "${OG}" -n 200 -t "${THREADS}" > "${STAT_OUT}/${G}.odgi_heaps.tsv"
  done

  log "step 3 finished"
}

# Step 4. Variant identification from the graph
count_alleles() {
  local VCF=$1
  bcftools view -H "${VCF}" \
    | awk -F'\t' -v f="$(basename "${VCF}")" '
        { total++; n = split($5, alts, ",")
          if (n == 1) a1++; else if (n == 2) a2++; else a3++ }
        END { print f, total, a1+0, a2+0, a3+0 }' \
    > "${VCF}.allele_count.txt"
}

variants() {
  log "step 4: variant identification"

  for G in ${GRAPHS}; do
    XG="${INDEX_OUT}/${G}/${G}.xg"
    [[ -s "${XG}" ]] || { log "  ${XG} not found, skipped"; continue; }

    OUT="${VAR_OUT}/${G}"
    mkdir -p "${OUT}"/{reference,samples,haplotypes}

    for P in ${REF_PATHS}; do
      vg deconstruct -P "${P}" -a -t "${THREADS}" "${XG}" \
        | bgzip -c > "${OUT}/reference/${G}.${P}.vcf.gz"
      bcftools index -f -t "${OUT}/reference/${G}.${P}.vcf.gz"
    done

    for P in ${REF_PATHS}; do
      V="${OUT}/reference/${G}.${P}"

      bcftools view -i 'strlen(REF)==1 && strlen(ALT)==1 && TYPE="snp"' -Oz "${V}.vcf.gz" \
        | bcftools norm --rm-dup all -Oz -o "${V}.SNP.vcf.gz"
      bcftools view -i 'TYPE="mnp"' -Oz "${V}.vcf.gz" \
        | bcftools norm --rm-dup all -Oz -o "${V}.MNP.vcf.gz"
      bcftools view -i 'TYPE="indel"' -Oz "${V}.vcf.gz" \
        | bcftools norm --rm-dup all -Oz -o "${V}.INDEL.vcf.gz"

      for T in vcf SNP.vcf MNP.vcf INDEL.vcf; do
        count_alleles "${V}.${T}.gz"
      done
    done

    while read -r P; do
      [[ -z "${P}" ]] && continue
      vg deconstruct -P "${P}" -a -t "${THREADS}" "${XG}" \
        | bgzip -c > "${OUT}/samples/${G}.${P}.vcf.gz"
      bcftools index -f -t "${OUT}/samples/${G}.${P}.vcf.gz"
    done < "${INDEX_OUT}/${G}/${G}.paths.samples"

    while read -r P; do
      [[ -z "${P}" ]] && continue
      vg deconstruct -P "${P}" -a -t "${THREADS}" "${XG}" \
        | bgzip -c > "${OUT}/haplotypes/${G}.${P}.vcf.gz"
      bcftools index -f -t "${OUT}/haplotypes/${G}.${P}.vcf.gz"
    done < "${INDEX_OUT}/${G}/${G}.paths.haplotypes"

    for V in "${OUT}"/reference/*.vcf.gz "${OUT}"/samples/*.vcf.gz "${OUT}"/haplotypes/*.vcf.gz; do
      case "${V}" in *.SNP.vcf.gz|*.MNP.vcf.gz|*.INDEL.vcf.gz|*var_sizes*) continue ;; esac

      bash "${SRC_DIR}/decon_split.sh" "${V}" "${SV_MIN_LEN}"

      P="${V%.vcf.gz}"
      for F in \
        "var_sizes.Smalls_list.no_mixed.SNP.SNP" \
        "var_sizes.Smalls_list.no_mixed.INS.INS" \
        "var_sizes.Smalls_list.no_mixed.DEL.DEL" \
        "var_sizes.SVs_list.all_50_alts.no_mixed.INS.INS" \
        "var_sizes.SVs_list.all_50_alts.no_mixed.DEL.DEL" ; do
        python3 "${SRC_DIR}/typecount.py" "${P}.${F}.vcf.gz"
      done
    done

    cat "${OUT}"/reference/*.allele_count.txt > "${STAT_OUT}/${G}.allele_counts.txt"
    { head -n 1 "$(find "${OUT}" -name '*.typecount' | head -n 1)"
      find "${OUT}" -name '*.typecount' -exec tail -n +2 {} \; ; } \
      > "${STAT_OUT}/${G}.variant_typecounts.tsv"
  done

  log "step 4 finished"
}

# Step 5. Pangenome growth curves and node frequency classes
growth() {
  log "step 5: pangenome growth curves"

  run_panacus() {                    # $1 = label, $2 = GFA
    local LABEL=$1 GFA=$2
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; return; }

    panacus histgrowth -t "${THREADS}" -c node -l 1,2,1,1 -q 0,0,0.05,0.95 -S \
      "${GFA}" > "${GROWTH_OUT}/${LABEL}.growth.node.tsv"
    panacus histgrowth -t "${THREADS}" -c bp -l 1,2,1,1 -q 0,0,0.05,0.95 -S \
      "${GFA}" > "${GROWTH_OUT}/${LABEL}.growth.bp.tsv"

    panacus hist -t "${THREADS}" -c node "${GFA}" > "${GROWTH_OUT}/${LABEL}.hist.tsv"

    awk -F'\t' 'BEGIN{OFS="\t"}
      /^#/ {next}
      NF>=2 && $1+0==$1 {cov[$1]=$2; if($1>max) max=$1}
      END{
        for (c in cov) {
          f = c / max
          if      (c == 1)    s="singleton"
          else if (c == 2)    s="doubleton"
          else if (f > 0.95)  s="core"
          else if (f >= 0.05) s="common"
          else                s="rare"
          n[s] += cov[c]
        }
        print "class", "nodes"
        for (s in n) print s, n[s]
      }' "${GROWTH_OUT}/${LABEL}.hist.tsv" \
      > "${GROWTH_OUT}/${LABEL}.node_classes.tsv"
  }

  for G in ${GRAPHS}; do
    run_panacus "${G}" "${INDEX_OUT}/${G}/${G}.v1.gfa"
  done

  run_panacus "hprc_v1.1" "${HPRC_GFA}"

  log "step 5 finished"
}

# Step 6. Allele mapping of the assemblies onto the graph
allele_mapping() {
  log "step 6: allele mapping"

  [[ -n "${ASM_DIR}" ]] || { log "  ASM_DIR is not set, skipped"; return; }

  for G in ${GRAPHS}; do
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; continue; }

    OUT="${CALL_OUT}/${G}"
    mkdir -p "${OUT}/call"

    for FA in "${ASM_DIR}"/*.fa; do
      B=$(basename "${FA}" .fa)
      minigraph -cxasm --call -t "${THREADS}" "${GFA}" "${FA}" \
        1> "${OUT}/call/${G}.${B}.call" \
        2> "${LOG_OUT}/${G}.${B}.call.err"
    done

    paste "${OUT}"/call/*.call \
      | awk '{
          numFiles = NF / 6
          merged = $6
          for (i = 2; i <= numFiles; i++) merged = merged "_" $(i*6)
          count = 0
          n = split(merged, arr, "_")
          for (i = 1; i <= n; i++) if (arr[i] != ".") count++
          print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" count "\t" merged
        }' > "${OUT}/${G}.minigraph.merged.call"
  done

  log "step 6 finished"
}

# Step 7. Bubbles and non-reference nodes
nonref_nodes() {
  log "step 7: bubbles and non-reference nodes"

  for G in ${GRAPHS}; do
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; continue; }

    OUT="${BUBBLE_OUT}/${G}"
    mkdir -p "${OUT}"

    gfatools bubble "${GFA}" > "${OUT}/${G}.bubble"

    awk '{ diff = ($3 - $2 < 0) ? -($3 - $2) : ($3 - $2); if (diff != 1) print }' \
      "${OUT}/${G}.bubble" > "${OUT}/${G}.bubble.noSNP.txt"

    awk -v m="${SV_MIN_LEN}" \
      '($13 != "*" || $14 != "*") && (length($13) > m || length($14) > m)' \
      "${OUT}/${G}.bubble.noSNP.txt" > "${OUT}/${G}.bubble.noSNP.SV.txt"

    awk -v m="${BUBBLE_LARGE}" \
      '($13 != "*" || $14 != "*") && (length($13) > m || length($14) > m)' \
      "${OUT}/${G}.bubble.noSNP.SV.txt" > "${OUT}/${G}.bubble.noSNP.SV.large.txt"

    cut -f -8 "${OUT}/${G}.bubble.noSNP.SV.large.txt" \
      > "${OUT}/${G}.bubble.noSNP.SV.large.cut8.txt"

    : > "${OUT}/${G}.reference_segments.txt"
    for P in ${REF_PATHS}; do
      awk -F'\t' -v ref="${P}" '$1=="W" && $2==ref {print $7}' "${GFA}" \
        | grep -o '[0-9]\+' >> "${OUT}/${G}.reference_segments.txt"
    done
    sort -u -o "${OUT}/${G}.reference_segments.txt" "${OUT}/${G}.reference_segments.txt"

    cut -f12 "${OUT}/${G}.bubble.noSNP.SV.txt" | tr ',' '\n' | sed '/^$/d' | sort -u \
      > "${OUT}/${G}.bubble_segments.txt"
    comm -23 "${OUT}/${G}.bubble_segments.txt" "${OUT}/${G}.reference_segments.txt" \
      > "${OUT}/${G}.nonreference_segments.txt"

    awk -F'\t' -v ids="${OUT}/${G}.nonreference_segments.txt" -v m="${SV_MIN_LEN}" '
      BEGIN{ while ((getline l < ids) > 0) keep[l]=1 }
      $1=="S" && ($2 in keep) && length($3) > m { print ">"$2"\n"$3 }' "${GFA}" \
      > "${OUT}/${G}.nonreference_SV_nodes.fa"

    MERGED="${CALL_OUT}/${G}/${G}.minigraph.merged.call"
    if [[ -s "${MERGED}" ]]; then
      N_ASM=$(awk 'NR==1{n=split($7,a,"_"); print n; exit}' "${MERGED}")

      awk -F'\t' -v n="${N_ASM}" 'BEGIN{OFS="\t"}
        { f = $6 / n
          if      ($6 == 1)   c = "singleton"
          else if ($6 == 2)   c = "doubleton"
          else if (f > 0.95)  c = "core"
          else if (f >= 0.05) c = "common"
          else                c = "rare"
          print $1, $2, $3, $6, c }' "${MERGED}" \
        > "${OUT}/${G}.bubble_frequency_class.tsv"

      cut -f5 "${OUT}/${G}.bubble_frequency_class.tsv" | sort | uniq -c \
        | awk '{print $2"\t"$1}' > "${STAT_OUT}/${G}.bubble_frequency_classes.tsv"
    else
      log "  ${MERGED} not found, frequency classes skipped"
    fi

    mkdir -p "${OUT}/repeatmasker_human" "${OUT}/repeatmasker_satellite"

    RepeatMasker -species human -pa "${THREADS}" -e ncbi -a -xsmall \
      -dir "${OUT}/repeatmasker_human" "${OUT}/${G}.nonreference_SV_nodes.fa"

    RepeatMasker -nolow -s -xsmall -e ncbi -pa "${THREADS}" \
      -lib "${SAT_LIB}" \
      -dir "${OUT}/repeatmasker_satellite" \
      "${OUT}/repeatmasker_human/${G}.nonreference_SV_nodes.fa.masked"
  done

  log "step 7 finished"
}

# Step 8. Annotation of the large bubbles
bubble_annotation() {
  log "step 8: bubble annotation"

  for G in ${GRAPHS}; do
    OUT="${BUBBLE_OUT}/${G}"
    BED="${OUT}/${G}.bubble.noSNP.SV.large.cut8.txt"
    [[ -s "${BED}" ]] || { log "  ${BED} not found, skipped"; continue; }

    bedtools intersect -a "${BED}" -b "${GENE_BED}" -wa -wb \
      > "${OUT}/${G}.large.annotGene.txt"
    bedtools intersect -a "${OUT}/${G}.large.annotGene.txt" -b "${EXON_BED}" -wa -wb \
      > "${OUT}/${G}.large.annotGene.annotExon.txt"

    awk -F"\t" 'FNR==NR {
        gsub(/^[ \t\r]+|[ \t\r]+$/, "", $1); genes[tolower($1)] = 1; next
      }
      {
        gsub(/^[ \t\r]+|[ \t\r]+$/, "", $12); gsub(/^[ \t\r]+|[ \t\r]+$/, "", $18)
        gene12 = tolower($12); gene18 = tolower($18)
        if (gene12 in genes)      match_gene = $12
        else if (gene18 in genes) match_gene = $18
        else                      match_gene = ""
        OFS = "\t"
        if (match_gene != "") { $22 = match_gene; $23 = "Medical" }
        else                  { $22 = "-";        $23 = "-" }
        print
      }' "${MEDICAL_GENES}" "${OUT}/${G}.large.annotGene.annotExon.txt" \
      > "${OUT}/${G}.large.annotGene.annotExon.Medical.txt"

    MERGED="${CALL_OUT}/${G}/${G}.minigraph.merged.call"
    [[ -s "${MERGED}" ]] || { log "  ${MERGED} not found, annotation stops here"; continue; }

    awk -F"\t" 'BEGIN{OFS="\t"}
      FNR==NR { key = $1 "\t" $2 "\t" $3; data1[key] = $6; data2[key] = $7; next }
      { key = $1 "\t" $2 "\t" $3
        if (key in data1) { $24 = data1[key]; $25 = data2[key] }
        print }' "${MERGED}" "${OUT}/${G}.large.annotGene.annotExon.Medical.txt" \
      > "${OUT}/${G}.large.annotated.txt"

    awk -F"\t" -v m="${ALLELE_MIN}" '($23=="Medical" && $24>m)' \
      "${OUT}/${G}.large.annotated.txt" > "${OUT}/${G}.target.medical.txt"

    awk -F"\t" -v m="${ALLELE_MAX}" 'BEGIN{OFS="\t"} ($24>=m){NF--; print}' \
      "${OUT}/${G}.large.annotated.txt" > "${OUT}/${G}.target.shared.txt"
  done

  log "step 8 finished"
}

# Step 9. Non-reference sequence carried by each sample
nonref_sequence() {
  log "step 9: non-reference sequence"

  for G in ${GRAPHS}; do
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; continue; }

    OUT="${SEQ_OUT}/${G}"
    mkdir -p "${OUT}"/{nodes,leninfo,nonref,classes}
    : > "${OUT}/${G}.sequence_per_sample.txt"
    : > "${OUT}/${G}.nonreference_sequence_per_sample.txt"

    awk -F'\t' '$1=="S"' "${GFA}" > "${OUT}/${G}.S_line"
    awk '{ if (NF >= 3) print $2, length($3) }' "${OUT}/${G}.S_line" \
      > "${OUT}/${G}.node_len.txt"

    python3 "${SRC_DIR}/wline_sample_nodes.py" "${GFA}" "${OUT}/nodes/${G}" \
      > "${OUT}/${G}.sample_node_counts.txt"

    for F in "${OUT}"/nodes/*.nodes.txt; do
      B=$(basename "${F}" .nodes.txt)
      awk 'NR==FNR {len[$1]=$2; next} {print $1, len[$1]}' \
        "${OUT}/${G}.node_len.txt" "${F}" > "${OUT}/leninfo/${B}.leninfo"
      awk -v s="${B}" '{sum+=$2} END {print s, sum+0}' \
        "${OUT}/leninfo/${B}.leninfo" >> "${OUT}/${G}.sequence_per_sample.txt"
    done

    : > "${OUT}/${G}.reference_nodes.txt"
    for P in ${REF_PATHS}; do
      R=$(find "${OUT}/leninfo" -name "${G}.${P}.leninfo" | head -n 1)
      [[ -n "${R}" ]] && cut -d' ' -f1 "${R}" >> "${OUT}/${G}.reference_nodes.txt"
    done
    sort -u -o "${OUT}/${G}.reference_nodes.txt" "${OUT}/${G}.reference_nodes.txt"

    for F in "${OUT}"/leninfo/*.leninfo; do
      B=$(basename "${F}" .leninfo)
      case " ${REF_PATHS} " in *" ${B#"${G}."} "*) continue ;; esac

      awk 'NR==FNR {exclude[$1]; next} !($1 in exclude) {print $1, $2}' \
        "${OUT}/${G}.reference_nodes.txt" "${F}" > "${OUT}/nonref/${B}.nonref.txt"
      awk -v s="${B}" '{sum+=$2} END {print s, sum+0}' \
        "${OUT}/nonref/${B}.nonref.txt" >> "${OUT}/${G}.nonreference_sequence_per_sample.txt"
    done

    awk '{ node_len[$1] = $2
           n = split(FILENAME, p, "/"); s = p[n]; sub(/\.nonref\.txt$/, "", s)
           samples[$1] = (samples[$1] ? samples[$1] ";" s : s)
           count[$1]++ }
         END { for (node in node_len)
                 print node "\t" node_len[node] "\t" count[node] "\t" samples[node] }' \
      "${OUT}"/nonref/*.nonref.txt > "${OUT}/${G}.nonreference_node_summary.txt"

    awk '{sum[$3] += $2} END {for (c in sum) print c, sum[c]}' \
      "${OUT}/${G}.nonreference_node_summary.txt" | sort -k1,1n \
      > "${STAT_OUT}/${G}.nonreference_sequence_by_sharing.txt"

    N_SAMPLE=$(find "${OUT}/nonref" -maxdepth 1 -type f -name '*.nonref.txt' | wc -l)
    (( N_SAMPLE > 0 )) || { printf 'No sample node lists for %s\n' "${G}" >&2; return 1; }

    awk -F'\t' -v n="${N_SAMPLE}" 'BEGIN{OFS="\t"}
      { f = $3 / n
        if      ($3 == 1)   c = "singleton"
        else if ($3 == 2)   c = "doubleton"
        else if (f > 0.95)  c = "core"
        else if (f >= 0.05) c = "common"
        else                c = "rare"
        print $1, $2, $3, c }' "${OUT}/${G}.nonreference_node_summary.txt" \
      > "${OUT}/${G}.nonreference_node_class.tsv"

    awk -F'\t' '{n[$4]++; bp[$4]+=$2}
      END {print "class\tnodes\tbp"; for (c in n) print c"\t"n[c]"\t"bp[c]}' \
      "${OUT}/${G}.nonreference_node_class.tsv" \
      > "${STAT_OUT}/${G}.nonreference_node_classes.tsv"

    for C in core common doubleton singleton; do
      awk -F'\t' -v c="${C}" 'NR==FNR {if ($4 == c) keys[$1]; next} $2 in keys' \
        "${OUT}/${G}.nonreference_node_class.tsv" "${OUT}/${G}.S_line" \
        > "${OUT}/classes/${G}.nonreference.${C}.S_line"

      gfatools gfa2fa "${OUT}/classes/${G}.nonreference.${C}.S_line" \
        > "${OUT}/classes/${G}.nonreference.${C}.fa"
    done

    awk -F'\t' 'NR==FNR {keys[$1]; next} $2 in keys' \
      "${OUT}/${G}.nonreference_node_summary.txt" "${OUT}/${G}.S_line" \
      > "${OUT}/${G}.nonreference.S_line"
    gfatools gfa2fa "${OUT}/${G}.nonreference.S_line" \
      > "${OUT}/${G}.nonreference.fa"
  done

  log "step 9 finished"
}

usage() {
  cat <<EOF
Usage: bash $(basename "$0") <step>

  construct           Minigraph-Cactus pangenome graph construction
  index_graph         vg and odgi indexing and format conversion
  graph_stats         odgi stats and pangenome size estimation
  variants            vg deconstruct, size and type split, classification
  growth              Panacus growth curves and node frequency classes
  allele_mapping      minigraph --call of each assembly against the graph
  nonref_nodes        bubbles, non-reference SV nodes and their repeat content
  bubble_annotation   gene, exon and medically relevant gene annotation
  nonref_sequence     non-reference sequence per sample and per frequency class
  all                 run every step in order
EOF
}

main() {
  case "${1:-}" in
    construct)         construct ;;
    index_graph)       index_graph ;;
    graph_stats)       graph_stats ;;
    variants)          variants ;;
    growth)            growth ;;
    allele_mapping)    allele_mapping ;;
    nonref_nodes)      nonref_nodes ;;
    bubble_annotation) bubble_annotation ;;
    nonref_sequence)   nonref_sequence ;;
    all)
      construct
      index_graph
      graph_stats
      variants
      growth
      allele_mapping
      nonref_nodes
      bubble_annotation
      nonref_sequence
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
