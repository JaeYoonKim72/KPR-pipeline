#!/usr/bin/env bash
#
# Constructing a Korean pangenome reveals distinct haplotypes at the amylase
# locus
#
# pangenome_construction_and_statistics.sh
#
# Usage: bash pangenome_construction_and_statistics.sh <step>
#        (see docs for the step list and for the required inputs)

set -euo pipefail

###############################################################################
# Configuration
###############################################################################

THREADS=${THREADS:-16}
CACTUS_MEM=${CACTUS_MEM:-800Gi}       # --consMemory
CACTUS_DISK=${CACTUS_DISK:-1.3Ti}     # --maxDisk

REF_DIR=${REF_DIR:?set REF_DIR to the directory containing the reference files}
SEQFILE_DIR=${SEQFILE_DIR:?set SEQFILE_DIR to the directory containing the seqfiles}
ASM_DIR=${ASM_DIR:-}                  # haplotype assemblies, needed by allele_mapping
OUT_DIR=${OUT_DIR:-$(pwd)/results}
TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}
SCRIPT_DIR=${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

# Graphs to build and analyse. Each name must have a matching seqfile at
# ${SEQFILE_DIR}/<name>.seqfile, prepared by the user.
#   kpr       : KOR haplotype assemblies only
#   hprc_kpr  : KOR haplotype assemblies together with HPRC assemblies
GRAPHS=${GRAPHS:-"kpr hprc_kpr"}

# Reference used as the backbone of the graph. cactus-pangenome accepts more
# than one name, the first being the primary reference.
REFERENCE=${REFERENCE:-GRCh38}

# Haplotype frequency filter passed to cactus-pangenome. Paths carrying an
# allele frequency below about 10% are removed with vg clip in the filtered
# graph, which is the graph used for short-read mapping.
FILTER=${FILTER:-9}

# Reference path names as they appear in the graph.
REF_PATHS=${REF_PATHS:-"GRCh38 CHM13"}

# Size thresholds.
SV_MIN_LEN=${SV_MIN_LEN:-50}          # small variants versus structural variants
BUBBLE_LARGE=${BUBBLE_LARGE:-10000}   # large bubbles retained for annotation

# Thresholds used when selecting annotated target bubbles.
ALLELE_MIN=${ALLELE_MIN:-5}           # minimum number of assemblies carrying an allele
ALLELE_MAX=${ALLELE_MAX:-29}          # bubbles carried by nearly every assembly

# --- reference files ---------------------------------------------------------
# GRCh38 primary assembly
#   https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
GRCH38=${GRCH38:-${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa}

# T2T-CHM13 v2.0
#   https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
CHM13=${CHM13:-${REF_DIR}/chm13v2.0.fa}

# HPRC v1.1 Minigraph-Cactus graph, reanalysed with the same parameters
#   https://data.humanpangenome.org/alignments
#   https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html
HPRC_GFA=${HPRC_GFA:-${REF_DIR}/hprc-v1.1-mc-grch38.gfa}

# Augmented CHM13 satellite library used for the repeat content of
# non-reference nodes
#   https://github.com/marbl/CHM13  (T2T-CHM13 repeat annotation resources)
SAT_LIB=${SAT_LIB:-${REF_DIR}/chm13_satellite_library.fa}

# BED files of gene and exon coordinates, with the gene name in column 4,
# derived from the GENCODE annotation
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
GENE_BED=${GENE_BED:-${REF_DIR}/gencode.GRCh38.genes.name.bed}
EXON_BED=${EXON_BED:-${REF_DIR}/gencode.GRCh38.exons.name.bed}

# List of medically relevant gene symbols, one per line
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

###############################################################################
# Step 1. Pangenome graph construction
###############################################################################

construct() {
  log "step 1: pangenome graph construction"

  # Minigraph-Cactus v2.9. The pipeline aligns the assemblies sequentially with
  # Minigraph v0.19 to build an SV graph, soft-masks centromeric and telomeric
  # regions with dna-brnn, remaps the assemblies to the SV graph, discards
  # soft-masked segments longer than 100 kb and alignments with MAPQ below 5,
  # splits the graph by chromosome, runs Cactus for base-level alignment,
  # converts the HAL output with hal2vg, removes unaligned paths longer than
  # 10 kb, normalises the graph with GFAffix and merges the chromosome graphs
  # into a whole-genome graph. These stages are Minigraph-Cactus defaults for
  # human data and are not overridden here. The order of the assemblies in the
  # seqfile determines the order in which they are aligned.
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

###############################################################################
# Step 2. Graph indexing and format conversion
###############################################################################

index_graph() {
  log "step 2: graph indexing"

  # vg v1.63.1, odgi v0.9.0, gfatools v0.5
  for G in ${GRAPHS}; do
    GBZ="${GRAPH_OUT}/${G}/${G}.gbz"
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GBZ}" ]] || { log "  ${GBZ} not found, skipped"; continue; }

    OUT="${INDEX_OUT}/${G}"
    mkdir -p "${OUT}"

    # 2-1) GBZ to vg and xg
    vg convert -t "${THREADS}" "${GBZ}" > "${OUT}/${G}.vg"
    vg index -x "${OUT}/${G}.xg" "${OUT}/${G}.vg" -b "${TMP_DIR}" -t "${THREADS}" -p

    # 2-2) r-index and haplotype index, used for haplotype sampling
    vg gbwt -r "${OUT}/${G}.ri" -Z "${GBZ}"
    vg haplotypes -v 2 -H "${OUT}/${G}.hapl" \
                  -d "${GRAPH_OUT}/${G}/${G}.dist" \
                  -r "${OUT}/${G}.ri" "${GBZ}"

    # 2-3) GBWT of the embedded paths
    vg gbwt -x "${OUT}/${G}.xg" -o "${OUT}/${G}.gbwt" -E

    # 2-4) path list, and the sample-level and haplotype-level path names
    #      derived from it (PanSN: SAMPLE#HAPLOTYPE#CONTIG)
    vg paths -L -x "${OUT}/${G}.xg" > "${OUT}/${G}.paths.txt"
    cut -d'#' -f1   "${OUT}/${G}.paths.txt" | sort -u > "${OUT}/${G}.paths.samples"
    cut -d'#' -f1,2 "${OUT}/${G}.paths.txt" | sort -u > "${OUT}/${G}.paths.haplotypes"

    # 2-5) GFA v1 and odgi graph
    vg convert -W -t "${THREADS}" -g "${GFA}" -f > "${OUT}/${G}.v1.gfa"
    odgi build -g "${OUT}/${G}.v1.gfa" -o "${OUT}/${G}.v1.og" -s -O -t "${THREADS}"
    odgi sort -i "${OUT}/${G}.v1.og" -o "${OUT}/${G}.v1.sort.og" -O

    # 2-6) node sequences in FASTA format
    gfatools gfa2fa -l 60 "${GFA}" > "${OUT}/${G}.nodes.fa"
  done

  log "step 2 finished"
}

###############################################################################
# Step 3. Graph statistics and pangenome size
###############################################################################

graph_stats() {
  log "step 3: graph statistics"

  for G in ${GRAPHS}; do
    OG="${INDEX_OUT}/${G}/${G}.v1.sort.og"
    [[ -s "${OG}" ]] || { log "  ${OG} not found, skipped"; continue; }

    # 3-1) nodes, edges, paths and steps, odgi v0.9.0
    odgi stats -i "${OG}" -S > "${STAT_OUT}/${G}.odgi_stats.txt"

    # 3-2) pangenome size and saturation curve over 200 permutations
    odgi heaps -i "${OG}" -n 200 -t "${THREADS}" > "${STAT_OUT}/${G}.odgi_heaps.tsv"
  done

  log "step 3 finished"
}

###############################################################################
# Step 4. Variant identification from the graph
###############################################################################

# Total number of variants and the number of sites with one, two or more
# alternative alleles.
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

    # 4-1) variant sites relative to each reference path, vg v1.63.1
    for P in ${REF_PATHS}; do
      vg deconstruct -P "${P}" -a -t "${THREADS}" "${XG}" \
        | bgzip -c > "${OUT}/reference/${G}.${P}.vcf.gz"
      bcftools index -f -t "${OUT}/reference/${G}.${P}.vcf.gz"
    done

    # 4-2) small variants of the reference-based call sets subclassified into
    #      SNPs, MNPs and indels with bcftools v1.16, duplicates removed
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

    # 4-3) variant sites relative to each sample and to each haplotype, used to
    #      count the variants contributed by a single assembly
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

    # 4-4) size and type split of every call set, then classification of each
    #      variant as singleton, doubleton, common or core
    for V in "${OUT}"/reference/*.vcf.gz "${OUT}"/samples/*.vcf.gz "${OUT}"/haplotypes/*.vcf.gz; do
      case "${V}" in *.SNP.vcf.gz|*.MNP.vcf.gz|*.INDEL.vcf.gz|*var_sizes*) continue ;; esac

      bash "${SCRIPT_DIR}/decon_split.sh" "${V}" "${SV_MIN_LEN}"

      P="${V%.vcf.gz}"
      for F in \
        "var_sizes.Smalls_list.no_mixed.SNP.SNP" \
        "var_sizes.Smalls_list.no_mixed.INS.INS" \
        "var_sizes.Smalls_list.no_mixed.DEL.DEL" \
        "var_sizes.SVs_list.all_50_alts.no_mixed.INS.INS" \
        "var_sizes.SVs_list.all_50_alts.no_mixed.DEL.DEL" ; do
        python3 "${SCRIPT_DIR}/typecount.py" "${P}.${F}.vcf.gz"
      done
    done

    # 4-5) collected counts
    cat "${OUT}"/reference/*.allele_count.txt > "${STAT_OUT}/${G}.allele_counts.txt"
    { head -n 1 "$(find "${OUT}" -name '*.typecount' | head -n 1)"
      find "${OUT}" -name '*.typecount' -exec tail -n +2 {} \; ; } \
      > "${STAT_OUT}/${G}.variant_typecounts.tsv"
  done

  log "step 4 finished"
}

###############################################################################
# Step 5. Pangenome growth curves and node frequency classes
###############################################################################

growth() {
  log "step 5: pangenome growth curves"

  # Panacus v0.4.1. The same parameters are applied to the HPRC v1.1 graph so
  # that the pangenomes are directly comparable.
  run_panacus() {                    # $1 = label, $2 = GFA
    local LABEL=$1 GFA=$2
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; return; }

    # 5-1) accumulation of novel sequence over the assemblies, in nodes and in
    #      base pairs
    panacus histgrowth -t "${THREADS}" -c node -l 1,2,1,1 -q 0,0,0.05,0.95 -S \
      "${GFA}" > "${GROWTH_OUT}/${LABEL}.growth.node.tsv"
    panacus histgrowth -t "${THREADS}" -c bp -l 1,2,1,1 -q 0,0,0.05,0.95 -S \
      "${GFA}" > "${GROWTH_OUT}/${LABEL}.growth.bp.tsv"

    # 5-2) node coverage histogram, used for the frequency classes
    panacus hist -t "${THREADS}" -c node "${GFA}" > "${GROWTH_OUT}/${LABEL}.hist.tsv"

    # 5-3) core (> 95% of the assemblies), common (5-95%), doubletons (two
    #      assemblies) and singletons (one assembly)
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

###############################################################################
# Step 6. Allele mapping of the assemblies onto the graph
###############################################################################

allele_mapping() {
  log "step 6: allele mapping"

  [[ -n "${ASM_DIR}" ]] || { log "  ASM_DIR is not set, skipped"; return; }

  for G in ${GRAPHS}; do
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; continue; }

    OUT="${CALL_OUT}/${G}"
    mkdir -p "${OUT}/call"

    # 6-1) allele carried by each haplotype assembly at every bubble of the
    #      graph, Minigraph v0.19
    for FA in "${ASM_DIR}"/*.fa; do
      B=$(basename "${FA}" .fa)
      minigraph -cxasm --call -t "${THREADS}" "${GFA}" "${FA}" \
        1> "${OUT}/call/${G}.${B}.call" \
        2> "${LOG_OUT}/${G}.${B}.call.err"
    done

    # 6-2) merge the per-assembly calls; the sixth column of each file holds
    #      the allele, and the number of assemblies carrying an allele is
    #      counted per bubble
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

###############################################################################
# Step 7. Bubbles and non-reference nodes
###############################################################################

nonref_nodes() {
  log "step 7: bubbles and non-reference nodes"

  for G in ${GRAPHS}; do
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; continue; }

    OUT="${BUBBLE_OUT}/${G}"
    mkdir -p "${OUT}"

    # 7-1) bubbles of the graph, gfatools v0.5
    gfatools bubble "${GFA}" > "${OUT}/${G}.bubble"

    # 7-2) bubbles that are not simple substitutions
    awk '{ diff = ($3 - $2 < 0) ? -($3 - $2) : ($3 - $2); if (diff != 1) print }' \
      "${OUT}/${G}.bubble" > "${OUT}/${G}.bubble.noSNP.txt"

    # 7-3) bubbles carrying an allele longer than the structural variant
    #      threshold, and the subset of large bubbles kept for annotation
    awk -v m="${SV_MIN_LEN}" \
      '($13 != "*" || $14 != "*") && (length($13) > m || length($14) > m)' \
      "${OUT}/${G}.bubble.noSNP.txt" > "${OUT}/${G}.bubble.noSNP.SV.txt"

    awk -v m="${BUBBLE_LARGE}" \
      '($13 != "*" || $14 != "*") && (length($13) > m || length($14) > m)' \
      "${OUT}/${G}.bubble.noSNP.SV.txt" > "${OUT}/${G}.bubble.noSNP.SV.large.txt"

    cut -f -8 "${OUT}/${G}.bubble.noSNP.SV.large.txt" \
      > "${OUT}/${G}.bubble.noSNP.SV.large.cut8.txt"

    # 7-4) segments carried by the reference walks; non-reference nodes are
    #      those absent from every reference
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

    # 7-5) sequences of the non-reference structural variant nodes
    awk -F'\t' -v ids="${OUT}/${G}.nonreference_segments.txt" -v m="${SV_MIN_LEN}" '
      BEGIN{ while ((getline l < ids) > 0) keep[l]=1 }
      $1=="S" && ($2 in keep) && length($3) > m { print ">"$2"\n"$3 }' "${GFA}" \
      > "${OUT}/${G}.nonreference_SV_nodes.fa"

    # 7-6) frequency of each bubble across the assemblies, taken from the
    #      merged allele mapping, and the resulting frequency classes
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

    # 7-7) repeat content of the non-reference structural variant nodes,
    #      RepeatMasker v4.1.6 with the human library followed by the augmented
    #      CHM13 satellite library, as for the assemblies
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

###############################################################################
# Step 8. Annotation of the large bubbles
###############################################################################

bubble_annotation() {
  log "step 8: bubble annotation"

  for G in ${GRAPHS}; do
    OUT="${BUBBLE_OUT}/${G}"
    BED="${OUT}/${G}.bubble.noSNP.SV.large.cut8.txt"
    [[ -s "${BED}" ]] || { log "  ${BED} not found, skipped"; continue; }

    # 8-1) genes and exons overlapping each bubble
    bedtools intersect -a "${BED}" -b "${GENE_BED}" -wa -wb \
      > "${OUT}/${G}.large.annotGene.txt"
    bedtools intersect -a "${OUT}/${G}.large.annotGene.txt" -b "${EXON_BED}" -wa -wb \
      > "${OUT}/${G}.large.annotGene.annotExon.txt"

    # 8-2) flag the bubbles overlapping a medically relevant gene
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

    # 8-3) attach the number of assemblies carrying an allele at each bubble
    MERGED="${CALL_OUT}/${G}/${G}.minigraph.merged.call"
    [[ -s "${MERGED}" ]] || { log "  ${MERGED} not found, annotation stops here"; continue; }

    awk -F"\t" 'BEGIN{OFS="\t"}
      FNR==NR { key = $1 "\t" $2 "\t" $3; data1[key] = $6; data2[key] = $7; next }
      { key = $1 "\t" $2 "\t" $3
        if (key in data1) { $24 = data1[key]; $25 = data2[key] }
        print }' "${MERGED}" "${OUT}/${G}.large.annotGene.annotExon.Medical.txt" \
      > "${OUT}/${G}.large.annotated.txt"

    # 8-4) bubbles in medically relevant genes carried by several assemblies,
    #      and bubbles carried by nearly every assembly
    awk -F"\t" -v m="${ALLELE_MIN}" '($23=="Medical" && $24>m)' \
      "${OUT}/${G}.large.annotated.txt" > "${OUT}/${G}.target.medical.txt"

    awk -F"\t" -v m="${ALLELE_MAX}" 'BEGIN{OFS="\t"} ($24>=m){NF--; print}' \
      "${OUT}/${G}.large.annotated.txt" > "${OUT}/${G}.target.shared.txt"
  done

  log "step 8 finished"
}

###############################################################################
# Step 9. Non-reference sequence carried by each sample
###############################################################################

nonref_sequence() {
  log "step 9: non-reference sequence"

  for G in ${GRAPHS}; do
    GFA="${GRAPH_OUT}/${G}/${G}.gfa"
    [[ -s "${GFA}" ]] || { log "  ${GFA} not found, skipped"; continue; }

    OUT="${SEQ_OUT}/${G}"
    mkdir -p "${OUT}"/{nodes,leninfo,nonref,classes}

    # 9-1) segment lines and their lengths
    awk -F'\t' '$1=="S"' "${GFA}" > "${OUT}/${G}.S_line"
    awk '{ if (NF >= 3) print $2, length($3) }' "${OUT}/${G}.S_line" \
      > "${OUT}/${G}.node_len.txt"

    # 9-2) nodes traversed by each sample, taken from the W lines
    python3 "${SCRIPT_DIR}/wline_sample_nodes.py" "${GFA}" "${OUT}/nodes/${G}" \
      > "${OUT}/${G}.sample_node_counts.txt"

    # 9-3) length of every node carried by each sample, and the total sequence
    #      each sample contributes to the graph
    for F in "${OUT}"/nodes/*.nodes.txt; do
      B=$(basename "${F}" .nodes.txt)
      awk 'NR==FNR {len[$1]=$2; next} {print $1, len[$1]}' \
        "${OUT}/${G}.node_len.txt" "${F}" > "${OUT}/leninfo/${B}.leninfo"
      awk -v s="${B}" '{sum+=$2} END {print s, sum+0}' \
        "${OUT}/leninfo/${B}.leninfo" >> "${OUT}/${G}.sequence_per_sample.txt"
    done

    # 9-4) nodes absent from the reference samples, and the non-reference
    #      sequence each sample carries
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

    # 9-5) every non-reference node with its length, the number of samples
    #      carrying it and their names
    awk '{ node_len[$1] = $2
           n = split(FILENAME, p, "/"); s = p[n]; sub(/\.nonref\.txt$/, "", s)
           samples[$1] = (samples[$1] ? samples[$1] ";" s : s)
           count[$1]++ }
         END { for (node in node_len)
                 print node "\t" node_len[node] "\t" count[node] "\t" samples[node] }' \
      "${OUT}"/nonref/*.nonref.txt > "${OUT}/${G}.nonreference_node_summary.txt"

    # 9-6) total non-reference sequence at each level of sharing
    awk '{sum[$3] += $2} END {for (c in sum) print c, sum[c]}' \
      "${OUT}/${G}.nonreference_node_summary.txt" | sort -k1,1n \
      > "${STAT_OUT}/${G}.nonreference_sequence_by_sharing.txt"

    # 9-7) frequency classes of the non-reference nodes, and the sequence of
    #      each class in FASTA format
    N_SAMPLE=$(awk -F'\t' '{if ($3 > max) max = $3} END {print max}' \
      "${OUT}/${G}.nonreference_node_summary.txt")

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

    # all non-reference nodes together
    awk -F'\t' 'NR==FNR {keys[$1]; next} $2 in keys' \
      "${OUT}/${G}.nonreference_node_summary.txt" "${OUT}/${G}.S_line" \
      > "${OUT}/${G}.nonreference.S_line"
    gfatools gfa2fa "${OUT}/${G}.nonreference.S_line" \
      > "${OUT}/${G}.nonreference.fa"
  done

  log "step 9 finished"
}

###############################################################################
# Entry point
###############################################################################

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
