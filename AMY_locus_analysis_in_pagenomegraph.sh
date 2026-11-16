#!/usr/bin/env bash
#
# Constructing a Korean pangenome reveals distinct haplotypes at the amylase
# locus
#
# AMY_locus_analysis_in_pagenomegraph.sh
#
# Usage: bash AMY_locus_analysis_in_pagenomegraph.sh <step>
#        (see docs for the step list and for the required inputs)

set -euo pipefail

###############################################################################
# Configuration
###############################################################################

THREADS=${THREADS:-8}

GRAPH_DIR=${GRAPH_DIR:?set GRAPH_DIR to the directory containing the pangenome graph}
REF_DIR=${REF_DIR:?set REF_DIR to the directory containing the reference files}
BUBBLE_DIR=${BUBBLE_DIR:-}            # bubble tables from the construction pipeline
OUT_DIR=${OUT_DIR:-$(pwd)/results}
LIST_DIR=${LIST_DIR:-$(pwd)/lists}
TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}
SCRIPT_DIR=${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

# The analysis is carried out on the combined HPRC+KPR graph.
GRAPH=${GRAPH:-hprc_kpr}

# Bubble selection: bubbles longer than this and carrying at least this many
# per-sample alleles are retained as candidate regions.
BUBBLE_MIN_LEN=${BUBBLE_MIN_LEN:-10000}
ALLELE_MIN=${ALLELE_MIN:-5}

# Genes defining the AMY cluster (structural / gene-tag analysis only).
AMY_CLUSTER=${AMY_CLUSTER:-"AMY1A AMY1B AMY1C AMYP1 AMY2B AMY2A"}

# Region of the amylase cluster in the coordinates of the reference embedded in
# the graph, used to extract the subgraph.
AMY_REGION=${AMY_REGION:-chr1:103,570,000-103,760,000}
REF_SAMPLE=${REF_SAMPLE:-GRCh38}

# --- reference files ---------------------------------------------------------
# BED file of protein-coding genes with the gene name in column 4, derived from
# the GENCODE annotation
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
GENE_BED=${GENE_BED:-${REF_DIR}/gencode.GRCh38.genes.name.bed}

# --- sample tables -----------------------------------------------------------
# Two columns, sample name and population group
POP_MAP=${POP_MAP:-${LIST_DIR}/sample_population.tsv}

BUBBLE_OUT=${OUT_DIR}/bubble_selection
SUB_OUT=${OUT_DIR}/subgraph
HAP_OUT=${OUT_DIR}/haplotypes
STAT_OUT=${OUT_DIR}/haplotype_distribution
LOG_OUT=${OUT_DIR}/logs

mkdir -p "${BUBBLE_OUT}" "${SUB_OUT}" "${HAP_OUT}" "${STAT_OUT}" "${LOG_OUT}" "${TMP_DIR}" "${LIST_DIR}"

log() { echo -e "[$(date '+%F %T')] $*" | tee -a "${LOG_OUT}/pipeline.log"; }

###############################################################################
# Step 1. Selection of the regions of interest
###############################################################################

select_regions() {
  log "step 1: bubble selection"

  GFA="${GRAPH_DIR}/${GRAPH}/${GRAPH}.gfa"
  [[ -s "${GFA}" ]] || { log "  ${GFA} not found"; return 1; }

  # 1-1) bubbles of the graph, gfatools v0.5
  if [[ -n "${BUBBLE_DIR}" && -s "${BUBBLE_DIR}/${GRAPH}.bubble" ]]; then
    cp "${BUBBLE_DIR}/${GRAPH}.bubble" "${BUBBLE_OUT}/${GRAPH}.bubble"
  else
    gfatools bubble "${GFA}" > "${BUBBLE_OUT}/${GRAPH}.bubble"
  fi

  # 1-2) per-sample alleles at every bubble, Minigraph v0.21 (--cxasm --call),
  #      merged so that the number of assemblies carrying an allele is known
  MERGED="${BUBBLE_OUT}/${GRAPH}.minigraph.merged.call"
  if [[ -n "${BUBBLE_DIR}" && -s "${BUBBLE_DIR}/${GRAPH}.minigraph.merged.call" ]]; then
    cp "${BUBBLE_DIR}/${GRAPH}.minigraph.merged.call" "${MERGED}"
  else
    log "  ${MERGED} not found; run the allele_mapping step of the construction pipeline"
    return 1
  fi

  # 1-3) bubbles longer than the threshold and carrying enough alleles
  awk -F'\t' -v m="${BUBBLE_MIN_LEN}" 'BEGIN{OFS="\t"}
    { span = ($3 - $2 < 0) ? $2 - $3 : $3 - $2; if (span > m) print }' \
    "${BUBBLE_OUT}/${GRAPH}.bubble" \
    | cut -f -8 > "${BUBBLE_OUT}/${GRAPH}.bubble.large.bed"

  awk -F'\t' -v a="${ALLELE_MIN}" 'BEGIN{OFS="\t"}
    FNR==NR { if ($6 >= a) keep[$1"\t"$2"\t"$3] = $6; next }
    { key = $1"\t"$2"\t"$3; if (key in keep) print $0, keep[key] }' \
    "${MERGED}" "${BUBBLE_OUT}/${GRAPH}.bubble.large.bed" \
    > "${BUBBLE_OUT}/${GRAPH}.bubble.selected.bed"

  # 1-4) protein-coding genes overlapping the selected bubbles
  bedtools intersect -a "${BUBBLE_OUT}/${GRAPH}.bubble.selected.bed" \
                     -b "${GENE_BED}" -wa -wb \
    > "${BUBBLE_OUT}/${GRAPH}.bubble.selected.genes.txt"

  # 1-5) the bubbles overlapping the amylase cluster
  GREP_ARGS=()
  for GENE in ${AMY_CLUSTER}; do GREP_ARGS+=(-e "${GENE}"); done

  grep -wF "${GREP_ARGS[@]}" "${BUBBLE_OUT}/${GRAPH}.bubble.selected.genes.txt" \
    > "${BUBBLE_OUT}/${GRAPH}.bubble.AMY.txt" || true

  log "  selected bubbles: $(wc -l < "${BUBBLE_OUT}/${GRAPH}.bubble.selected.bed")"
  log "  overlapping the amylase cluster: $(wc -l < "${BUBBLE_OUT}/${GRAPH}.bubble.AMY.txt")"
  log "step 1 finished"
}

###############################################################################
# Step 2. Extraction and visualisation of the subgraph
###############################################################################

extract_subgraph() {
  log "step 2: subgraph extraction"

  GFA="${GRAPH_DIR}/${GRAPH}/${GRAPH}.gfa"
  [[ -s "${GFA}" ]] || { log "  ${GFA} not found"; return 1; }

  # 2-1) index the whole graph and cut out the region, gfabase v0.5
  gfabase load "${GFA}" "${SUB_OUT}/${GRAPH}.gfab" --view

  gfabase sub "${SUB_OUT}/${GRAPH}.gfab" \
    --range "${REF_SAMPLE}#0#${AMY_REGION}" \
    --view --cutpoints 1 \
    -o "${SUB_OUT}/${GRAPH}.AMY.gfa"

  # 2-2) build, sort and position the subgraph, odgi v0.9.0
  odgi build -g "${SUB_OUT}/${GRAPH}.AMY.gfa" -o "${SUB_OUT}/${GRAPH}.AMY.og" \
             -t "${THREADS}" -P
  odgi sort -i "${SUB_OUT}/${GRAPH}.AMY.og" -o "${SUB_OUT}/${GRAPH}.AMY.sorted.og" \
            -O -t "${THREADS}"

  odgi extract -i "${SUB_OUT}/${GRAPH}.AMY.sorted.og" \
               -r "${REF_SAMPLE}#0#${AMY_REGION}" \
               -o "${SUB_OUT}/${GRAPH}.AMY.extract.og" -t "${THREADS}"

  odgi view -i "${SUB_OUT}/${GRAPH}.AMY.extract.og" -g \
    > "${SUB_OUT}/${GRAPH}.AMY.extract.gfa"

  # 2-3) position of every node on the reference path, used to attach the gene
  #      annotation to the nodes
  odgi position -i "${SUB_OUT}/${GRAPH}.AMY.extract.og" \
                -r "${REF_SAMPLE}#0#$(echo "${AMY_REGION}" | cut -d: -f1)" \
                -g > "${SUB_OUT}/${GRAPH}.AMY.node_positions.tsv"

  # 2-4) node to gene table, written both as TSV for the tagging step and as
  #      CSV for colouring the graph in Bandage
  awk -F'\t' 'NR>1 && $2 != "" {print $1"\t"$2"\t"$3}' \
    "${SUB_OUT}/${GRAPH}.AMY.node_positions.tsv" \
    | sort -k2,2 -k3,3n > "${TMP_DIR}/node_pos.tsv"

  awk 'BEGIN{OFS="\t"} {print $2, $3-1, $3, $1}' "${TMP_DIR}/node_pos.tsv" \
    | sed "s/^${REF_SAMPLE}#0#//" | sort -k1,1 -k2,2n > "${TMP_DIR}/node_pos.bed"

  bedtools intersect -a "${TMP_DIR}/node_pos.bed" -b "${GENE_BED}" -wa -wb \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $4, $8}' | sort -u \
    > "${SUB_OUT}/${GRAPH}.AMY.node_gene.tsv"

  { echo "Name,Gene"
    awk -F'\t' '{print $1","$2}' "${SUB_OUT}/${GRAPH}.AMY.node_gene.tsv"
  } > "${SUB_OUT}/${GRAPH}.AMY.node_gene.csv"

  # 2-5) figure of the subgraph, Bandage v0.8.1
  Bandage image "${SUB_OUT}/${GRAPH}.AMY.extract.gfa" \
    "${SUB_OUT}/${GRAPH}.AMY.bandage.png" \
    --colours "${SUB_OUT}/${GRAPH}.AMY.node_gene.csv" \
    --height 2000 \
    2> "${LOG_OUT}/bandage.err" || \
    log "  Bandage needs a display; run it interactively if this failed"

  log "step 2 finished"
}

###############################################################################
# Step 3. Haplotype structure of the region
###############################################################################

tag_haplotypes() {
  log "step 3: haplotype tagging"

  SUBGFA="${SUB_OUT}/${GRAPH}.AMY.extract.gfa"
  NODEGENE="${SUB_OUT}/${GRAPH}.AMY.node_gene.tsv"
  [[ -s "${SUBGFA}" && -s "${NODEGENE}" ]] || { log "  subgraph not found"; return 1; }

  # 3-1) genes carried by each haplotype, in walk order and with the
  #      orientation of the nodes (+ forward, - reverse)
  python3 "${SCRIPT_DIR}/amy_haplotype_tags.py" \
    "${SUBGFA}" "${NODEGENE}" "${HAP_OUT}/${GRAPH}.AMY"

  # 3-2) haplotypes sharing the same gene structure are one haplotype class
  awk -F'\t' 'NR>1 {count[$3]++; if (!($3 in first)) first[$3] = $1}
    END {print "haplotype_class\tn_haplotypes\texample\tstructure"
         i = 0
         for (s in count) {i++; printf "H%d\t%d\t%s\t%s\n", i, count[s], first[s], s}}' \
    "${HAP_OUT}/${GRAPH}.AMY.haplotype_structure.tsv" \
    | sort -k2,2nr > "${HAP_OUT}/${GRAPH}.AMY.haplotype_classes.tsv"

  # 3-3) distribution of the haplotypes across the population groups
  if [[ -s "${POP_MAP}" ]]; then
    awk -F'\t' 'BEGIN{OFS="\t"}
      FNR==NR {pop[$1] = $2; next}
      FNR==1 {print "haplotype", "population", "structure"; next}
      { split($1, a, "#"); p = ($0 && (a[1] in pop)) ? pop[a[1]] : "NA"
        print $1, p, $3 }' \
      "${POP_MAP}" "${HAP_OUT}/${GRAPH}.AMY.haplotype_structure.tsv" \
      > "${HAP_OUT}/${GRAPH}.AMY.haplotype_population.tsv"

    awk -F'\t' 'NR>1 {count[$2"\t"$3]++}
      END {print "population\tstructure\tn"; for (k in count) print k"\t"count[k]}' \
      "${HAP_OUT}/${GRAPH}.AMY.haplotype_population.tsv" \
      | sort -k1,1 -k3,3nr \
      > "${STAT_OUT}/${GRAPH}.AMY.haplotype_distribution.tsv"
  else
    log "  ${POP_MAP} not found, population distribution skipped"
  fi

  log "step 3 finished"
}

###############################################################################
# Step 4. Phylogeny of the tagged sequences
###############################################################################

haplotype_tree() {
  log "step 4: haplotype phylogeny"

  SUBGFA="${SUB_OUT}/${GRAPH}.AMY.extract.gfa"
  TAGS="${HAP_OUT}/${GRAPH}.AMY.gene_tags.tsv"
  [[ -s "${SUBGFA}" && -s "${TAGS}" ]] || { log "  inputs not found"; return 1; }

  # 4-1) sequence of every tagged haplotype, assembled from the nodes it
  #      carries in the subgraph
  python3 - "${SUBGFA}" "${TAGS}" "${HAP_OUT}/${GRAPH}.AMY.haplotypes.fa" <<'PY'
import sys
from collections import OrderedDict

gfa, tags, out = sys.argv[1:4]

seq = {}
with open(gfa) as handle:
    for line in handle:
        if line.startswith("S"):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                seq[f[1]] = f[2]

def revcomp(s):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(table)[::-1]

order = OrderedDict()
with open(tags) as handle:
    next(handle)
    for line in handle:
        f = line.rstrip("\n").split("\t")
        order.setdefault(f[0], []).append((f[2], f[3]))

with open(out, "w") as handle:
    for name, steps in order.items():
        pieces = []
        for node, strand in steps:
            s = seq.get(node, "")
            pieces.append(s if strand == "+" else revcomp(s))
        handle.write(">{}\n{}\n".format(name, "".join(pieces)))
PY

  # 4-2) multiple sequence alignment, MAFFT v7
  mafft --auto --thread "${THREADS}" \
    "${HAP_OUT}/${GRAPH}.AMY.haplotypes.fa" \
    > "${HAP_OUT}/${GRAPH}.AMY.haplotypes.aln.fa" \
    2> "${LOG_OUT}/mafft.haplotypes.err"

  # 4-3) phylogenetic tree, IQ-TREE 2
  iqtree2 -s "${HAP_OUT}/${GRAPH}.AMY.haplotypes.aln.fa" \
          -m MFP -B 1000 -T "${THREADS}" \
          --prefix "${HAP_OUT}/${GRAPH}.AMY.haplotypes"

  # 4-4) tree and gene structures side by side
  Rscript "${SCRIPT_DIR}/plot_haplotypes.R" \
    "${TAGS}" \
    "${HAP_OUT}/${GRAPH}.AMY.haplotypes.treefile" \
    "${HAP_OUT}/${GRAPH}.AMY.haplotypes.pdf" \
    2> "${LOG_OUT}/plot_haplotypes.err" || \
    log "  the R packages ape, ggplot2, gggenes and dplyr are required for the figure"

  log "step 4 finished"
}

###############################################################################
# Entry point
###############################################################################

usage() {
  cat <<EOF
Usage: bash $(basename "$0") <step>

  select_regions     bubble selection and gene annotation of the candidates
  extract_subgraph   subgraph extraction, node-gene annotation and Bandage figure
  tag_haplotypes     gene tags per haplotype and haplotype distribution
  haplotype_tree     MAFFT alignment, IQ-TREE 2 phylogeny and figure
  all                run every step in order
EOF
}

main() {
  case "${1:-}" in
    select_regions)   select_regions ;;
    extract_subgraph) extract_subgraph ;;
    tag_haplotypes)   tag_haplotypes ;;
    haplotype_tree)   haplotype_tree ;;
    all)
      select_regions
      extract_subgraph
      tag_haplotypes
      haplotype_tree
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
