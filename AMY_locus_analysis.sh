#!/usr/bin/env bash
#
# Constructing a Korean pangenome reveals distinct haplotypes at the amylase
# locus
#
# AMY_locus_analysis.sh
#
# Usage: bash AMY_locus_analysis.sh <step>
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

# Genes of the amylase cluster. AMY1A, AMY1B and AMY1C share about 99% sequence
# identity and are treated together as AMY1.
AMY_GENES=${AMY_GENES:-"AMY1 AMY2A AMY2B AMYP1"}
AMY_CLUSTER=${AMY_CLUSTER:-"AMY1A AMY1B AMY1C AMYP1 AMY2B AMY2A"}

# Region of the amylase cluster in the coordinates of the reference embedded in
# the graph, used to extract the subgraph.
AMY_REGION=${AMY_REGION:-chr1:103,570,000-103,760,000}
REF_SAMPLE=${REF_SAMPLE:-GRCh38}

# Population groups used for the differentiation analysis. EUR and SAS are
# represented by a single individual each and are excluded.
POPULATIONS=${POPULATIONS:-"AFR AMR EAS KOR"}

# Filters applied to the SNPs called from the gene-copy alignments.
MAF=${MAF:-0.01}

# Windowed FST parameters.
FST_WINDOW=${FST_WINDOW:-2000}
FST_STEP=${FST_STEP:-1000}

# --- reference files ---------------------------------------------------------
# BED file of protein-coding genes with the gene name in column 4, derived from
# the GENCODE annotation
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
GENE_BED=${GENE_BED:-${REF_DIR}/gencode.GRCh38.genes.name.bed}

# FASTA holding the GRCh38 sequence of each amylase gene, one record per gene
# named after the gene, used as the reference of the gene-copy alignments
AMY_REF_FA=${AMY_REF_FA:-${REF_DIR}/AMY_genes.GRCh38.fa}

# --- sample tables -----------------------------------------------------------
# Two columns, sample name and population group
POP_MAP=${POP_MAP:-${LIST_DIR}/sample_population.tsv}

BUBBLE_OUT=${OUT_DIR}/bubble_selection
SUB_OUT=${OUT_DIR}/subgraph
HAP_OUT=${OUT_DIR}/haplotypes
COPY_OUT=${OUT_DIR}/gene_copies
SNP_OUT=${OUT_DIR}/gene_copy_snps
STAT_OUT=${OUT_DIR}/population_statistics
LOG_OUT=${OUT_DIR}/logs

mkdir -p "${BUBBLE_OUT}" "${SUB_OUT}" "${HAP_OUT}" "${COPY_OUT}" "${SNP_OUT}" \
         "${STAT_OUT}" "${LOG_OUT}" "${TMP_DIR}" "${LIST_DIR}"

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
# Step 5. Decomposition of the haplotypes into gene copies
###############################################################################

gene_copies() {
  log "step 5: gene copy extraction"

  SUBGFA="${SUB_OUT}/${GRAPH}.AMY.extract.gfa"
  TAGS="${HAP_OUT}/${GRAPH}.AMY.gene_tags.tsv"
  [[ -s "${SUBGFA}" && -s "${TAGS}" ]] || { log "  inputs not found"; return 1; }

  # Every haplotype carries a variable number of copies of each amylase gene.
  # Each copy is written as a separate sequence and is treated as one
  # independent observation, assigned to the population of its source
  # individual.
  python3 - "${SUBGFA}" "${TAGS}" "${COPY_OUT}" "${AMY_GENES}" <<'PY'
import os
import sys

gfa, tags, outdir, genes = sys.argv[1:5]
genes = genes.split()

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

copies = {g: [] for g in genes}
counter = {}

with open(tags) as handle:
    next(handle)
    for line in handle:
        f = line.rstrip("\n").split("\t")
        haplotype, node, strand, gene = f[0], f[2], f[3], f[4]

        # AMY1A, AMY1B and AMY1C share about 99% identity and are one gene here
        if gene.startswith("AMY1"):
            gene = "AMY1"
        if gene not in copies:
            continue

        key = (haplotype, gene)
        counter[key] = counter.get(key, 0) + 1
        name = "{}|{}|copy{}".format(haplotype, gene, counter[key])

        s = seq.get(node, "")
        copies[gene].append((name, s if strand == "+" else revcomp(s)))

os.makedirs(outdir, exist_ok=True)
for gene, records in copies.items():
    path = os.path.join(outdir, "{}.copies.fa".format(gene))
    with open(path, "w") as handle:
        for name, s in records:
            handle.write(">{}\n{}\n".format(name, s))
    print("{}\t{} copies".format(gene, len(records)))
PY

  # Number of copies per gene and per population group
  if [[ -s "${POP_MAP}" ]]; then
    {
      echo -e "gene\tpopulation\tcopies"
      for G in ${AMY_GENES}; do
        grep "^>" "${COPY_OUT}/${G}.copies.fa" 2>/dev/null \
          | sed 's/^>//' \
          | awk -F'|' '{split($1, a, "#"); print a[1]}' \
          | awk -F'\t' -v g="${G}" 'BEGIN{OFS="\t"}
              FNR==NR {pop[$1]=$2; next}
              {p = ($1 in pop) ? pop[$1] : "NA"; count[p]++}
              END {for (p in count) print g, p, count[p]}' "${POP_MAP}" -
      done
    } > "${STAT_OUT}/${GRAPH}.AMY.copies_per_population.tsv"
  fi

  log "step 5 finished"
}

###############################################################################
# Step 6. SNPs of the gene copies
###############################################################################

copy_snps() {
  log "step 6: SNP calling on the gene copies"

  [[ -s "${AMY_REF_FA}" ]] || { log "  ${AMY_REF_FA} not found"; return 1; }

  {
    echo -e "gene\tcopies\tsites\tlength_bp"

    for G in ${AMY_GENES}; do
      FA="${COPY_OUT}/${G}.copies.fa"
      [[ -s "${FA}" ]] || { log "  ${FA} not found, skipped"; continue; }

      # 6-1) the copies are aligned together with the GRCh38 sequence of the
      #      gene, which provides the coordinates of the call set
      samtools faidx "${AMY_REF_FA}" "${G}" > "${TMP_DIR}/${G}.ref.fa"
      cat "${TMP_DIR}/${G}.ref.fa" "${FA}" > "${TMP_DIR}/${G}.input.fa"

      mafft --auto --thread "${THREADS}" "${TMP_DIR}/${G}.input.fa" \
        > "${SNP_OUT}/${G}.aln.fa" 2> "${LOG_OUT}/mafft.${G}.err"

      # 6-2) substitutions relative to the reference sequence, one haploid
      #      genotype per gene copy
      python3 "${SCRIPT_DIR}/msa_to_vcf.py" \
        "${SNP_OUT}/${G}.aln.fa" "${G}" "${G}" "${SNP_OUT}/${G}.raw.vcf"

      bgzip -f "${SNP_OUT}/${G}.raw.vcf"
      bcftools index -f -t "${SNP_OUT}/${G}.raw.vcf.gz"

      # 6-3) multiallelic records are split, a minor allele frequency filter is
      #      applied and only biallelic sites are kept
      bcftools norm -m -any "${SNP_OUT}/${G}.raw.vcf.gz" -Oz \
        -o "${TMP_DIR}/${G}.split.vcf.gz"
      bcftools index -f -t "${TMP_DIR}/${G}.split.vcf.gz"

      bcftools view -m2 -M2 -v snps -q "${MAF}:minor" \
        "${TMP_DIR}/${G}.split.vcf.gz" -Oz -o "${SNP_OUT}/${G}.biallelic.vcf.gz"
      bcftools index -f -t "${SNP_OUT}/${G}.biallelic.vcf.gz"

      N_COPY=$(bcftools query -l "${SNP_OUT}/${G}.biallelic.vcf.gz" | wc -l)
      N_SITE=$(bcftools index -n "${SNP_OUT}/${G}.biallelic.vcf.gz")
      N_BP=$(awk '/^>/{next} {n += length($0)} END {print n}' "${TMP_DIR}/${G}.ref.fa")
      echo -e "${G}\t${N_COPY}\t${N_SITE}\t${N_BP}"
    done
  } > "${STAT_OUT}/${GRAPH}.AMY.snp_summary.tsv"

  # 6-4) the four genes together
  bcftools concat -a $(for G in ${AMY_GENES}; do echo "${SNP_OUT}/${G}.biallelic.vcf.gz"; done) \
    -Oz -o "${SNP_OUT}/AMY_all.biallelic.vcf.gz" 2> "${LOG_OUT}/concat.err" || \
  bcftools concat $(for G in ${AMY_GENES}; do echo "${SNP_OUT}/${G}.biallelic.vcf.gz"; done) \
    -Oz -o "${SNP_OUT}/AMY_all.biallelic.vcf.gz"
  bcftools index -f -t "${SNP_OUT}/AMY_all.biallelic.vcf.gz"

  log "step 6 finished"
}

###############################################################################
# Step 7. Population differentiation and diversity
###############################################################################

population_stats() {
  log "step 7: population statistics"

  [[ -s "${POP_MAP}" ]] || { log "  ${POP_MAP} not found"; return 1; }

  # 7-1) list of the gene copies belonging to each population group; a copy is
  #      assigned to the population of the individual it comes from
  for P in ${POPULATIONS}; do
    bcftools query -l "${SNP_OUT}/AMY_all.biallelic.vcf.gz" \
      | awk -F'|' -v p="${P}" 'BEGIN{OFS="\t"}
          FNR==NR {pop[$1] = $2; next}
          {split($1, a, "#"); if ((a[1] in pop) && pop[a[1]] == p) print $0}' \
        "${POP_MAP}" - > "${STAT_OUT}/${P}.copies.list"
    log "  ${P}: $(wc -l < "${STAT_OUT}/${P}.copies.list") gene copies"
  done

  # 7-2) pairwise weighted FST in sliding windows, VCFtools v0.1.16. The
  #      weighted estimator combines the per-site numerators and denominators
  #      as a ratio of sums rather than averaging the per-site ratios.
  for G in ${AMY_GENES} all; do
    VCF="${SNP_OUT}/${G}.biallelic.vcf.gz"
    [[ "${G}" == "all" ]] && VCF="${SNP_OUT}/AMY_all.biallelic.vcf.gz"
    [[ -s "${VCF}" ]] || continue

    for P1 in ${POPULATIONS}; do
      for P2 in ${POPULATIONS}; do
        [[ "${P1}" < "${P2}" ]] || continue

        vcftools --gzvcf "${VCF}" \
          --weir-fst-pop "${STAT_OUT}/${P1}.copies.list" \
          --weir-fst-pop "${STAT_OUT}/${P2}.copies.list" \
          --fst-window-size "${FST_WINDOW}" \
          --fst-window-step "${FST_STEP}" \
          --out "${STAT_OUT}/fst.${G}.${P1}_${P2}" \
          2> "${LOG_OUT}/fst.${G}.${P1}_${P2}.log"
      done
    done
  done

  # 7-3) AMY-wide FST with Hudson's estimator, a ratio of averages that is
  #      robust to the unequal group sizes
  {
    echo -e "pop1\tpop2\tsites\tHudson_FST"
    for P1 in ${POPULATIONS}; do
      for P2 in ${POPULATIONS}; do
        [[ "${P1}" < "${P2}" ]] || continue
        python3 "${SCRIPT_DIR}/hudson_fst.py" \
          "${SNP_OUT}/AMY_all.biallelic.vcf.gz" \
          "${STAT_OUT}/${P1}.copies.list" "${STAT_OUT}/${P2}.copies.list" \
          "${P1}" "${P2}" | tail -n +2
      done
    done
  } > "${STAT_OUT}/${GRAPH}.AMY.hudson_fst.tsv"

  # 7-4) nucleotide diversity per group, normalised by the length of the
  #      calling reference of each gene
  {
    echo -e "gene\tpopulation\tsites\tpi_sum\tlength_bp\tpi_per_bp"
    for G in ${AMY_GENES}; do
      VCF="${SNP_OUT}/${G}.biallelic.vcf.gz"
      [[ -s "${VCF}" ]] || continue
      N_BP=$(awk -F'\t' -v g="${G}" '$1==g {print $4}' \
        "${STAT_OUT}/${GRAPH}.AMY.snp_summary.tsv")

      for P in ${POPULATIONS}; do
        vcftools --gzvcf "${VCF}" \
          --keep "${STAT_OUT}/${P}.copies.list" \
          --site-pi --out "${STAT_OUT}/pi.${G}.${P}" \
          2> "${LOG_OUT}/pi.${G}.${P}.log"

        awk -F'\t' -v g="${G}" -v p="${P}" -v n="${N_BP}" 'BEGIN{OFS="\t"}
          NR>1 {sum += $3; sites++}
          END {print g, p, sites+0, sum+0, n, (n > 0) ? sum / n : "NA"}' \
          "${STAT_OUT}/pi.${G}.${P}.sites.pi"
      done
    done
  } > "${STAT_OUT}/${GRAPH}.AMY.nucleotide_diversity.tsv"

  # 7-5) linkage disequilibrium as the mean phased r-squared between the SNP
  #      pairs of each gene, PLINK v2.0
  {
    echo -e "gene\tpairs\tmean_r2"
    for G in ${AMY_GENES}; do
      VCF="${SNP_OUT}/${G}.biallelic.vcf.gz"
      [[ -s "${VCF}" ]] || continue

      plink2 --vcf "${VCF}" --double-id --allow-extra-chr \
             --r2-phased --ld-window-r2 0 \
             --out "${STAT_OUT}/ld.${G}" \
             > "${LOG_OUT}/ld.${G}.log" 2>&1 || true

      LD="${STAT_OUT}/ld.${G}.vcor"
      [[ -s "${LD}" ]] || LD="${STAT_OUT}/ld.${G}.ld"
      [[ -s "${LD}" ]] || { echo -e "${G}\tNA\tNA"; continue; }

      awk 'NR>1 {sum += $NF; n++} END {print (n ? n : 0)"\t"(n ? sum/n : "NA")}' "${LD}" \
        | awk -v g="${G}" 'BEGIN{OFS="\t"} {print g, $1, $2}'
    done
  } > "${STAT_OUT}/${GRAPH}.AMY.linkage_disequilibrium.tsv"

  # AMY-wide linkage disequilibrium is the mean over the four genes
  awk -F'\t' 'NR>1 && $3 != "NA" {sum += $3; n++}
    END {if (n) printf "AMY_wide\t%d\t%.6f\n", n, sum / n}' \
    "${STAT_OUT}/${GRAPH}.AMY.linkage_disequilibrium.tsv" \
    >> "${STAT_OUT}/${GRAPH}.AMY.linkage_disequilibrium.tsv"

  log "step 7 finished"
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
  gene_copies        decomposition of the haplotypes into single gene copies
  copy_snps          alignment of the gene copies and SNP calling
  population_stats   FST, nucleotide diversity and linkage disequilibrium
  all                run every step in order
EOF
}

main() {
  case "${1:-}" in
    select_regions)   select_regions ;;
    extract_subgraph) extract_subgraph ;;
    tag_haplotypes)   tag_haplotypes ;;
    haplotype_tree)   haplotype_tree ;;
    gene_copies)      gene_copies ;;
    copy_snps)        copy_snps ;;
    population_stats) population_stats ;;
    all)
      select_regions
      extract_subgraph
      tag_haplotypes
      haplotype_tree
      gene_copies
      copy_snps
      population_stats
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
