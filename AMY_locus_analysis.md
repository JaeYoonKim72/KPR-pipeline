# AMY_locus_analysis.sh

Structural characterization of the AMY locus in the combined HPRC+KPR graph:
candidate bubbles, subgraph extraction, gene-order haplotype tagging, and
haplotype phylogeny. This pipeline **does not** estimate population-genetic
statistics for individual AMY1A/AMY1B/AMY1C or AMY2A/AMY2B copies.

Population-genetic tests are performed on two unique flanking regions instead,
using [step 5](5.Population_genetics_in_AMY_flanking_locus.md) to avoid
paralog-specific short-read mapping and gene-conversion confounding.

## Setup

The graph is in `GRAPH_DIR/GRAPH/GRAPH.gfa`; `BUBBLE_DIR` optionally points
to bubble tables from the pangenome construction pipeline. Provide `REF_DIR`
containing `GENE_BED` (BED with gene names in column 4), and `POP_MAP`
if population labels are desired for **haplotype-structure counts only**.
These counts are not estimates of per-paralog allele frequencies or FST.

```bash
export GRAPH_DIR=/path/to/graphs
export REF_DIR=/path/to/reference
export BUBBLE_DIR=/path/to/bubbles
export OUT_DIR=/path/to/results
export LIST_DIR=/path/to/lists
export THREADS=8
bash AMY_locus_analysis.sh all
```

## Steps

| Step | Description |
|---|---|
| `select_regions` | Identify candidate structural bubbles overlapping AMY annotations |
| `extract_subgraph` | Extract AMY subgraph; annotate and visualize its nodes |
| `tag_haplotypes` | Recover AMY copy-order/orientation structures from assembly paths |
| `haplotype_tree` | Align locus haplotypes, infer phylogeny, plot gene order |
| `all` | Run these four steps in order |

`bash AMY_locus_analysis.sh <step>` selects one step. The default
`AMY_REGION=chr1:103,570,000-103,760,000` defines the structural subgraph
extraction and is **not** the distal/proximal flanking analysis interval.
Flanking regions are independently configurable in step 5.

## Outputs

- `bubble_selection/`: candidate bubbles and annotations
- `subgraph/`: AMY GFA, node positions and gene tags, Bandage image
- `haplotypes/`: gene-order structures, classes, sequences and tree
- `haplotype_distribution/`: descriptive counts of structural haplotype
  classes by sample group, when `POP_MAP` is provided
- `logs/`: processing logs

No copy-derived SNP calling, per-copy FST, per-copy nucleotide diversity,
or per-copy LD is run. Obsolete step names `gene_copies`, `copy_snps`
and `population_stats` are intentionally no longer accepted.

## Dependencies

gfatools, Minigraph, gfabase, odgi, Bandage, MAFFT, IQ-TREE 2,
bedtools, Python 3, R and `ape`, `ggplot2`, `gggenes`, `dplyr`.
