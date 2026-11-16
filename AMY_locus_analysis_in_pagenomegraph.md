# AMY_locus_analysis_in_pagenomegraph.sh

Structural analysis of the AMY region in the combined HPRC+KPR pangenome graph.
This step identifies graph bubbles, extracts the locus subgraph, reconstructs
AMY gene-order haplotypes and infers a haplotype phylogeny.

This script **does not** treat individual AMY gene copies as independent diploid
SNP loci. Copy-level FST, nucleotide diversity and LD were excluded from the
revision. For locus-level short-read SNP analyses of distal and proximal
unique flanks, see [Population genetics of AMY flanks](5.Population_genetics_in_AMY_flanking_locus.md).

Part of the [KPR-pipeline](README.md).

---

## 1. Set up

### 1.1 Software

Make all tools under [Software](#software) available on `PATH`. The tree and
gene-structure figure requires R and `ape`, `ggplot2`, `gggenes`, `dplyr`.

### 1.2 Pangenome graph

The script expects `GRAPH_DIR/GRAPH/GRAPH.gfa` from
[pangenome construction](pangenome_construction_and_statistics.md).
The default `GRAPH` is `hprc_kpr` (combined HPRC+KPR).

### 1.3 Bubble tables and references

The construction pipeline's `BUBBLE_DIR` should contain
`GRAPH.bubble` and `GRAPH.minigraph.merged.call`. A missing bubble file
can be regenerated; the per-sample allele-call table must already exist.

`GENE_BED` is a GENCODE-derived BED of protein-coding genes with gene
names in column 4. `REF_DIR` is required and defaults `GENE_BED` to
`REF_DIR/gencode.GRCh38.genes.name.bed`.

### 1.4 Sample population table

`POP_MAP` (default `LIST_DIR/sample_population.tsv`) optionally maps graph
sample identifiers to population labels. This yields **descriptive counts
of graph haplotype architecture classes**, not paralog-specific
genotype-based population statistics.

### 1.5 Environment

~~~bash
export GRAPH_DIR=/path/to/graphs
export REF_DIR=/path/to/reference
export BUBBLE_DIR=/path/to/bubbles
export OUT_DIR=/path/to/results
export LIST_DIR=/path/to/lists
export THREADS=8
~~~

---

## 2. Run

~~~bash
bash AMY_locus_analysis_in_pagenomegraph.sh          # show steps
bash AMY_locus_analysis_in_pagenomegraph.sh all
~~~

Steps can also be run separately:

~~~bash
bash AMY_locus_analysis_in_pagenomegraph.sh select_regions
bash AMY_locus_analysis_in_pagenomegraph.sh extract_subgraph
bash AMY_locus_analysis_in_pagenomegraph.sh tag_haplotypes
bash AMY_locus_analysis_in_pagenomegraph.sh haplotype_tree
~~~

The structural subgraph's default reference interval is
`AMY_REGION=chr1:103,570,000-103,760,000`, which is **independent** of
the Bolognini distal/proximal regions used by step 7.

---

## 3. Steps

| Step | Description |
|---|---|
| `select_regions` | Select large structural bubbles intersecting AMY genes |
| `extract_subgraph` | Extract graph region with gfabase and odgi; annotate nodes; Bandage figure |
| `tag_haplotypes` | Determine gene order and orientation of every haplotype, group structural classes |
| `haplotype_tree` | Assemble per-haplotype sequences, MAFFT alignment, IQ-TREE 2 and R plotting |
| `all` | Run the four structural steps in order |

The former steps `gene_copies`, `copy_snps` and `population_stats`
are intentionally no longer supported. Existing supporting scripts may remain
in this repository for provenance, but are **not** invoked by `all`.

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `GRAPH_DIR` | required | Directory containing graph folders |
| `REF_DIR` | required | Reference and gene annotation directory |
| `BUBBLE_DIR` | unset | Existing bubble tables |
| `OUT_DIR` | `./results` | Output root |
| `LIST_DIR` | `./lists` | Optional population map location |
| `TMP_DIR` | `OUT_DIR/tmp` | Temporary files |
| `GRAPH` | `hprc_kpr` | Combined graph name |
| `GENE_BED` | `REF_DIR/gencode.GRCh38.genes.name.bed` | Gene BED with names |
| `POP_MAP` | `LIST_DIR/sample_population.tsv` | Structural haplotype labels |
| `AMY_CLUSTER` | `AMY1A AMY1B AMY1C AMYP1 AMY2B AMY2A` | Gene names for tagging |
| `AMY_REGION` | `chr1:103,570,000-103,760,000` | Structural subgraph extraction |
| `REF_SAMPLE` | `GRCh38` | Reference graph path |
| `BUBBLE_MIN_LEN` | `10000` | Structural bubble minimum length |
| `ALLELE_MIN` | `5` | Minimum number of sample alleles |
| `THREADS` | `8` | Threads for supported tools |

---

## 5. Output

~~~text
results/
├── bubble_selection/          selected structural bubbles, gene annotations
├── subgraph/                  AMY GFA, node annotation, Bandage figure
├── haplotypes/                gene order, structure classes, sequences, tree
├── haplotype_distribution/    optional structural haplotype-class counts
└── logs/                      pipeline and per-tool stderr
~~~

The results characterize **graph-based locus architecture** and are not
the same as π, θW, Tajima's D, LD, nSL, Hudson FST or XP-nSL from
the separate flank-level population pipeline.

---

## Software

gfatools 0.5; Minigraph 0.21; gfabase 0.5; odgi 0.9.0;
Bandage 0.8.1; MAFFT 7; IQ-TREE 2; bedtools; Python 3;
R with ape, ggplot2, gggenes and dplyr.