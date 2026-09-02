# AMY_locus_analysis.sh

Selection of the structurally variable regions of the pangenome, extraction and
visualisation of the amylase locus subgraph, reconstruction of the haplotype
structures, and population differentiation analysis of the amylase gene copies.

Part of the [KPR-pipeline](README.md) repository.

The analysis runs on the combined HPRC+KPR graph. Candidate regions are chosen
from the bubbles of the graph, the amylase locus is extracted as a subgraph,
each haplotype is described by the genes it carries and their orientation, and
the individual gene copies are then used to compute FST, nucleotide diversity
and linkage disequilibrium between population groups.

---

## 1. Set up

### 1.1 Software

Install the tools listed under [Software](#software) and make sure they are on
`PATH`. The figure of step 4 additionally needs the R packages `ape`,
`ggplot2`, `gggenes` and `dplyr`.

### 1.2 Pangenome graph

`GRAPH_DIR` holds one subdirectory per graph, as produced by
[`pangenome_construction_and_statistics.sh`](pangenome_construction_and_statistics.md):

```
graphs/
└── hprc_kpr/
    └── hprc_kpr.gfa
```

### 1.3 Bubble tables

Step 1 reuses two files written by the construction pipeline, which are looked
up in `BUBBLE_DIR`:

```
bubbles/
├── hprc_kpr.bubble                     gfatools bubble output
└── hprc_kpr.minigraph.merged.call      merged per-sample allele calls
```

If the bubble file is absent it is regenerated, but the merged allele table is
not: run the `allele_mapping` step of the construction pipeline first.

### 1.4 Reference files

- `GENE_BED` — BED file of protein-coding genes with the gene name in column 4,
  derived from the GENCODE annotation
- `AMY_REF_FA` — FASTA holding the GRCh38 sequence of each amylase gene, one
  record per gene named after the gene (`AMY1`, `AMY2A`, `AMY2B`, `AMYP1`),
  indexed with `samtools faidx`. It provides the coordinates of the SNP call
  set of each gene.

### 1.5 Population map

`LIST_DIR/sample_population.tsv`, two columns, the sample name as it appears in
the graph and its population group:

```
HG00733	AMR
HG02257	AFR
KOR-101	KOR
```

Population statistics are computed for AFR, AMR, EAS and KOR. EUR and SAS are
represented by a single individual each and are excluded; EAS comprises CHN and
KHV.

### 1.6 Environment

```bash
export GRAPH_DIR=/path/to/graphs         # required
export REF_DIR=/path/to/reference        # required
export BUBBLE_DIR=/path/to/bubbles       # bubble tables from the construction run
export OUT_DIR=/path/to/results          # default: ./results
export LIST_DIR=/path/to/lists           # default: ./lists
export THREADS=8
```

---

## 2. Run

Print the available steps:

```bash
bash AMY_locus_analysis.sh
```

Run everything in order:

```bash
bash AMY_locus_analysis.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash AMY_locus_analysis.sh select_regions
bash AMY_locus_analysis.sh extract_subgraph
bash AMY_locus_analysis.sh tag_haplotypes
```

Override individual settings on the command line:

```bash
BUBBLE_MIN_LEN=10000 ALLELE_MIN=5 bash AMY_locus_analysis.sh select_regions
AMY_REGION=chr1:103,570,000-103,760,000 bash AMY_locus_analysis.sh extract_subgraph
FST_WINDOW=5000 FST_STEP=1000 bash AMY_locus_analysis.sh population_stats
```

The same script analyses another locus by changing `AMY_REGION` and the gene
lists; the steps themselves are not specific to the amylase cluster.

Bandage draws the graph through a graphical toolkit. On a headless machine the
image step fails with a message in `logs/bandage.err`; the extracted GFA can be
opened in Bandage interactively instead, and the rest of the pipeline is
unaffected.

---

## 3. Steps

| Step | Description |
|---|---|
| `select_regions` | Bubbles longer than 10 kb carrying at least five per-sample alleles, annotated with the protein-coding genes they overlap |
| `extract_subgraph` | Subgraph of the region with `gfabase`, node ordering and positions with `odgi`, node-gene annotation, Bandage figure |
| `tag_haplotypes` | Genes carried by each haplotype in walk order with node orientation, haplotype classes, distribution across the population groups |
| `haplotype_tree` | Sequence of each haplotype, MAFFT alignment, IQ-TREE 2 phylogeny, tree and gene structure figure |
| `gene_copies` | Decomposition of each haplotype into its individual gene copies |
| `copy_snps` | Alignment of the gene copies against the GRCh38 gene sequence and SNP calling |
| `population_stats` | Windowed and AMY-wide FST, nucleotide diversity, linkage disequilibrium |
| `all` | All of the above, in order |

Dependencies: every step after the first needs the output of the previous one;
`select_regions` needs the bubble tables of the construction pipeline.

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `GRAPH_DIR` | required | pangenome graph directory |
| `REF_DIR` | required | reference file directory |
| `BUBBLE_DIR` | unset | bubble tables from the construction pipeline |
| `OUT_DIR` | `./results` | output directory |
| `LIST_DIR` | `./lists` | sample tables |
| `SCRIPT_DIR` | script directory | helper scripts |
| `GRAPH` | `hprc_kpr` | graph to analyse |
| `BUBBLE_MIN_LEN` | `10000` | minimum bubble length |
| `ALLELE_MIN` | `5` | minimum number of per-sample alleles |
| `AMY_REGION` | `chr1:103,570,000-103,760,000` | region extracted as a subgraph |
| `REF_SAMPLE` | `GRCh38` | reference path of the graph |
| `AMY_GENES` | `AMY1 AMY2A AMY2B AMYP1` | genes analysed |
| `AMY_CLUSTER` | `AMY1A AMY1B AMY1C AMYP1 AMY2B AMY2A` | gene names of the cluster |
| `POPULATIONS` | `AFR AMR EAS KOR` | groups compared |
| `MAF` | `0.01` | minor allele frequency filter |
| `FST_WINDOW`, `FST_STEP` | `2000`, `1000` | sliding window for FST |
| `GENE_BED`, `AMY_REF_FA`, `POP_MAP` | see script | annotation inputs |
| `THREADS` | `8` | threads per job |

---

## 5. What the steps do

### Region selection

The bubbles of the graph are filtered to those longer than 10 kb carrying at
least five per-sample alleles, and intersected with the protein-coding gene
annotation. The bubbles overlapping the amylase cluster are written separately.
The per-sample alleles come from `minigraph -cxasm --call`, run in the
construction pipeline.

### Subgraph and figure

`gfabase load` indexes the whole graph and `gfabase sub` cuts out the region.
`odgi build`, `sort`, `extract`, `position` and `view` order the subgraph and
give the position of every node on the reference path, which is intersected
with the gene annotation to produce the node-to-gene table. The table is
written both as TSV, used by the tagging step, and as CSV, used to colour the
nodes in the Bandage figure.

### Haplotype structure

`amy_haplotype_tags.py` reads the W lines of the subgraph and the
node-to-gene table and writes, for each haplotype, the ordered list of genes it
carries with the orientation of the nodes (`+` forward, `-` reverse).
Consecutive nodes of the same gene in the same orientation are one copy.
Haplotypes sharing a structure form one haplotype class, and the classes are
counted per population group.

### Phylogeny

The sequence of each haplotype is assembled from the nodes it traverses,
aligned with MAFFT and used to infer a tree with IQ-TREE 2.
`plot_haplotypes.R` draws the tree next to the gene structures with
`ape` and `gggenes`.

### Gene copies and SNPs

Each haplotype carries a variable number of copies of each amylase gene. Every
copy is written as a separate sequence and treated as one independent
observation assigned to the population of its source individual. AMY1A, AMY1B
and AMY1C share about 99% sequence identity and are treated together as AMY1.

Per gene, the copies are aligned with MAFFT together with the GRCh38 sequence
of that gene, which provides the coordinates.
`msa_to_vcf.py` turns the alignment into a VCF of substitutions, with
one haploid genotype per gene copy. Multiallelic records are split, a minor
allele frequency filter of 1% is applied and only biallelic SNPs are kept.

### Population statistics

Pairwise weighted FST is computed with VCFtools in 2 kb windows with a 1 kb
step. The weighted estimator combines the per-site numerators and denominators
as a ratio of sums rather than averaging the per-site ratios.

AMY-wide FST across the four genes uses Hudson's estimator, implemented in
`hudson_fst.py`. It is also a ratio of averages and is robust to
unequal sample sizes, which matters here because the groups differ by more than
fourfold. For each biallelic site with frequencies p1, p2 and sample sizes n1,
n2 the numerator is `(p1-p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)` and the
denominator is `p1(1-p2) + p2(1-p1)`; the sums are divided at the end.

Nucleotide diversity is computed per group with VCFtools on the same biallelic
SNP set and normalised by the length of the calling reference of each gene.
Linkage disequilibrium is the mean phased r² between the SNP pairs within a
gene, computed with PLINK, and the AMY-wide value is the mean across the four
genes.

---

## 6. Output

```
results/
├── bubble_selection/       selected bubbles and their gene annotation
├── subgraph/               subgraph GFA, node positions, node-gene tables, Bandage figure
├── haplotypes/             gene tags, haplotype structures, alignment, tree, figure
├── gene_copies/            one FASTA of copies per gene
├── gene_copy_snps/         alignments and biallelic SNP call sets
├── population_statistics/  FST, diversity, linkage disequilibrium
└── logs/
```

Main tables:

- `haplotypes/GRAPH.AMY.haplotype_structure.tsv` — gene structure of every
  haplotype
- `haplotypes/GRAPH.AMY.haplotype_classes.tsv` — distinct structures and how
  many haplotypes carry each
- `population_statistics/GRAPH.AMY.haplotype_distribution.tsv` — structures per
  population group
- `population_statistics/GRAPH.AMY.copies_per_population.tsv` — number of gene
  copies per gene and group
- `population_statistics/GRAPH.AMY.snp_summary.tsv` — copies, sites and length
  per gene
- `population_statistics/fst.*.windowed.weir.fst` — windowed pairwise FST
- `population_statistics/GRAPH.AMY.hudson_fst.tsv` — AMY-wide pairwise FST
- `population_statistics/GRAPH.AMY.nucleotide_diversity.tsv`
- `population_statistics/GRAPH.AMY.linkage_disequilibrium.tsv`

---

## Software

gfatools 0.5 · Minigraph 0.21 · gfabase 0.5 · odgi 0.9.0 · Bandage 0.8.1 ·
MAFFT 7 · IQ-TREE 2 · VCFtools 0.1.16 · PLINK 2.0 · bedtools · bcftools ·
samtools · Python 3 · R with ape, ggplot2, gggenes and dplyr

All executables must be available on `PATH`.
