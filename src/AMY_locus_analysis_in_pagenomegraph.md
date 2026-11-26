# 7.AMY_locus_analysis_in_pagenomegraph.sh

Identification and visualization of structural haplotypes in the AMY locus of the augmented HPRC+KPR pangenome (Figures 5, S8; Tables S28–S29).

Part of the [KPR-pipeline](../README.md) repository.

## 1. Set up

### 1.1 Software

Make all tools under [Software and versions](#software-and-versions) available on `PATH`. The tree and
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


## 2. Run

~~~bash
bash 7.AMY_locus_analysis_in_pagenomegraph.sh          # show steps
bash 7.AMY_locus_analysis_in_pagenomegraph.sh all
~~~

Steps can also be run separately:

~~~bash
bash 7.AMY_locus_analysis_in_pagenomegraph.sh select_regions
bash 7.AMY_locus_analysis_in_pagenomegraph.sh extract_subgraph
bash 7.AMY_locus_analysis_in_pagenomegraph.sh tag_haplotypes
bash 7.AMY_locus_analysis_in_pagenomegraph.sh haplotype_tree
~~~

The structural subgraph's default reference interval is
`AMY_REGION=chr1:103,570,000-103,760,000`, which is **independent** of
the Bolognini distal/proximal regions analyzed in stage 8.


## 3. Steps

| Step | Description |
|---|---|
| `select_regions` | Select large structural bubbles intersecting AMY genes |
| `extract_subgraph` | Extract graph region with gfabase and odgi; annotate nodes; Bandage figure |
| `tag_haplotypes` | Determine gene order and orientation of every haplotype, group structural classes |
| `haplotype_tree` | Assemble per-haplotype sequences, MAFFT alignment, IQ-TREE 2 and R plotting |
| `all` | Run the four structural steps in order |

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
| `POP_MAP` | `LIST_DIR/sample_population.tsv` | Sample-to-population mapping |
| `AMY_CLUSTER` | `AMY1A AMY1B AMY1C AMYP1 AMY2B AMY2A` | Gene names for tagging |
| `AMY_REGION` | `chr1:103,570,000-103,760,000` | Structural subgraph extraction |
| `REF_SAMPLE` | `GRCh38` | Reference graph path |
| `BUBBLE_MIN_LEN` | `10000` | Structural bubble minimum length |
| `ALLELE_MIN` | `5` | Minimum number of sample alleles |
| `THREADS` | `8` | Threads for supported tools |
| `HAPLOTYPE_MAP` | unset | TSV with structure and haplotype_class columns for study nomenclature |

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

The principal outputs support manuscript Figure 4 and
Supplementary Tables S28–S29: the node/gene structure
map; distinct HxAy-Sz structures; sequence tree;
AMY1 copy-number group counts; and distribution
across KOR, HPRC-EAS and HPRC-nEAS assembly panels.
Check that the total panel count reconciles to
28 + 94 + 2 reference paths, that population
statistics exclude the two references, and that
multiple haplotype paths from the same donor are
not described as independent people. Structural
counts may differ if the graph, assembly panel,
gene interval or class definition changes.

`haplotypes/GRAPH.AMY.haplotype_assignments.tsv` records the class and structure for each path. `haplotype_classes.tsv` contains structure group counts.

## 6. Selection of structural variation loci

We identified bubbles in HPRC+KPR using gfatools v0.5 and obtained per-assembly alleles with Minigraph v0.21 (--cxasm --call). Bubbles >10 kb with ≥5 alleles were retained and overlapping protein-coding genes annotated. The AMY gene-containing GRCh38 interval was chr1:103,554,644–103,758,692; set AMY_REGION to the graph extraction interval used for the subgraph.

## 7. AMY subgraph extraction and gene tagging

We extracted the subgraph with gfabase v0.5, processed node–gene annotations with odgi v0.9.0 and visualized the graph using Bandage v0.8.1. Gene tags record AMY1, AMY2A, AMY2B and AMYP1 copies and the +/− orientation of their nodes along haplotype walks. Inspect the full gene-tag structures and reconcile AMY1 paralog names to the AMY1 copy group before structural classification.

## 8. Structural haplotype classification

We analyzed 124 assemblies/paths (28 KOR, 94 HPRC, GRCh38 and CHM13) and classified structural haplotypes by AMY1 copy number (H1–H7), AMY2A copy number (A0–A2) and orientation/AMYP1 configurations (S1–S6), following the published AMY classification. Supply a two-column TSV with header structure and haplotype_class that maps the gene-tag strings to the verified HxAy–Sz labels in S28. The wrapper otherwise assigns deterministic Structure IDs; these are not the biological H-group names.

~~~bash
HAPLOTYPE_MAP=/path/to/AMY.structure_labels.tsv \
  bash 7.AMY_locus_analysis_in_pagenomegraph.sh tag_haplotypes
~~~

## 9. Sequence alignment and phylogenetic analysis

Align the extracted haplotype sequences using MAFFT, infer the phylogeny using IQ-TREE 2 and visualize with R. We compared classifications of shared assemblies with the previously reported HPRC AMY haplotypes (Table S29).

## 10. Population distribution of AMY haplotypes

Exclude GRCh38 and CHM13 from population counts. Compare 28 KOR, 10 HPRC-EAS and 84 HPRC-nEAS assemblies for structural haplotype membership (Figure 5c). For Figure 5d, group the 122 population assemblies into AFR, AMR, CHN, KOR, EUR, EAS and PJL as in S28; calculate each population's contribution within each H group using the total number of population assemblies in that H group as denominator. KOR-private and CHN-private labels refer to configurations absent from the other sampled populations in this panel. Use the verified HxAy–Sz assignments and population map for the final counts.

## References

1. Kim JY, et al. *Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus.* KPR manuscript, revision 2026; Figure 4 and Tables S28–S29.
2. Bolognini D, et al. (2024). Recurrent evolution and selection shape structural diversity at the amylase locus. *Nature* 634:617–625. https://doi.org/10.1038/s41586-024-07911-1
3. Yilmaz F, et al. (2024). Reconstruction of the human amylase locus reveals ancient duplications seeding modern-day variation. *Science*. https://doi.org/10.1126/science.adn0609
4. Usher CL, et al. (2015). Structural forms of the human amylase locus and their relationships to SNPs, haplotypes and obesity. *Nature Genetics*. https://doi.org/10.1038/ng.3340

## Software and versions

Versions are aligned with the structural haplotype Methods
in the manuscript.

| Task | Software (manuscript version) |
|---|---|
| Structural bubbles and allele calls | gfatools 0.5; Minigraph 0.21 |
| Locus subgraph extraction | gfabase 0.5; odgi 0.9.0 |
| Structural graph visualization | Bandage 0.8.1 |
| Gene-coordinate intersections | bedtools; GENCODE v38 annotation |
| Haplotype sequence alignment | MAFFT v7 |
| Haplotype phylogeny | IQ-TREE 2 |
| Haplotype and tree plotting | R with ape, ggplot2, gggenes and dplyr |
| Supporting scripts | Python 3 |

Preserve the software release, input graph version and sequence-label
conventions for each analysis; graph path changes can affect
structural haplotype classification without any change to the
visualization program.
