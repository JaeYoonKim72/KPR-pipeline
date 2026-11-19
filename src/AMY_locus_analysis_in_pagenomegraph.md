# 7.AMY_locus_analysis_in_pagenomegraph.sh

## Scope and rationale

The chromosome 1 amylase (AMY) locus consists of the salivary
amylase gene family AMY1, the pancreatic amylase genes AMY2A
and AMY2B, and the AMYP1 pseudogene. Large near-identical
tandem and segmental duplications create a broad range of
copy-number states, gene orientations and structural
haplotypes. The **structurally variable region (SVR)** can
change substantially in size between haplotypes and
is difficult to reconstruct with one linear reference.
Long-read haplotype assemblies and graph-based paths
allow the physical ordering and orientation of the
repeated segments to be examined directly [1–3].

This stage identifies AMY-related structural bubbles
in the augmented HPRC+KPR graph, extracts the
corresponding subgraph, annotates AMY gene
segments, reconstructs ordered gene-copy patterns
and summarizes the structural haplotypes represented
by the input assemblies. MAFFT alignment and IQ-TREE 2
provide a sequence-based phylogenetic description
of the resulting haplotypes. The tree organizes
sequence relationships; it does not itself prove
the timing, selection history or functional impact
of a duplication event.

**Analysis boundary.** This script does not
treat individual AMY gene copies as independent
diploid SNP loci. Copy-level FST, nucleotide
diversity and LD from separated paralogs are
not part of the revised analysis because
variable copy number and potential interlocus
gene conversion complicate independent allelic
sampling. For locus-level population statistics
using unique distal/proximal sequences, see
[Population genetics of AMY flanks](Population_genetics_in_AMY_flanking_locus.md).

Part of the [KPR-pipeline](../README.md).

---

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

---

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

---

## 3. Steps

| Step | Description |
|---|---|
| `select_regions` | Select large structural bubbles intersecting AMY genes |
| `extract_subgraph` | Extract graph region with gfabase and odgi; annotate nodes; Bandage figure |
| `tag_haplotypes` | Determine gene order and orientation of every haplotype, group structural classes |
| `haplotype_tree` | Assemble per-haplotype sequences, MAFFT alignment, IQ-TREE 2 and R plotting |
| `all` | Run the four structural steps in order |

### 3.1 Selection of a structurally complex locus

Large bubbles provide candidates for multi-allelic
or copy-number-variable sequence structures.
The pipeline selects bubbles >10 kb with at
least five per-sample allele observations,
intersects them with protein-coding gene
annotations and prioritizes the AMY gene
cluster. These are *graph structural*
selection criteria, not variant-allele
frequency thresholds in a diploid population.
The selected bubble interval and reference
path coordinate must be recorded together:
the script's extraction interval
`chr1:103570000-103760000` is an adjustable
subgraph boundary, not the unique distal/
proximal windows of Bolognini et al. [2].

### 3.2 Structural subgraph reconstruction

The graph is extracted with gfabase and
converted for odgi analysis. Node/edge
orientation, path traversal and
gene overlaps determine the order of
AMY segments for each haplotype.
Inversions change gene orientation,
whereas duplications alter repeat
unit counts; both must be retained
when defining a structural class.
A path missing part of the requested
region or an ambiguous gene annotation
requires manual review, rather than
automatic assignment to the nearest
copy-number class. Bandage provides
a structural overview, not a genome
assembly validity check.

### 3.3 AMY structural haplotype naming

The manuscript grouped haplotypes
first by AMY1 copy number (H1–H7),
then by AMY2A copy number (A0–A2),
and subsequently by gene orientation
and AMYP1 presence (S subtypes).
For example, H5A1-S5 denotes a
particular AMY1/AMY2A copy-count
combination with a defined
substructure; the subtype must be
read from the gene-order annotation.
Structural configurations that share
AMY1 counts should not be collapsed
if AMY2A content or gene direction
differs. The manuscript reported
no H6 assembly in its 124-path panel;
that observation does not imply
that H6 cannot exist in other samples.

### 3.4 Haplotype sequence relationships

MAFFT aligns the reconstructed
haplotype-associated sequences, and
IQ-TREE 2 estimates a nucleotide-based
phylogenetic tree. Tree tips must
retain the source path, population,
sample identifier and structural
HxAy-Sz label. Similar structural
arrangements can evolve independently
on distinct SNP backgrounds,
as documented by Bolognini et al.;
consequently, a structural tree and
a single-copy flanking SNP tree
address related but different
evolutionary questions [2,3].
The manuscript observed matching
AMY1 H-group placement among
72 of 85 shared haplotypes in a
comparison to the prior HPRC analysis,
which supports broad comparability
without implying identical
subtype definitions.

### 3.5 Population-panel description and limitations

Population annotations summarize
how frequently each **assembly path**
carries an observed structural class.
The relevant denominators are
the number of haplotype assemblies
in each panel, not the number of
individual gene copies. Reference
paths GRCh38 and CHM13 are included
in structure comparison but excluded
from population-panel frequencies.
The small HPRC East Asian panel
and the absence of independently
measured function prevent claims
that a structural class is fixed,
uniquely Korean or itself under
diet-associated selection.
The separate population-genetics
pipeline evaluates diploid SNPs
in the unique SVR flanks.

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

The principal outputs support manuscript Figure 4 and
Supplementary Tables S23–S24: the node/gene structure
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

---

## Study structural haplotype results

The KPR manuscript evaluated **124 assembly/reference
paths**: 28 Korean haplotype assemblies, 94 HPRC
v1.0 haplotype assemblies, and the GRCh38 and
CHM13 reference sequences. This analysis identified
31 structural architectures, 17 of which were
represented among the 28 Korean assemblies.
Six haplotype architectures were observed only
in the sampled Korean assembly panel, not in
the HPRC comparisons. Their absence from those
relatively small comparison panels should not
be interpreted as population-wide exclusivity.
The structural panel also remains distinct
from the separate 302-person diploid SNP
cohort used for flanking population genetics [1].

---

## References

1. Kim JY, et al. *Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus.* KPR manuscript, revision 2026; Figure 4 and Tables S23–S24.
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
