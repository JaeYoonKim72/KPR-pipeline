# 3.pangenome_construction_and_statistics.sh

## Scope and rationale

A single linear reference captures only one sequence arrangement
at a given locus. Differences in population haplotypes,
insertions, deletions and segmental duplications are therefore
incompletely represented by one coordinate path. A pangenome
graph instead represents sequences as nodes and connections
as edges, with individual haplotypes traversing graph paths.
Graph construction makes shared and alternative sequence
arrangements directly comparable, while preserving the
provenance of each input assembly [1–3].

This workflow builds two Minigraph–Cactus graphs from
quality-controlled, haplotype-resolved assemblies: **KPR**,
containing 28 Korean haplotype assemblies from 14 individuals;
and **HPRC+KPR**, combining those 28 assemblies with the
94 HPRC v1.0 haplotype assemblies, for 122 population
haplotypes. CHM13 and GRCh38 reference paths may also be
present in graph files, and must not be mistaken for
additional sampled people. The CHM13 sequence serves as
the main graph backbone for comparison; GRCh38 provides
a familiar annotation and variant-comparison coordinate
system. Input sample names, path names and reference
path membership must remain consistent [1].

The graph analyses quantify how additional KOR sequences
change the available representation of genomic variation:
graph nodes, edges and path complexity; SNPs, MNPs,
indels and structural variants; shared versus
population-panel-specific variants; pangenome growth;
non-reference sequence accumulation; frequency classes
and repeat composition. Graph-based estimates provide
a representation of observed assemblies, not a direct
estimate of population allele frequency unless sampling
and denominators are explicitly accounted for.

Part of the [KPR-pipeline](../README.md) repository.

---

## 1. Set up

### 1.1 Software

Install the tools listed under [Software and versions](#software-and-versions) and make sure they are on
`PATH`. Cactus is distributed as a virtual environment and must be activated
before the `construct` step:

```bash
source /path/to/cactus-bin-v2.9.0/venv-cactus-v2.9.0/bin/activate
cactus-pangenome --help     # check that it runs
```

### 1.2 Reference files

Download the reference files into one directory. The URLs are given in the
configuration block at the top of the script.

```bash
mkdir -p reference && cd reference

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip -d chm13v2.0.fa.gz && samtools faidx chm13v2.0.fa

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38.gfa.gz
gzip -d hprc-v1.1-mc-grch38.gfa.gz
```

The `bubble_annotation` step additionally needs two BED files with the gene
name in column 4, and a list of medically relevant gene symbols:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gzip -d gencode.v38.annotation.gtf.gz

awk -F'\t' '$3=="gene"  {match($9,/gene_name "[^"]+"/); g=substr($9,RSTART+11,RLENGTH-12);
                         print $1"\t"$4-1"\t"$5"\t"g}' gencode.v38.annotation.gtf \
  | sort -k1,1 -k2,2n > gencode.GRCh38.genes.name.bed

awk -F'\t' '$3=="exon"  {match($9,/gene_name "[^"]+"/); g=substr($9,RSTART+11,RLENGTH-12);
                         print $1"\t"$4-1"\t"$5"\t"g}' gencode.v38.annotation.gtf \
  | sort -k1,1 -k2,2n > gencode.GRCh38.exons.name.bed
```

### 1.3 Seqfiles

Seqfiles are prepared by the user, one per graph, named after the graph
(`kpr.seqfile`, `hprc_kpr.seqfile`) and placed in `SEQFILE_DIR`. Each line
gives an assembly name and the path to its FASTA file; the two haplotypes of a
sample are distinguished by a `.1` and `.2` suffix:

```
CHM13	/path/to/chm13v2.0.fa
GRCh38	/path/to/GRCh38.fa
KOR-101.1	/path/to/KOR-101.hap1.fa
KOR-101.2	/path/to/KOR-101.hap2.fa
```

Contig names inside each FASTA follow the PanSN convention
(`SAMPLE#HAPLOTYPE#CONTIG`). Assemblies are aligned in the order in which they
are listed, so order the seqfile by population group (AFR, EUR, AMR, EAS, KOR),
broadly reflecting the order of human dispersal.

### 1.4 Assemblies

`ASM_DIR` holds the haplotype assembly FASTA files used by `allele_mapping`,
one file per haplotype. Every `*.fa` file in that directory is processed.

### 1.5 Environment

```bash
export SEQFILE_DIR=/path/to/seqfiles     # required
export REF_DIR=/path/to/reference        # required
export ASM_DIR=/path/to/assemblies       # required by allele_mapping
export OUT_DIR=/path/to/results          # default: ./results
export THREADS=16
```

---

## 2. Run

Print the available steps:

```bash
bash 3.pangenome_construction_and_statistics.sh
```

Run everything in order:

```bash
bash 3.pangenome_construction_and_statistics.sh all
```

Run one step at a time, which is what we recommend:

```bash
bash 3.pangenome_construction_and_statistics.sh construct
bash 3.pangenome_construction_and_statistics.sh index_graph
bash 3.pangenome_construction_and_statistics.sh graph_stats
```

Graph construction is long-running: on a single node with 27 CPU cores and
1.5 TB of RAM, the KPR graph took about 11 days and the combined HPRC+KPR graph
about 22 days of wall-clock time. Detach it from the terminal:

```bash
nohup bash 3.pangenome_construction_and_statistics.sh construct > construct.log 2>&1 &
```

Build and analyse a single graph by restricting `GRAPHS`:

```bash
GRAPHS=kpr bash 3.pangenome_construction_and_statistics.sh construct
GRAPHS=hprc_kpr bash 3.pangenome_construction_and_statistics.sh construct
```

Override individual settings on the command line:

```bash
THREADS=32 CACTUS_MEM=1200Gi bash 3.pangenome_construction_and_statistics.sh construct
SV_MIN_LEN=50 bash 3.pangenome_construction_and_statistics.sh variants
```

Every step writes to `$OUT_DIR/logs/pipeline.log` and skips a graph whose input
files are missing, so a step can be rerun after fixing a problem without
redoing the graph. Cactus keeps its job store under `$TMP_DIR`; remove it only
once construction has finished.

---

## 3. Steps

| Step | Description |
|---|---|
| `construct` | Minigraph–Cactus graph construction |
| `index_graph` | vg conversion and indexing (vg, xg, r-index, haplotype index, GBWT, path lists), odgi build and sort, node FASTA |
| `graph_stats` | Nodes, edges, paths and steps from `odgi stats`; pangenome size and saturation curve from `odgi heaps` (200 permutations) |
| `variants` | `vg deconstruct` against the reference paths and against every sample and haplotype, size and type splitting, allele counts, variant classification |
| `growth` | Panacus growth curves and node frequency classes, applied to both graphs and to the HPRC v1.1 graph |
| `allele_mapping` | `minigraph --call` of each haplotype assembly against the graph, merged into a per-bubble allele table |
| `nonref_nodes` | Bubbles from `gfatools bubble`, non-reference SV nodes, frequency classes and repeat content |
| `bubble_annotation` | Gene, exon and medically relevant gene annotation of the large bubbles, and target selection |
| `nonref_sequence` | Sequence carried by each sample, the non-reference part of it, and the sequence of each frequency class |
| `all` | All of the above, in order |

Dependencies between steps: `index_graph`, `graph_stats`, `growth`,
`allele_mapping`, `nonref_nodes` and `nonref_sequence` all need `construct`;
`variants` needs `index_graph`; the frequency classes in `nonref_nodes` and the
whole of `bubble_annotation` need `allele_mapping`.

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `SEQFILE_DIR` | required | directory holding the seqfiles |
| `REF_DIR` | required | reference files |
| `ASM_DIR` | unset | haplotype assemblies, required by `allele_mapping` |
| `OUT_DIR` | `./results` | output directory |
| `TMP_DIR` | `$OUT_DIR/tmp` | Cactus job store and working directory |
| `SCRIPT_DIR` | script directory | helper scripts |
| `GRAPHS` | `kpr hprc_kpr` | graphs to build and analyse |
| `REFERENCE` | `GRCh38` | backbone reference; several names may be given |
| `REF_PATHS` | `GRCh38 CHM13` | reference path names in the graph |
| `FILTER` | `9` | haplotype frequency filter passed to `cactus-pangenome` |
| `SV_MIN_LEN` | `50` | SV if absolute REF/ALT length difference >50 bp; small variant if ≤50 bp |
| `BUBBLE_LARGE` | `10000` | large bubbles retained for annotation |
| `ALLELE_MIN` | `5` | minimum assemblies carrying an allele, medical targets |
| `ALLELE_MAX` | `29` | bubbles carried by nearly every assembly |
| `THREADS` | `16` | cores for each Cactus phase and the downstream tools |
| `CACTUS_MEM` | `800Gi` | `--consMemory` |
| `CACTUS_DISK` | `1.3Ti` | `--maxDisk` |
| `GENE_BED`, `EXON_BED`, `MEDICAL_GENES` | see script | annotation inputs |

Soft-masking of centromeric and telomeric regions with dna-brnn, the MAPQ < 5
alignment filter, the exclusion of soft-masked segments longer than 100 kb, the
removal of unaligned paths longer than 10 kb, `hal2vg` conversion and GFAffix
normalisation are internal Minigraph–Cactus stages and are not overridden.

---

## 5. What the analysis steps do

### 5.1 Input identity and graph construction checks

Use the two QC-approved haplotype FASTA files for each
KOR individual from the preceding assembly workflow.
The HPRC haplotypes are treated as additional assembly
paths, not as separate individuals when the two paths
belong to the same diploid donor. Confirm that the
seqfile contains the intended populations, both haplotypes
where expected, and one entry for each assembly path.
The sample name and haplotype suffix must agree with
the names encoded in the FASTA. A mismatched name
changes downstream sample grouping and sharing
statistics even when graph construction succeeds.

Minigraph–Cactus aligns input paths to the reference-backed
graph, resolves nucleotide-level differences and emits
normalized graph/VCF representations. The input sequence
ordering is held constant across study graphs for
reproducibility, but the order reflects an analysis
configuration and must not be interpreted as an
evolutionary branching order. The reference backbone,
software versions, path ordering and filtering choices
are part of graph provenance.

### 5.2 Graph structure and saturation

`odgi stats` reports nodes, edges, paths and steps;
graph length is not simply the summed sizes of all
input haplotypes, because shared sequence is
represented once in the graph. `odgi heaps`
estimates how pangenome content changes with the
number of sampled paths using 200 permutations.
The graph statistics should be compared on identical
node/edge conventions and comparable filtered or
unfiltered graph states. Growth trajectories require
the same reference definition and sampling unit
when KPR and HPRC are compared.

### 5.3 Population-panel and sample-level comparisons

Small variants and structural variants are compared
between HPRC and KPR paths within the combined graph.
A variant absent from a panel is **panel-specific
within the sampled assemblies**, not necessarily
private to that population worldwide. Per-sample
variant counts distinguish an individual/haplotype
record from a variant site shared by multiple paths.
African samples may contribute larger numbers of
observed variants for demographic and ascertainment
reasons; compare populations with appropriate
denominators and document unequal panel sizes.

### Variant classification

`vg deconstruct -P <path> -a` is run three ways: against each reference path,
against each sample, and against each haplotype. The path lists come from
`vg paths -L` in the `index_graph` step.

Reference-based call sets are also split with `bcftools` into SNPs, MNPs and
indels with duplicate records removed, and the number of sites carrying one,
two or more alternative alleles is counted.

`decon_split.sh` splits call sets by allele size and type. An SV allele
has an absolute REF/ALT length difference **>50 bp**; alleles with a
difference **≤50 bp** are small variants. A site enters the all-SV class
only when every alternative allele meets the SV threshold. Sites mixing
substitutions with length-changing alleles are placed in `only_mixed` files.
The `all_50_alts` output suffix is retained for compatibility and does not
mean that 50-bp differences are included as SVs. The final files
are:

```
*.var_sizes.Smalls_list.no_mixed.SNP.SNP.vcf.gz
*.var_sizes.Smalls_list.no_mixed.INS.INS.vcf.gz
*.var_sizes.Smalls_list.no_mixed.DEL.DEL.vcf.gz
*.var_sizes.SVs_list.all_50_alts.no_mixed.INS.INS.vcf.gz
*.var_sizes.SVs_list.all_50_alts.no_mixed.DEL.DEL.vcf.gz
```

`typecount.py` classifies every record of these files by the number of
alternative alleles carried across all samples: singleton (one), doubleton
(two), common (more than two, frequency below 95%), core (more than two,
frequency of 95% or above). The output column `polymorphism` is defined
in this helper as the sum of doubleton, common and core **records**;
it excludes singletons and is not a count of all variable sites.

Either script can also be run on a single file:

```bash
bash src/decon_split.sh results/graph_variants/kpr/reference/kpr.CHM13.vcf.gz 50
python3 src/typecount.py results/graph_variants/kpr/reference/kpr.CHM13.var_sizes.Smalls_list.no_mixed.SNP.SNP.vcf.gz
```

### Bubbles and non-reference nodes

`gfatools bubble` lists the bubbles of the graph. Bubbles spanning a single
position are discarded as substitutions, and those carrying an allele longer
than 50 bp are kept as structural. Nodes absent from every reference walk are
the non-reference SV nodes; their sequences are annotated for repeat content
with the same two RepeatMasker rounds used for the assemblies.

`minigraph -cxasm --call` records the allele each haplotype assembly carries at
every bubble. Merging the per-assembly calls gives, for each bubble, the number
of assemblies carrying an allele, which is used both for the frequency classes
and for the annotation step.

The large bubbles are intersected with the gene and exon BED files, flagged
against the list of medically relevant genes, joined to the merged allele
table, and filtered into the target sets.

### 5.4 Panacus growth curves and node sharing

Panacus evaluates graph growth as assembly paths are
added, and `odgi heaps` provides an independent
size-related summary. Nodes are partitioned by
frequency among the selected input paths:
singletons are observed in one path, doubletons
in two, common nodes in 5–95% of paths and
core nodes in more than 95%. These categories
describe **assembly-path prevalence** and depend
on whether a panel has one or two haplotypes
per individual. Avoid describing two haplotypes
from one donor as two independent individuals.
A continued rise in singleton content suggests
that rare non-reference sequence remains
under-sampled; it does not imply that common
variant discovery is equally incomplete.

### 5.5 Non-reference bubbles, annotation and repeats

Non-reference nodes are identified relative to
the explicitly designated reference paths.
A node can be absent from GRCh38 while present
in CHM13, so report which reference is used.
Large structural bubbles are intersected with
GENCODE gene/exon annotations and the curated
medically relevant gene list. RepeatMasker
characterizes the associated sequence by
repeat class. Counts of bubbles, nodes, alleles
and physical base pairs are distinct quantities
and should not be pooled without a stated
unit and threshold. The manuscript compares
the frequency and content of repeat-rich
non-reference components between KPR and HPRC.

### Non-reference sequence

`wline_sample_nodes.py` reads the W lines of the GFA and writes the set
of nodes traversed by each sample. Node lengths come from the S lines, so the
total sequence a sample contributes to the graph is the sum of the lengths of
its nodes. Removing the nodes carried by the reference samples leaves the
non-reference sequence of each sample.

Pooling the per-sample lists gives, for every non-reference node, its length
and the number of samples carrying it. Nodes are then classified as core (more
than 95% of the samples), common (5–95%), doubletons (two samples) or
singletons (one sample), and the sequence of each class is written in FASTA
format.

The helper can also be run on its own:

```bash
python3 src/wline_sample_nodes.py results/graph/kpr/kpr.gfa nodes/kpr
```

---

## 6. Output

```
results/
├── graph/                    Minigraph–Cactus output (GFA, GBZ, VCF, odgi, viz)
├── index/                    vg, xg, r-index, haplotype index, GBWT, GFA v1, odgi, path lists
├── graph_statistics/         odgi stats, odgi heaps, allele counts, type counts
├── graph_variants/           deconstructed VCFs by reference, sample and haplotype
├── growth/                   Panacus growth curves and node frequency classes
├── allele_mapping/           per-assembly minigraph calls and the merged table
├── bubbles/                  bubbles, non-reference SV nodes, repeats, annotation
├── nonreference_sequence/    per-sample nodes, non-reference sequence, class FASTA
└── logs/
```

Main tables:

- `graph_statistics/GRAPH.odgi_stats.txt` — nodes, edges, paths and steps
- `graph_statistics/GRAPH.odgi_heaps.tsv` — pangenome size over 200 permutations
- `graph_statistics/GRAPH.allele_counts.txt` — variants and multi-allelic sites
- `graph_statistics/GRAPH.variant_typecounts.tsv` — singleton, doubleton, common
  and core counts for every call set
- `graph_statistics/GRAPH.nonreference_sequence_by_sharing.txt` — base pairs of
  non-reference sequence at each level of sharing
- `graph_statistics/GRAPH.nonreference_node_classes.tsv` — nodes and base pairs
  per frequency class
- `growth/GRAPH.growth.node.tsv`, `growth/GRAPH.node_classes.tsv`
- `bubbles/GRAPH/GRAPH.large.annotated.txt` — annotated large bubbles
- `bubbles/GRAPH/GRAPH.target.medical.txt`, `GRAPH.target.shared.txt`
- `nonreference_sequence/GRAPH/GRAPH.nonreference_sequence_per_sample.txt`
- `nonreference_sequence/GRAPH/classes/*.fa` — sequence of each frequency class

The indexed graph in `index/` and the GBZ in `graph/` are the input for
`6.pangenome_short_read_mapping.sh`. The HPRC+KPR graph also supplies the
AMY locus analysis described in
[AMY_locus_analysis_in_pagenomegraph.md](AMY_locus_analysis_in_pagenomegraph.md).

### Interpretation and audit points

The main checks for manuscript Figure 3 are the
relative graph complexities and content (Table S19);
per-sample SNP/indel/SV counts (Table S20), and KPR/HPRC
private and shared variants at locus and allele levels
(Figure 3b–d, Table S21); the accumulation of
non-reference sequence (Figure 3e, Table S22);
and the repeat composition of panel-specific
non-reference sequences (Figure 3f–g, Table S23).
Keep the exact sample/path counts, frequency
denominators and reference assignment with each table.
Variant records at the same coordinate with
different alleles may be counted separately.
Multiple nodes in one bubble are not necessarily
independent SV events. Graph representation and
variant discovery are dependent on alignment,
deconstruction, filtering and backbone choice.

---

## Study results

In our cohort, the combined HPRC+KPR graph contained
approximately 87.3 million nodes and 120.9 million edges,
with an estimated graph size near 3.2 Gb. Adding KOR
assemblies contributed approximately 1.6 million SNPs
and 0.2 million indels beyond HPRC alone, while the
non-reference sequence component grew from about
130 Mb in HPRC to about 152 Mb in the augmented graph.
These values describe the assembly panel and graph configuration used in our study.

---

## References

1. Kim JY, et al. *Constructing a Korean pangenome reveals distinct haplotypes at the amylase locus.* KPR manuscript, revision 2026; sections "Integrating Korean diversity into a pangenome framework" and "Assessment of the pangenome graph".
2. Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nature Biotechnology*. https://doi.org/10.1038/s41587-023-01793-w
3. Liao WW, et al. (2023). A draft human pangenome reference. *Nature*. https://doi.org/10.1038/s41586-023-05896-x
4. Guarracino A, et al. (2022). ODGI: understanding pangenome graphs. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btac308
5. Parmigiani L, et al. (2024). Panacus: fast and exact pangenome growth and core size estimation. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btae720

## Software and versions

The following versions are taken from the manuscript Methods and the
corresponding commands in the workflow.

| Operation | Software (manuscript version) |
|---|---|
| Graph construction | Minigraph–Cactus 2.9; Minigraph 0.19 in the construction pipeline |
| Graph alignment/allele extraction | Minigraph 0.21 for later allele calls |
| Graph indexing and variant deconstruction | vg 1.63.1 |
| Graph structure, paths and heaps | odgi 0.9.0 |
| Growth and core sequence | Panacus 0.4.1 |
| Bubble and graph extraction | gfatools 0.5; gfabase 0.5 |
| Repetitive sequence annotation | RepeatMasker 4.1.6 |
| Supporting utilities | bcftools (workflow specifies 1.16), bedtools, samtools, Python 3 |

The manuscript describes CHM13 as the backbone for its
pangenome comparison. The public wrapper also exposes
`REFERENCE` and `REF_PATHS` as configurable settings
(defaults are recorded in section 4); these must match
the intended reproducibility target and graph metadata.
