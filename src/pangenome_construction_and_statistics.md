# 3.pangenome_construction_and_statistics.sh

Construction and assessment of KPR and HPRC+KPR pangenome graphs (Figures 3a–b, 3e–h, S7; Tables S19–S20, S22–S23).

Part of the [KPR-pipeline](../README.md) repository.

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


## 3. Steps

| Step | Description |
|---|---|
| `construct` | Minigraph–Cactus graph construction |
| `index_graph` | vg conversion and indexing (vg, xg, r-index, haplotype index, GBWT, path lists), odgi build and sort, node FASTA |
| `graph_stats` | Nodes, edges, paths and steps from `odgi stats`; pangenome size and saturation curve from `odgi heaps` (200 permutations) |
| `variants` | `vg deconstruct` against the reference paths and against every sample and haplotype, size and type splitting, allele counts, variant classification |
| `growth` | Panacus growth curves and node frequency classes, for the configured graphs and supplied HPRC graph |
| `allele_mapping` | `minigraph --call` of each haplotype assembly against the graph, merged into a per-bubble allele table |
| `nonref_nodes` | Bubbles from `gfatools bubble`, non-reference SV nodes, frequency classes and repeat content |
| `bubble_annotation` | Gene, exon and medically relevant gene annotation of the large bubbles, and target selection |
| `nonref_sequence` | Sequence carried by each sample, the non-reference part of it, and the sequence of each frequency class |
| `all` | All of the above, in order |

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

## 5. Output

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

## 6. Pangenome graph construction

We constructed KPR from 28 haplotype assemblies of 14 KOR individuals and HPRC+KPR from these assemblies and 94 HPRC v1.0 haplotype assemblies. GRCh38 and CHM13 reference paths were not additional population samples. Assemblies were ordered by population (AFR, EUR, AMR, EAS, KOR, as applicable). Minigraph–Cactus v2.9 used Minigraph v0.19 with CHM13 as the backbone, soft-masked repetitive regions, excluded MAPQ <5 and soft-masked segments >100 kb, refined chromosome-level graphs by base-level alignment, normalized and merged them with the human-data defaults. The graphs were indexed using vg v1.63.1.

For study-specific construction, set REFERENCE to the actual CHM13 backbone identifier in the seqfile and retain the GRCh38 reference path for stage 4. The configured reference identifiers and seqfile ordering determine the graph; the wrapper default remains GRCh38 for other applications. The study used 27 CPU cores on a node with 1.5 TB physical RAM. Cactus memory and disk arguments are configured separately.

~~~bash
THREADS=27 REFERENCE="CHM13 GRCh38" REF_PATHS="CHM13 GRCh38" \
  bash 3.pangenome_construction_and_statistics.sh construct
~~~

Use the exact reference labels from the seqfile when running this command.

## 7. Graph statistics and variant identification

We compared KPR, the published HPRC graph and HPRC+KPR. odgi v0.9.0 stat (odgi stats) provided graph topology and sequence statistics; odgi heaps used 200 permutations (Table S19). vg deconstruct generated variant records against CHM13 for Figure 3b and the GRCh38-referenced input of stage 4. Small variants and SVs were separated by the 50-bp REF–ALT length-difference threshold and subdivided by type and occurrence. Per-individual variant counts in HPRC+KPR were summarized by population, keeping singleton, doubleton and polymorphic categories (Figure 3e; Table S20). The locus and allele comparison in Figure 3c–d uses stage 4 rather than subtracting counts from separate graphs.

## 8. Pangenome growth curves

We used Panacus v0.4.1 to calculate cumulative non-reference sequence growth at the sample level, grouping the two haplotypes of each individual. The main growth analysis used the augmented HPRC+KPR graph; Figure S7 used HPRC (47 individuals) and KPR (14 individuals) subsets within that same graph. The combined panel comprised 61 individuals and 122 population haplotype assemblies. Exclude reference-only paths from the population grouping. Use the ordered population path lists for the full and subset runs, preserving the same graph, reference-node definition and frequency classification.

Nodes were classified as core (>95%), common (5–95%), doubletons (two assemblies) or singletons (one assembly), with the study's frequency denominator and grouping applied consistently. The growth wrapper runs the configured graph files independently; for the subset analysis, supply the HPRC+KPR graph and the corresponding path-selection/grouping input to Panacus for each panel. Record the actual grouping and frequency denominators when preparing S22. The wrapper retains a rare category for nodes below 5% with more than two carriers; the original study used its four reported categories.

## 9. Non-reference sequence and repeat analysis

Extract nodes traversed by each sample and identify non-reference nodes relative to the reference paths used in the study. Within HPRC+KPR, construct the union of HPRC non-reference nodes and the union of KPR non-reference nodes. HPRC-private nodes are the HPRC union minus the KPR union, and KPR-private nodes the reverse. Calculate cumulative non-reference bases from node lengths without counting a node twice. Use all panel members, including members with no non-reference node in a given class, for frequency denominators.

Annotate non-reference sequences by repeat type and frequency class, following the repeat procedure in stage 1. Summarize HPRC, KPR, HPRC-private and KPR-private repeat counts and annotated lengths separately (Figures 3g–h; Table S23). The nonref_nodes and nonref_sequence steps provide node, sequence and repeat-annotation inputs; the panel-specific unions, private subsets and final repeat summaries are calculated within the augmented graph.

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

