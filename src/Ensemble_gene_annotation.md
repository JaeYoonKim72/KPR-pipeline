# 2.Ensemble_gene_annotation.sh

The ensemble gene annotation workflow identifies protein-coding genes in
haplotype-resolved assemblies while reducing false-positive predictions.
It follows the hybrid gene-annotation strategy described for Asian
Pangenome Project phase 1 (APGp1), combining `de novo` prediction,
protein-to-genome alignment and reference annotation transfer.
The workflow consists of seven stages using RepeatMasker, AUGUSTUS,
Exonerate, Liftoff, BLASTp and HMM-Flagger classification data.

AUGUSTUS-derived candidates are first supported by human Swiss-Prot
proteins and compared with GENCODE v46 annotations transferred from
GRCh38. Genes not represented in the Liftoff annotation can be added
after protein homology screening. Only loci whose entire genomic spans
are classified as haploid by the corresponding assembly-specific
HMM-Flagger BED are retained. The separate copy-number gain (CNG)
tables summarize the resulting annotation across sample groups.

This script can execute the native annotation tools (`RUN_MODE=tools`)
or validate and integrate previously generated tool outputs
(`RUN_MODE=outputs`). It does not independently run HMM-Flagger;
the per-assembly classification BED is required as an input.
The workflow figure below corresponds to Figure SX1 in the manuscript.

![Figure SX1. Ensemble gene annotation workflow](../gene_annotation_workflow.png)

Part of the [KPR-pipeline](../README.md) repository.

---

## 1. Set up

### 1.1 Software

Use the following versions for reproducing the described workflow:

| Tool | Version | Role |
|---|---|---|
| RepeatMasker | 4.1.6 | Repetitive-sequence soft masking |
| AUGUSTUS | 3.5.0 | `de novo` gene prediction |
| Exonerate | 2.4.0 | Swiss-Prot protein-to-genome alignment |
| Liftoff | 1.6.3 | GENCODE annotation transfer and extra copies |
| BLASTp (BLAST+) | 2.15.0 | Protein homology search |
| HMM-Flagger | 1.0.0 | Precomputed assembly-reliability classifications |

The executable also requires Bash, Python ≥3.9, Biopython,
and (where applicable) RepeatMasker libraries and AUGUSTUS parameters.
`makeblastdb` must be available for native BLASTp execution.
The code reads supplied HMM-Flagger classification BEDs; it does not
invoke HMM-Flagger. The versions listed above are those in the
manuscript, rather than a claim that the script installs them.

### 1.2 Reference data

Prepare the GRCh38 reference FASTA, GRCh38-based GENCODE v46
annotation, and **human UniProtKB/Swiss-Prot release 2024_01**
protein sequences. These resources must correspond to the versions
used in the study. For GENCODE:

~~~bash
mkdir -p reference
wget -c \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz \
  -O reference/gencode.v46.annotation.gtf.gz
~~~

Reference paths are specified by `GENCODE_GTF`, `REFERENCE_FA` and
`SWISSPROT_FA`. A GENCODE-ID-compatible CHM13 annotation is optional
for the supplementary annotation-ratio table (`CHM13_ANNOTATION`).
An additional RepeatMasker satellite library (`SATELLITE_LIBRARY`)
and a Swiss-Prot protein/gene mapping table (`PROTEIN_GENE_MAP`)
are also optional.

### 1.3 Haplotype assemblies and manifest

The input TSV contains one row per **haplotype assembly** and must
include these columns:

| Column | Meaning |
|---|---|
| `sample` | Unique haplotype-assembly ID |
| `group` | `KPR`, `HPRC-EAS` or `HPRC-nEAS` |
| `fasta` | Assembly FASTA (plain or gzipped) |
| `flagger_bed` | Full classified BED for the same assembly; Hap/Dup/Col/Err/Unk labels |

Example (paths and sample identifiers are illustrative):

~~~tsv
sample	group	fasta	flagger_bed
KOR_EXAMPLE_1	KPR	assemblies/KOR_EXAMPLE_1.fa.gz	flagger/KOR_EXAMPLE_1.classified.bed
HPRC_EXAMPLE_1	HPRC-EAS	assemblies/HPRC_EXAMPLE_1.fa.gz	flagger/HPRC_EXAMPLE_1.classified.bed
~~~

All relative paths in the manifest are resolved against the
manifest's directory. Optional input columns are `masked_fasta`,
`augustus_gff`, `exonerate_gff`, `liftoff_gff`, `blastp_tsv`,
`contig_map`, `augustus_hints` and `extrinsic_cfg`.
`RUN_MODE=outputs` requires the five cached-result columns
`masked_fasta`, `augustus_gff`, `exonerate_gff`,
`liftoff_gff` and `blastp_tsv` to be populated for every row.
An empty, existing BLASTp TSV is allowed when no hits were found.

In the manuscript group comparison, there were 28 KPR,
10 HPRC-EAS and 84 HPRC-nEAS haplotype assemblies.
The actual manifest defines the sample counts for a new run;
do not substitute these counts as fabricated input records.

### 1.4 Environment

~~~bash
export REF_DIR=/path/to/reference
export GENCODE_GTF="$REF_DIR/gencode.v46.annotation.gtf.gz"
export REFERENCE_FA="$REF_DIR/GRCh38.fa"
export SWISSPROT_FA="$REF_DIR/uniprot_sprot.human.2024_01.fasta"
export SAMPLE_TABLE=/path/to/lists/samples.tsv
export OUT_DIR=/path/to/results/gene_annotation
export THREADS=8
export RUN_MODE=tools
~~~

---

## 2. Run

### 2.1 Check inputs and tools

~~~bash
bash 2.Ensemble_gene_annotation.sh --help
bash 2.Ensemble_gene_annotation.sh check
~~~

`check` validates the manifest and required files, and verifies
the tools needed for the selected execution mode. It does not
produce a finished annotation.

### 2.2 Run the full workflow

~~~bash
bash 2.Ensemble_gene_annotation.sh all
~~~

Stages can also run individually, after their preceding stages:

~~~bash
bash 2.Ensemble_gene_annotation.sh softmask
bash 2.Ensemble_gene_annotation.sh predict_orfs
bash 2.Ensemble_gene_annotation.sh protein_support
bash 2.Ensemble_gene_annotation.sh gene_annotation
bash 2.Ensemble_gene_annotation.sh compare_annotations
bash 2.Ensemble_gene_annotation.sh search_unannotated
bash 2.Ensemble_gene_annotation.sh filter_flagger
bash 2.Ensemble_gene_annotation.sh tables
~~~

Command-line aliases `mask`, `augustus`, `exonerate`,
`liftoff`, `compare`, `blastp` and `flagger` are supported.

### 2.3 Reuse existing tool outputs

~~~bash
RUN_MODE=outputs \
  SAMPLE_TABLE=/path/to/lists/samples_with_cached_outputs.tsv \
  OUT_DIR=/path/to/results/gene_annotation_cached \
  bash 2.Ensemble_gene_annotation.sh all
~~~

Reuse does not mean that reference files or HMM-Flagger BEDs
can be omitted. Each stage records input hashes and completed
outputs; a changed script, input or parameter requires a new
output directory. This portable workflow generates new analyses
and is not automatically identical to an earlier completed
supplementary-table integration.

---

## 3. Steps

| Step | Shell command | Main operation | Main software |
|---|---|---|---|
| 1 | `softmask` | Soft-mask repetitive sequence without shifting coordinates | RepeatMasker |
| 2 | `predict_orfs` | Predict protein-coding transcripts and exon/intron/CDS structures | AUGUSTUS |
| 3 | `protein_support` | Align Swiss-Prot proteins and retain AUGUSTUS models with amino acid identity ≥95% | Exonerate |
| 4 | `gene_annotation` | Transfer GENCODE v46 annotation; identify extra gene copies | Liftoff |
| 5 | `compare_annotations` | Require reciprocal CDS overlap ≥80% and intron concordance ≥80% | Python validation |
| 6 | `search_unannotated` | Review loci without Liftoff annotations; BLASTp identity ≥90% | BLASTp |
| 7 | `filter_flagger` | Keep gene spans fully covered by Hap regions | HMM-Flagger BED |
| Additional | `tables` | Generate annotation ratios and CNG summaries | Python |

---

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `THREADS` | `8` | Requested parallel threads |
| `REF_DIR` | `./reference` | Base reference directory |
| `LIST_DIR` | `./lists` | Default manifest directory |
| `SAMPLE_TABLE` | `$LIST_DIR/samples.tsv` | Haplotype manifest |
| `GENCODE_GTF` | `$REF_DIR/gencode.v46.annotation.gtf.gz` | GRCh38 GENCODE v46 |
| `SWISSPROT_FA` | `$REF_DIR/uniprot_sprot.human.fasta` | Human Swiss-Prot proteins; use release 2024_01 |
| `REFERENCE_FA` | `$REF_DIR/GRCh38.fa` if found | GRCh38 FASTA for native Liftoff |
| `OUT_DIR` | `./results/gene_annotation` | Output directory |
| `RUN_MODE` | `tools` | Run native tools or read cached outputs |
| `AUGUSTUS_SPECIES` | `human` | AUGUSTUS species model |
| `MIN_PROTEIN_IDENTITY` | `95` | AUGUSTUS vs Swiss-Prot protein identity (%) |
| `MIN_BLASTP_IDENTITY` | `90` | BLASTp identity for new gene annotations (%) |
| `MIN_EXONERATE_IDENTITY` | `90` | Additional Exonerate alignment identity (%) in source pipeline |
| `MIN_PROTEIN_COVERAGE` | `80` | Minimum protein alignment coverage (%) |
| `MIN_LIFTOFF_OVERLAP` | `80` | Reciprocal coding-interval support (%) |
| `MIN_EXON_SUPPORT` | `80` | Supported AUGUSTUS coding exons (%) |
| `MIN_INTRON_SUPPORT` | `80` | Supported AUGUSTUS introns (%) |
| `MIN_INTRON_MATCH` | `80` | AUGUSTUS–Liftoff intron concordance (%) |
| `MIN_UNANNOTATED_OVERLAP` | `95` | Strict Exonerate–AUGUSTUS CDS overlap for new loci (%) |
| `EXONERATE_BATCH_BP` | `4000000` | Approximate batch size for protein-to-genome alignment |
| `CHM13_ANNOTATION` | Unset | Optional GENCODE-ID-compatible CHM13 annotation |
| `SATELLITE_LIBRARY` | Unset | Optional additional repeat library |
| `PROTEIN_GENE_MAP` | Unset | Optional protein-to-gene mapping table |

The key identity thresholds are **95% for AUGUSTUS–Swiss-Prot
candidate validation** and **90% for BLASTp validation of genes
not overlapping Liftoff annotations**. The Exonerate command's
`--percent 80` controls alignment reporting and is distinct
from either amino acid identity cutoff.

---

## 5. Output

Each haplotype is processed separately under
`$OUT_DIR/samples/SAMPLE/`. Results are organized into
`mask/`, `augustus/`, `exonerate/`, `liftoff/`, `compare/`,
`blastp/` and `flagger/` folders.

| Output | Contents |
|---|---|
| `samples/SAMPLE/mask/` | Masked assembly FASTA and input audit |
| `samples/SAMPLE/augustus/orf_review.tsv` | Validated AUGUSTUS coding transcript candidates |
| `samples/SAMPLE/exonerate/protein_support_review.tsv` | Protein identity, coding coverage and alignment checks |
| `samples/SAMPLE/liftoff/liftoff.gff3` | GRCh38/GENCODE transferred annotation |
| `samples/SAMPLE/compare/liftoff_structure_review.tsv` | Reciprocal CDS overlap and intron concordance |
| `samples/SAMPLE/blastp/candidate_review.tsv` | UniProt-supported candidates not annotated by Liftoff |
| `samples/SAMPLE/flagger/locus_review.tsv` | Hap/Dup/Col/Err classification decisions |
| `samples/SAMPLE/flagger/final_annotation.gff3` | Final assembly-specific annotation |
| `tables/TableS13.annotation_ratios.tsv` | Optional CHM13-normalized GENCODE annotation ratios |
| `tables/TableS14.tsv` | Per-haplotype CNG counts |
| `tables/TableS15.tsv`, `tables/TableS15.two_headers.csv` | CNG frequency distributions |
| `tables/TableS16.all_genes.tsv` | Gene frequencies and between-group differences |
| `tables/TableS17.membership.tsv` | CNG membership among KPR/HPRC-EAS/HPRC-nEAS |
| `tables/metadata.json` | Method parameters and summary metadata |

The standard tables are written as TSV and CSV where applicable,
while `TableS15.two_headers.csv` and `TableS17.seven_columns.csv`
use layout-specific formats. The program also writes checkpoint
and input-provenance files.

---

## 6. Analysis procedures

### 6.1 Repeat soft masking

To reduce interference from repetitive sequence during gene annotation,
repeats in each assembly are soft-masked using **RepeatMasker v4.1.6**.
Soft masking preserves the original sequence coordinates. The soft-masked
assembly is used for coding transcript prediction.

### 6.2 AUGUSTUS coding transcript prediction

**AUGUSTUS v3.5.0** is used to identify coding transcript/ORF candidates.
Predicted CDSs are translated, and candidates with intact start and stop
codons and no internal stop codons are selected during coding-model validation.
The source prediction command and coordinate mapping are retained with each run.

The de novo prediction branch corresponds to this step. Existing extra-copy
validation runs that used Exonerate-derived AUGUSTUS hints retain that original
provenance; their predictions are not relabeled as independent de novo results.

### 6.3 Exonerate alignment and candidate selection

Human UniProtKB/Swiss-Prot proteins are aligned to candidate genomic regions
using **Exonerate v2.4.0**:

```text
--model protein2genome --percent 80
```

AUGUSTUS predictions are compared with the Exonerate alignments to retain
protein-supported candidates with amino acid identity >=95%. For the
Liftoff-unannotated candidate branch, the recorded postprocessing also requires
matching strand and >=95% CDS overlap in both directions between the AUGUSTUS
prediction and Exonerate alignment.

The native Exonerate command also specifies:

```text
--softmasktarget no
--showvulgar no
--showalignment no
--showquerygff no
--showtargetgff yes
```

The `--percent 80` reporting option is separate from the subsequent >=95%
identity filter. The recorded `--softmasktarget no` setting is preserved when
reproducing this run.

### 6.4 Liftoff gene annotation and extra-copy detection

**Liftoff v1.6.3** maps GRCh38-based GENCODE v46 annotations to each assembly
and detects additional gene copies using:

```text
-copies -sc 0.95 -polish -exclude_partial
```

The polished annotation and its gene/transcript hierarchy are retained.
Loci with explicit `extra_copy_number >= 1` are treated as Liftoff extra-copy
candidates. Different transcripts of the same gene locus are not counted as
independent copies.

### 6.5 Concordance with Liftoff annotation

The protein-supported AUGUSTUS transcripts are compared with the corresponding
Liftoff extra-copy transcripts. Extra copies are retained when both
**reciprocal CDS overlap and intron structure concordance are >=80%**.
The selected AUGUSTUS coding transcript must also satisfy >=95% amino acid
identity, intact start/stop codons and absence of internal stop codons.

The Python validation engine retains these additional screening parameters:

```text
--min-identity 90
--min-augustus-identity 95
--min-query-coverage 80
--min-pred-coverage 80
--min-structure-overlap 80
--min-exon-support 80
--min-intron-support 80
--min-intron-chain-match 80
```

Here, `--min-identity 90` belongs to the earlier Exonerate-support stage;
`--min-augustus-identity 95` controls the final AUGUSTUS protein identity.
The >=80% Liftoff concordance criteria are distinct from the >=95%
AUGUSTUS/Exonerate CDS-overlap criterion used for unannotated-region candidates.

### 6.6 BLASTp review of candidates without Liftoff annotation

Predicted proteins from candidates supported by AUGUSTUS and Exonerate but
lacking a corresponding Liftoff annotation are searched against human
UniProtKB/Swiss-Prot using **BLASTp v2.15.0**:

```text
-task blastp-fast
-comp_based_stats 0
-seg yes
-evalue 1e-5
-max_target_seqs 20
-max_hsps 1
-outfmt "6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send evalue bitscore stitle"
```

The BLASTp review applies the following thresholds (the additional coverage and ORF filters are retained from the source pipeline):

- Amino acid identity >=90% (the manuscript's BLASTp threshold).
- Query coverage >=95%.
- Reference coverage >=80%.
- Complete ORFs with no internal stop codons or ambiguous translations.

Candidates with >=50% lowercase CDS bases are excluded from the eligible set
as repeat-dominated predictions. The saved review also excludes the specified
LINE-1 ORF1 reference and quarantines the designated DUX-family and uncertain
REXO1L1P reference matches. Hits within 99% of the best bit score are checked
for conflicting gene symbols; ambiguous assignments are not automatically added.

Passing candidates receive homologous protein information and enter the
additional-annotation candidate set. All eligible loci in `all_loci.tsv` are
considered during integration. No minimum of three new gene copies is imposed by the current annotation integration.

### 6.7 HMM-Flagger filtering and final integration

The gene sets obtained from steps 5 and 6 are screened against the matching
assembly-specific **HMM-Flagger classification tracks**.

- Use the original Liftoff gene span for an existing extra copy.
- Use the projected AUGUSTUS gene span for a new candidate without Liftoff annotation.
- Exclude a locus if it overlaps a Dup, Col or Err region by at least 1 bp.
- Retain newly assessed loci only when their entire spans are covered by Haploid regions.
- Record unknown classification, incomplete BED coverage or unmatched contigs separately.

One-based inclusive annotation coordinates are converted to zero-based
half-open BED coordinates. Haplotype identity is checked when matching contigs,
and overlapping BED intervals are evaluated by covered-base union.

Before adding a new annotation, compare it with the full existing Liftoff gene
annotation. The current integration quarantines candidates overlapping an
existing Liftoff gene span rather than automatically adding duplicate or
conflicting loci. Overlapping same-symbol/same-strand candidate spans are
merged conservatively; adjacent nonoverlapping loci remain separate.

---

## 7. Copy-number gain summaries and annotation ratios

An added annotation does not automatically establish CNG. For newly assessed
candidates, the confirmed primary loci, existing validated extra loci and new
independent loci are considered together. A new CNG call is added only when
the confirmed locus count exceeds the GRCh38 GENCODE protein-coding gene-locus
count for that symbol. The previously established baseline CNG calls are
preserved.

For each gene, group-specific CNG frequency is the proportion of haplotype
assemblies retaining at least one extra copy. The denominators are 28 for KPR,
10 for HPRC-EAS and 84 for HPRC-nEAS. Group-specific gene lists describe CNG
evidence within these datasets, not absence of the gene from other populations.

Annotation ratios are calculated separately for protein-coding genes,
noncoding genes, protein-coding transcripts and noncoding transcripts by
dividing each assembly's unique GENCODE ID count by the corresponding CHM13
count. Swiss-Prot-supported additions do not automatically create GENCODE IDs
or change these ratios.


The CNG tables are generated after the HMM-Flagger filter.
A transferred transcript isoform is not treated as an independent
gene copy. The tables reflect the input haplotypes, and
the software does not assume that variants or copy-number patterns
observed in one panel are exclusive to its source population.

---

## 8. Provenance and reproducibility

`RUN_MODE=outputs` uses precomputed prediction/alignment files
but runs the same validation and output integration logic.
Each manifest entry is matched to its original haplotype sequence
and HMM-Flagger classification BED. Haplotype coordinate spaces
are not interchangeable. Where an HMM-Flagger BED is shared
across multiple assemblies without an explicit contig mapping,
the code does not treat unassigned segments as reliable Hap coverage.

The original source pipeline also includes additional screening
of translated ORFs, ambiguous Swiss-Prot gene assignments, repeats
and incomplete alignments. These extra filters are documented in
the detailed stages above and may make its results more restrictive
than the threshold summary in the main manuscript Methods.
Output files from the portable implementation should not be
presented as historical manuscript results without a corresponding
run and audit.

---

## Software and versions

RepeatMasker v4.1.6; AUGUSTUS v3.5.0; Exonerate v2.4.0;
Liftoff v1.6.3; BLASTp/BLAST+ v2.15.0;
precomputed HMM-Flagger v1.0.0 classifications;
Python 3.9+ and Biopython.

The annotation strategy follows APGp1. Source IDs in the
manuscript: APGp1 [R1], RepeatMasker [25], AUGUSTUS [R2],
Exonerate [R3], UniProtKB/Swiss-Prot [R4], GENCODE [R5],
Liftoff [59], BLASTp [R6], and HMM-Flagger [R7].
