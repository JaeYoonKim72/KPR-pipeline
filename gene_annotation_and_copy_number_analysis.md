# Gene annotation and copy-number gain analysis

This workflow combines repeat soft masking, coding transcript prediction,
protein alignment, GENCODE-based annotation, protein homology searches and
assembly-error filtering to reduce false positives in gene prediction and
copy-number gain (CNG) analysis.

Part of the [KPR-pipeline](README.md) repository.

The workflow uses six tools in seven steps, following the analysis layout
adapted from the APG ensemble annotation approach. KPR, HPRC-EAS and
HPRC-nEAS are compared at the haplotype-assembly level. T2T-CHM13 v2.0 is
used for the annotation-ratio comparison, and GRCh38 GENCODE v46 provides
the reference annotation for extra-copy analysis.

---

## 1. Set up

### 1.1 Software

Install the tools listed under [Software](#software) and make them available
on `PATH`. Python 3 is required for the filtering and table-generation scripts;
Biopython is required by the existing transcript-support validation engine.

### 1.2 Reference files

Prepare the following reference files:

- The GRCh38 FASTA used in the Liftoff analysis.
- GRCh38-based GENCODE v46 annotation in GTF/GFF format.
- The human UniProtKB/Swiss-Prot protein FASTA used by both Exonerate and BLASTp.
- Assembly-specific full classified HMM-Flagger BED files.
- The T2T-CHM13 v2.0 annotation used for annotation-ratio normalization.

The standard GENCODE v46 reference-chromosome annotation can be downloaded as:

```bash
mkdir -p reference
wget -c \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz \
  -O reference/gencode.v46.annotation.gtf.gz.download

gzip -t reference/gencode.v46.annotation.gtf.gz.download && \
  mv reference/gencode.v46.annotation.gtf.gz.download reference/gencode.v46.annotation.gtf.gz
```

Use the same reference scope as the original Liftoff analysis. Do not replace
GRCh38 reference annotation with a GTF already mapped to CHM13 or a KPR assembly.
Retain the original human Swiss-Prot FASTA and its checksum when reproducing
protein searches.

### 1.3 Assemblies and sample groups

| Group | Haplotype assemblies |
|---|---:|
| KPR | 28 |
| HPRC-EAS | 10 |
| HPRC-nEAS | 84 |

Each haplotype is processed separately. Preserve sample, haplotype, contig and
coordinate information throughout the workflow. Assembly coordinates from
different haplotypes must not be combined by numeric position alone.

---

## 2. Seven-step workflow

### 1) Repeat soft masking

To reduce interference from repetitive sequence during gene annotation,
repeats in each assembly are soft-masked using **RepeatMasker v4.1.6**.
Soft masking preserves the original sequence coordinates. The soft-masked
assembly is used for coding transcript prediction.

### 2) AUGUSTUS coding transcript prediction

**AUGUSTUS v3.5.0** is used to identify coding transcript/ORF candidates.
Predicted CDSs are translated, and candidates with intact start and stop
codons and no internal stop codons are selected during coding-model validation.
The source prediction command and coordinate mapping are retained with each run.

The de novo prediction branch corresponds to this step. Existing extra-copy
validation runs that used Exonerate-derived AUGUSTUS hints retain that original
provenance; their predictions are not relabeled as independent de novo results.

### 3) Exonerate alignment and candidate selection

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

The supplied unannotated-region Exonerate run used these additional options:

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

### 4) Liftoff gene annotation and extra-copy detection

**Liftoff v1.6.3** maps GRCh38-based GENCODE v46 annotations to each assembly
and detects additional gene copies using:

```text
-copies -sc 0.95 -polish -exclude_partial
```

The polished annotation and its gene/transcript hierarchy are retained.
Loci with explicit `extra_copy_number >= 1` are treated as Liftoff extra-copy
candidates. Different transcripts of the same gene locus are not counted as
independent copies.

### 5) Concordance with Liftoff annotation

The protein-supported AUGUSTUS transcripts are compared with the corresponding
Liftoff extra-copy transcripts. Extra copies are retained when both
**reciprocal CDS overlap and intron structure concordance are >=80%**.
The selected AUGUSTUS coding transcript must also satisfy >=95% amino acid
identity, intact start/stop codons and absence of internal stop codons.

For the existing identity95 validation run, the recorded settings are:

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

### 6) BLASTp review of candidates without Liftoff annotation

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

The completed review uses:

- Amino acid identity >=95%.
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
considered during integration. The historical exploratory `min_copies=3`
setting does not impose a three-new-copy requirement on this integration.

### 7) HMM-Flagger filtering and final integration

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

## 3. CNG classification and annotation ratios

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

---

## 4. Reuse completed results

The completed-output integration is run with the previously supplied
`CNG_BLASTp_Flagger_update` package. It reads existing tool outputs rather than
rerunning the native programs:

```bash
bash CNG_BLASTp_Flagger_update/run.sh \
  --inputs /path/to/matching_input_snapshot \
  --reference-gtf /path/to/reference/gencode.v46.annotation.gtf.gz \
  --out-parent /path/to/results
```

The input snapshot records source files, sample groups, candidate-review paths
and Flagger checksums. Missing source files or changed BED checksums stop the
integration. A new output directory is created without overwriting the baseline.

The supplied update includes new-candidate reviews for 113 haplotypes:
28 KPR, 4 HPRC-EAS and 81 HPRC-nEAS, comprising 434 eligible independent loci
before Flagger filtering. Baseline calls for the nine assemblies without a new
review remain unchanged. Coverage is recorded in `sample_coverage.csv`; these
updated tables do not represent complete supplemental screening of all 122
assemblies.

---

## 5. Output

| File | Description |
|---|---|
| `TableS14.csv` | Per-haplotype CNG-gene counts before and after integration, including added genes |
| `TableS14_group_summary.csv` | Group union, mean and range summaries |
| `TableS15.two_headers.csv` | Gene counts by the number of CNG-positive haplotypes |
| `TableS15.named_protein.two_headers.csv` | Named protein-coding version of the frequency table |
| `TableS16.same_selected_genes.csv` | Updated group frequencies and differences for the original selected genes |
| `TableS16.all_genes.csv` | Frequencies and differences for all retained CNG genes |
| `TableS17.venn_counts.csv` | Seven mutually exclusive CNG-membership counts before and after integration |
| `TableS17.membership.csv` | Per-gene group membership and contributing haplotypes |
| `TableS17.KPR_only.csv` | CNG genes observed only in KPR among the analysed groups |
| `new_locus_flagger_review.csv` | New-locus Flagger and existing-annotation overlap decisions |
| `new_CNG_classification.csv` | Evidence distinguishing annotation additions from new CNG calls |
| `sample_coverage.csv` | Availability of the additional candidate review for each haplotype |
| `input_audit.csv`, `metadata.json` | Input provenance, checksums and filtering criteria |

TSV copies are also generated. CSV files can be opened in Excel; an XLSX
workbook is generated when `openpyxl` is available.

---

## Software

RepeatMasker 4.1.6 · AUGUSTUS 3.5.0 · Exonerate 2.4.0 · Liftoff 1.6.3 ·
BLAST+ 2.15.0 · HMM-Flagger classification results · Python 3 · Biopython
