# Software and workflow registration metadata

This document prepares KPR-pipeline v1.0.0 for scholarly registries. Identifiers are not issued by this repository. Confirm authorship, ORCID information and licenses before publishing a permanent DOI.

## Shared metadata

| Field | Value |
|---|---|
| Name | KPR-pipeline |
| Title | KPR-pipeline: Korean Pangenome Reference construction and analysis |
| Version | 1.0.0 |
| Homepage and source | https://github.com/JaeYoonKim72/KPR-pipeline |
| GitHub release tag | v1.0.0 (create after enabling Zenodo integration) |
| License | MIT |
| Repository contact | Jae-Yoon Kim |
| Affiliation | Genomic Medicine Research Center, Korea Research Institute of Bioscience and Biotechnology (KRIBB), Republic of Korea |
| Operating system | Linux |
| Programming languages | Bash, Python, R |
| Resource type | Research software; command-line computational workflows |

**Description:** KPR-pipeline contains eight bioinformatics workflows for long-read genome assembly and QC, ensemble gene annotation and copy-number analysis, pangenome graph construction and characterization, pangenome variant comparisons, KOR/GIAB variant benchmarking and trio genotype inheritance inconsistency analysis, graph-based short-read mapping and variant calling, and AMY-locus structural and population-genetic analyses. Documentation describes inputs, methods, outputs and implementation limitations. Genomic data, pangenome graph binaries and reference resources are distributed separately.

**Keywords:** Korean pangenome, pangenomics, genome assembly, long-read sequencing, structural variation, graph genome, variant benchmarking, GIAB, AMY locus, population genetics.

## Zenodo (software DOI)

1. Sign in at https://zenodo.org/ and link GitHub in **Profile → Linked accounts**.
2. From the Zenodo GitHub page, click **Sync now** and enable `JaeYoonKim72/KPR-pipeline`.
3. Verify the `CITATION.cff` creator attribution. It intentionally does **not** contain a DOI or release date before release. No `.zenodo.json` is used because CITATION.cff alone is sufficient.
4. **After enabling the repository in Zenodo**, run the manually triggered GitHub Action **Publish KPR-pipeline v1.0.0 release** on the **main** branch. Alternatively create the GitHub Release manually with tag `v1.0.0`, source `main`, and notes in `docs/RELEASE_NOTES_v1.0.0.md`.
5. Confirm the released GitHub tag and verify the Zenodo archive and actual DOI. Keep the version-specific DOI for the manuscript and do not invent an identifier.

GitHub release URL (once published): https://github.com/JaeYoonKim72/KPR-pipeline/releases/tag/v1.0.0

Zenodo instructions: https://help.zenodo.org/docs/github/enable-repository/ and https://help.zenodo.org/docs/github/archive-software/github-upload/.

## WorkflowHub (citable workflows)

Open https://workflowhub.eu/ and import Git repository `https://github.com/JaeYoonKim72/KPR-pipeline`. Choose release tag `v1.0.0`. For one overarching workflow-suite record, the primary file can be `README.md` (WorkflowHub accepts it for a series of Bash scripts). Suggested title: **KPR-pipeline: Korean Pangenome Reference analysis workflows**. Use the shared Description, MIT license, appropriate maturity and verified creator(s). Make it public, freeze the version, and mint a DOI.

The same repository and tag can support separate workflow entries with individual executable files:

| Suggested workflow title | Main executable |
|---|---|
| KPR genome assembly and quality assessment | `1.genome_assembly_and_statistics.sh` |
| KPR ensemble gene annotation | `2.Ensemble_gene_annotation.sh` |
| KPR pangenome graph construction | `3.pangenome_construction_and_statistics.sh` |
| KPR variant-count comparison | `4.compare_pangenome_variant_counts.sh` |
| KPR variant benchmarking and trio inheritance assessment | `5.Variant_benchmark.sh` |
| KPR short-read graph mapping and variant calling | `6.pangenome_short_read_mapping.sh` |
| KPR AMY locus structural haplotype analysis | `7.AMY_locus_analysis_in_pagenomegraph.sh` |
| KPR AMY flanking-region population genetics | `8.Population_genetics_in_AMY_flanking_locus.sh` |

WorkflowHub documentation: https://about.workflowhub.eu/docs/registering_workflows/adding-files/ and https://about.workflowhub.eu/docs/citable/.

## bio.tools

Search https://bio.tools/ for duplicates, then use https://bio.tools/register to create one combined software entry. Name: **KPR-pipeline**; homepage: GitHub repository URL; version: **1.0.0**; type: command-line pipeline/workflows; license: **MIT**; description: shared Description; documentation: `README.md`; source code: repository URL; creator/contact: Jae-Yoon Kim (verify authorship and optional ORCID). Function descriptions include genome assembly, gene annotation, pangenome construction, variant calling/benchmarking and population genetics. Record only the **actual issued biotoolsID**.

Instructions: https://biotools.readthedocs.io/en/latest/quickstart_guide.html.

## SciCrunch RRID

At https://scicrunch.org/create/resource, search for an existing **KPR-pipeline** entry. If absent, submit a software-resource suggestion with the shared Description, homepage, institution and MIT license. A suggestion does **not** immediately confer an RRID; a SciCrunch curator must accept it. Record only the resulting identifier.

Instructions: https://scicrunch.org/scicrunch/about/Curation%20Guidelines.

## Manuscript

The GitHub repository and MIT license are public. Zenodo DOI, WorkflowHub DOI, bio.tools ID and SciCrunch RRID must remain **unreported as issued** until their actual registry records are available. Keep the Zenodo software-archive DOI and the WorkflowHub workflow DOI separate.
