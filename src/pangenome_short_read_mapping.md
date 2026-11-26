# 6.pangenome_short_read_mapping.sh

Short-read mapping, graph-based and linear-reference variant calling, and cross-reference comparisons (Figure 4; Tables S24–S27).

Part of the [KPR-pipeline](../README.md) repository.

## 1. Set up

Supply paired-end FASTQ files as SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz in DATA_DIR, with one sample ID per line in LIST_DIR/samples.list. Prepare the GRCh38 FASTA, bwa-mem2 index, dictionary and indexed known-sites VCFs for the GATK BQSR branch. Graph directories contain matching GRAPH.gbz, GRAPH.xg, GRAPH.hapl and GRAPH.snarls for kpr, hprc and hprc_kpr. Use vg v1.77.0 for giraffe mapping and the compatible indexes; the graph construction stage used vg v1.63.1. DeepVariant v1.6 runs through the configured container.

~~~bash
export DATA_DIR=/path/to/fastq
export REF_DIR=/path/to/reference
export GRAPH_DIR=/path/to/graphs
export LIST_DIR=/path/to/lists
export OUT_DIR=/path/to/results
export BQSR_KNOWN_SITES="/path/to/known_snps.vcf.gz /path/to/known_indels.vcf.gz"
~~~

For the leave-one-sample-out benchmark supply the sample-specific graph/index set, as described in section 7. Use the same reference coordinate system for all comparisons.

## 2. Run

~~~bash
bash 6.pangenome_short_read_mapping.sh all
~~~

For the study's graph BAM preprocessing, select rmdup rather than the GATK BQSR option:

~~~bash
GRAPH_BAM_MODE=rmdup bash 6.pangenome_short_read_mapping.sh bam_cleanup
GRAPH_BAM_MODE=rmdup bash 6.pangenome_short_read_mapping.sh small_variants
~~~

Run a single reference or sample by setting GRAPHS and SAMPLE_LIST. Supply the matching sample-specific graph directory for leave-one-out runs. Run grch38_baseline independently for linear-reference comparisons.

## 3. Steps

| Step | Description |
|---|---|
| `prepare_reference` | Reference paths, FASTA, dictionary, contig rename table and snarls, extracted from each graph |
| `read_prep` | Read name clean-up and per-sample read lists |
| `kmer_index` | kmc k-mer index used for haplotype sampling |
| `giraffe` | Read mapping to the pangenome graph, output in GAF |
| `surject` | Projection of the alignments onto the linear reference |
| `bam_cleanup` | Coordinate sorting, duplicate removal and indexing; optional GATK BQSR mode |
| `small_variants` | DeepVariant calling, PASS filtering, contig renaming |
| `sv_variants` | `vg pack` and `vg call` genotyping, PASS filtering |
| `classify` | Classification into SNP, bi-allelic SNP, indel and SV, genotype filtering and counting |
| `grch38_baseline` | Alignment, duplicate removal, GATK BQSR, DeepVariant, Manta and Delly against GRCh38 |
| `all` | All of the above, in order |

## 4. Configuration

| Variable | Default | Description |
|---|---|---|
| `DATA_DIR` | required | input FASTQ directory |
| `REF_DIR` | required | reference file directory |
| `GRAPH_DIR` | required | pangenome graph directory |
| `OUT_DIR` | `./results` | output directory |
| `LIST_DIR` | `./lists` | sample list directory |
| `TMP_DIR` | `$OUT_DIR/tmp` | temporary files |
| `GRAPHS` | `kpr hprc hprc_kpr` | references to evaluate |
| `REF_SAMPLE` | `GRCh38` | linear reference embedded in the graph |
| `KMC_K` | `29` | k-mer length |
| `KMC_MEM` | `128` | kmc memory in GB |
| `SV_MIN_LEN` | `50` | indels versus structural variants |
| `MIN_DP` | `7` | minimum genotype depth |
| `MIN_GQ` | `10` | minimum genotype quality |
| `MIN_QUAL` | `30` | minimum site quality |
| `DV_IMAGE` | Parabricks 4.4.0-1 | container image providing DeepVariant v1.6 |
| `BQSR_KNOWN_SITES` | required for the GATK BQSR branch | Space-separated indexed GRCh38 known-sites VCF(s), used by GATK BQSR |
| `BQSR_JAVA_MEM` | `16g` | Java heap for GATK BaseRecalibrator and ApplyBQSR |
| `GRCH38` | `$REF_DIR/GRCh38.fa` | linear reference for the baseline |
| `THREADS` | `8` | threads per job |
| `GRAPH_BAM_MODE` | `gatk_bqsr` | `rmdup` for samtools duplicate removal; `gatk_bqsr` for GATK removal and recalibration |

## 5. Output

```
results/
├── reads/                cleaned FASTQ files and read lists
├── kmer_index/           per-sample KFF files
├── reference/            per-graph reference FASTA, index, dictionary, paths
├── gaf/                  graph alignments
├── bam/                  surjected BAMs, duplicate metrics and BQSR BAMs
├── small_variants/       DeepVariant calls, PASS and NoZero sets
├── structural_variants/  vg pack and vg call output, PASS and NoZero sets
├── variant_classes/      classified and filtered call sets, count tables
├── grch38/               linear reference baseline (including BQSR BAM/table)
└── logs/
```

Main tables:

- `variant_classes/GRAPH.variant_counts.tsv` — number of records in every
  classified and filtered call set, per reference
- `small_variants/GRAPH/GRAPH.SAMPLE.PASS.vcf.gz` — small variants
- `structural_variants/GRAPH/GRAPH.SAMPLE.sv.PASS.vcf.gz` — structural variants
- `grch38/small_variants/`, `grch38/sv/` — the linear reference baseline

The final cross-reference variant-calling results for 50 Korean
short-read samples are reported in Supplementary Table S27. These
sample-level call sets are pipeline outputs; the formatted S27
spreadsheet is not automatically generated by this wrapper.

## 6. Graph-based short-read mapping

We used vg giraffe v1.77.0 with k-mer frequency-based haplotype sampling and minimizer-based mapping. The graph alignments were projected onto GRCh38 with vg surject, sorted and deduplicated using samtools rmdup. The rmdup mode retains the deduplicated BAM for DeepVariant v1.6; the alternative GATK mode performs duplicate removal and BQSR for general applications.

## 7. Leave-one-sample-out analysis

For the KOR seven-sample benchmark, remove both haplotype paths of the focal KOR individual from the KPR and HPRC+KPR GBWT using vg gbwt -R, and combine the modified GBWT with the original graph to construct a sample-specific GBZ. Generate the compatible haplotype and mapping indexes for this GBZ, run k-mer sampling and map with vg giraffe. Apply the same exclusion to HG002 and HG005 for HPRC and HPRC+KPR, including the corresponding trio analyses. The wrapper accepts the prepared graph directory; create the personalized GBZ and indexes before each sample run. Record the exact graph sample identifier and the retained paths.

## 8. Small and structural variant calling

DeepVariant v1.6 called SNPs and small indels from the projected BAMs. vg pack (-Q 5) computed read support from graph alignments, and vg call (-a -z) genotyped alleles present in the GBZ. Retain PASS calls and rename reference paths to GRCh38 chromosome names for the benchmark. The genotype-called VCFs retaining reference genotypes are used for trio analysis; non-reference call sets are used for count and accuracy comparisons.

## 9. GRCh38-based variant calling

Align SRS reads to GRCh38 and call small variants using DeepVariant v1.6. For the study's SV benchmark, retain precise Manta diploidSV INS/DUP calls and precise Delly DEL calls, with PASS and the size restrictions applied during benchmarking. The wrapper offers bwa-mem2, GATK duplicate removal/BQSR, DeepVariant, Manta and Delly for the linear baseline. Normalize DUPs to the insertion representation used by the benchmark before supplying the stage 5 insertion VCF.

## 10. Variant counts and sharing across references

For the independent 50 Korean SRS samples from the Korean Reference Genome Information project, use the four references GRCh38, KPR, HPRC and HPRC+KPR. Classify sample-level calls into SNP, 1–50-bp indel, SV INS and SV DEL; compute mean counts across the 50 individuals (Table S27; Figure 4f). Compare variant sharing among all four references after bringing the VCFs into the common GRCh38 coordinate system. Preserve sample IDs and the variant matching rules for each class when constructing the intersection counts for Figure 4g. The wrapper produces per-reference count tables and call sets; the cross-reference sharing table is a separate analysis of these outputs. This cohort has no corresponding truth sets and is used for variant-count and sharing comparisons, not recall estimation.

## References

- Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nature Biotechnology* 36:875–879. https://doi.org/10.1038/nbt.4227
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nature Biotechnology* 42:663–673. https://doi.org/10.1038/s41587-023-01793-w

## Software and versions

The variant-calling and benchmarking methods in the manuscript
used the versions below.

| Step | Software (manuscript version) |
|---|---|
| Graph short-read alignment and surjection | vg 1.77.0 (`giraffe` for study mapping; compatible `surject`, `pack`, `call`) |
| Linear reference alignment | bwa-mem2 2.2.1 |
| Small-variant calling | DeepVariant 1.6 |
| Linear-reference structural variants | Manta 1.6; Delly 0.9.1 |
| Duplicate removal and base-quality recalibration | GATK 4.3.0.0 (`MarkDuplicates`, `BaseRecalibrator`, `ApplyBQSR`) |
| Supporting read processing | samtools, bcftools, KMC; Docker/Singularity as required |

