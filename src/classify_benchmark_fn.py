#!/usr/bin/env python3
"""Classify benchmark false negatives; retain unresolved cases when graph evidence is unavailable."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path

from summarize_mendelian import read_vcf, variants_key


def classify_small(path):
    for r in read_vcf(path):
        truth = r["genotypes"].get("TRUTH", {})
        if truth.get("BD") != "FN":
            continue
        bk = truth.get("BK", ".")
        reason = ("detection_failure" if bk in (".", "") else
                  "matching_failure" if bk in ("lm", "am") else
                  "unresolved")
        yield r, reason, bk


def nearby_query(variant, calls, distance=500):
    return any(q["chrom"] == variant["chrom"] and abs(q["pos"] - variant["pos"]) <= distance
               for q in calls)


def classify_sv(fn_path, tpb_path, query_path, dataset, evidence=None):
    queries = list(read_vcf(query_path))
    represented = None
    if evidence:
        represented = {variants_key(r) for r in read_vcf(evidence)}
    for r in read_vcf(tpb_path):
        mismatch = r["info"].get("GTMatch")
        if mismatch is not None and str(mismatch) != "0":
            yield r, "matching_failure", "genotype_discordance"
    for r in read_vcf(fn_path):
        if not nearby_query(r, queries):
            yield r, "detection_failure", "no_query_within_500bp"
        elif dataset in ("kpr", "hprc_kpr", "hprc"):
            if represented is None:
                yield r, "unresolved", "personalized_GBZ_not_tested"
            elif variants_key(r) not in represented:
                yield r, "detection_failure", "truth_allele_absent_from_GBZ_evidence"
            else:
                yield r, "matching_failure", "query_allele_or_genotype_discordance"
        else:
            yield r, "matching_failure", "query_allele_or_genotype_discordance"


def run(outdir, samples, datasets, levels, gbz_evidence_dir):
    root = Path(outdir)
    lines = []
    for sample in samples:
        for ds in list(datasets) + ["deepvariant"]:
            for variant_type in ("SNP", "INS", "DEL"):
                path = root / "bench" / "small" / sample / f"{ds}.{variant_type}.vcf.gz"
                if not path.is_file():
                    continue
                for rec, status, evidence in classify_small(path):
                    lines.append((sample, ds, "small", "-", variant_type,
                                  rec["chrom"], rec["pos"], rec["id"], status, evidence))
        for dataset in list(datasets) + ["manta", "delly"]:
            for level in levels:
                for typ in ("INS", "DEL"):
                    branch = root / "bench" / "sv" / sample / dataset / level / typ
                    if not (branch / "fn.vcf.gz").is_file():
                        continue
                    qpath = root / "query" / "sv" / sample / f"{dataset}.{typ}.vcf.gz"
                    evidence = None
                    if gbz_evidence_dir and dataset in datasets:
                        candidate = Path(gbz_evidence_dir) / sample / f"{dataset}.{typ}" / "tp-base.vcf.gz"
                        if candidate.is_file():
                            evidence = candidate
                    for rec, status, note in classify_sv(
                            branch / "fn.vcf.gz", branch / "tp-base.vcf.gz",
                            qpath, dataset, evidence):
                        lines.append((sample, dataset, "SV", level, typ,
                                      rec["chrom"], rec["pos"], rec["id"], status, note))
    path = root / "summary" / "fn_failure_classes.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(("sample", "dataset", "variant_group", "support_level",
                         "type", "chrom", "pos", "id", "classification", "evidence"))
        writer.writerows(lines)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", required=True)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--datasets", nargs="+", default=["kpr", "hprc_kpr", "hprc"])
    p.add_argument("--support-levels", nargs="+", default=["ge1", "ge2", "ge3"])
    p.add_argument("--gbz-evidence-dir", default="")
    a = p.parse_args()
    print(run(a.outdir, a.samples, a.datasets, a.support_levels,
              a.gbz_evidence_dir))


if __name__ == "__main__":
    main()
