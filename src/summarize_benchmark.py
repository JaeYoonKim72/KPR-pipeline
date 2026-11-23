#!/usr/bin/env python3
"""Summarize genotype-aware hap.py and Truvari benchmarks across samples."""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import statistics


METRICS = ("Precision", "Recall", "F1", "Presence", "GT_concordance")


def ratio(num, den):
    return num / den if den else None


def f1_score(precision, recall):
    if precision is None or recall is None:
        return None
    return ratio(2 * precision * recall, precision + recall) if precision + recall else 0.0


def vcf_lines(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def hap_annotated(path):
    if not path.is_file():
        raise FileNotFoundError(f"hap.py annotated VCF not found: {path}. Run hap.py with -V")
    counts = {"am": 0, "lm": 0, "nocall": 0, "truth_fn": 0}
    with vcf_lines(path) as inp:
        sample_idx = None
        for line in inp:
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if "TRUTH" not in header:
                    raise ValueError(f"TRUTH sample missing from {path}")
                sample_idx = header.index("TRUTH")
                continue
            if line.startswith("#"):
                continue
            if sample_idx is None:
                raise ValueError(f"#CHROM header missing: {path}")
            fields = line.rstrip("\n").split("\t")
            fmt = fields[8].split(":")
            val = fields[sample_idx].split(":")
            data = dict(zip(fmt, val))
            if data.get("BK") == "am":
                counts["am"] += 1
            if data.get("BD") != "FN":
                continue
            counts["truth_fn"] += 1
            bk = data.get("BK", ".")
            if bk == "lm":
                counts["lm"] += 1
            elif bk == "am":
                pass
            elif bk in (".", ""):
                counts["nocall"] += 1
    return counts


def hap_row(path, typ):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    kind = "SNP" if typ == "SNP" else "INDEL"
    selected = [r for r in rows if r.get("Type") == kind and r.get("Filter") == "ALL"]
    if len(selected) != 1:
        raise ValueError(f"Expected one Type={kind}, Filter=ALL row in {path}")
    r = selected[0]
    qtotal, ttotal = int(r["QUERY.TOTAL"]), int(r["TRUTH.TOTAL"])
    ann = hap_annotated(path.with_name(path.name.replace(".summary.csv", ".vcf.gz")))
    presence_numerator = int(r["TRUTH.TP"]) + ann["am"]
    if ttotal and presence_numerator > ttotal:
        raise ValueError("hap.py presence count exceeds truth total; check BK annotations")
    presence = ratio(presence_numerator, ttotal)
    gt = ratio(int(r["TRUTH.TP"]), presence_numerator)
    return {
        "Truth": ttotal, "Query": qtotal, "TP": int(r["TRUTH.TP"]),
        "FN": int(r["TRUTH.FN"]), "FP": int(r["QUERY.FP"]),
        "Precision": float(r["METRIC.Precision"]),
        "Recall": float(r["METRIC.Recall"]),
        "F1": float(r["METRIC.F1_Score"]),
        "Presence": presence, "GT_concordance": gt,
        "FN_detection_no_query": ann["nocall"],
        "FN_matching_allele": ann["lm"],
        "FN_matching_genotype": ann["am"],
    }


def truvari_row(path):
    with path.open(encoding="utf-8") as stream:
        d = json.load(stream)
    for k in ("base cnt", "comp cnt", "TP-base", "TP-comp", "TP-base_TP-gt",
              "TP-comp_TP-gt", "TP-base_FP-gt", "TP-comp_FP-gt", "FP", "FN"):
        if k not in d:
            raise ValueError(f"Missing Truvari summary field {k} in {path}")
    truth, query = d["base cnt"], d["comp cnt"]
    tp_truth, tp_query = d["TP-base_TP-gt"], d["TP-comp_TP-gt"]
    recall = ratio(tp_truth, truth)
    precision = ratio(tp_query, query)
    return {
        "Truth": truth, "Query": query, "TP": tp_truth,
        "FN": d["FN"] + d["TP-base_FP-gt"],
        "FP": d["FP"] + d["TP-comp_FP-gt"],
        "Precision": precision, "Recall": recall,
        "F1": f1_score(precision, recall),
        "Presence": ratio(d["TP-base"], truth),
        "GT_concordance": d.get("gt_concordance"),
        "Matched_alleles": d["TP-base"],
        "GT_discordant_truth": d["TP-base_FP-gt"],
        "GT_discordant_query": d["TP-comp_FP-gt"],
    }


def write(path, rows):
    if not rows:
        return
    fields = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output:
        w = csv.DictWriter(output, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else
                            ("NA" if v is None else v)) for k, v in row.items()})


def means(rows, keys):
    groups = {}
    for row in rows:
        groups.setdefault(tuple(row[k] for k in keys), []).append(row)
    out = []
    for label, samples in sorted(groups.items()):
        entry = dict(zip(keys, label), N=len(samples))
        for key in ("Truth", "Query", "TP", "FN", "FP"):
            entry[key + "_mean"] = statistics.mean(r[key] for r in samples)
        for metric in METRICS:
            valid = [r[metric] for r in samples if r[metric] is not None]
            entry[metric + "_mean"] = statistics.mean(valid) if valid else None
            entry[metric + "_sd"] = (statistics.stdev(valid) if len(valid) > 1 else
                                      0.0 if valid else None)
        out.append(entry)
    return out


def run(outdir, samples, support_levels):
    bench, dest = Path(outdir) / "bench", Path(outdir) / "summary"
    small = []
    sv = []
    for sample in samples:
        for path in sorted((bench / "small" / sample).glob("*.summary.csv")):
            dataset, typ = path.name.removesuffix(".summary.csv").rsplit(".", 1)
            if typ not in ("SNP", "DEL", "INS"):
                continue
            small.append(dict(Sample=sample, Dataset=dataset, Type=typ, **hap_row(path, typ)))
        for path in sorted((bench / "sv" / sample).glob("*/*/*/summary.json")):
            typ, level, dataset = path.parent.name, path.parent.parent.name, path.parent.parent.parent.name
            if level not in support_levels:
                continue
            sv.append(dict(Sample=sample, Dataset=dataset, Support=level, Type=typ, **truvari_row(path)))
    if not small and not sv:
        raise ValueError(f"No benchmark summaries under {bench}")
    write(dest / "kor_small_detail.tsv", small)
    write(dest / "kor_small_mean.tsv", means(small, ("Type", "Dataset")))
    write(dest / "kor_sv_detail.tsv", sv)
    write(dest / "kor_sv_mean.tsv", means(sv, ("Support", "Type", "Dataset")))
    return len(small), len(sv)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--graphs", nargs="*", default=[])
    ap.add_argument("--support-levels", nargs="+", default=["ge1", "ge2", "ge3"])
    args = ap.parse_args()
    counts = run(args.outdir, args.samples, args.support_levels)
    print(f"Summarized {counts[0]} small-variant and {counts[1]} SV benchmark results")


if __name__ == "__main__":
    main()
