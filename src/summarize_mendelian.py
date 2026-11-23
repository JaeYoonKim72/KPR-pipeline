#!/usr/bin/env python3
"""Summarize RTG Mendelian calls and truth-concordance among child variants."""

import argparse
from collections import Counter
import csv
import gzip
from pathlib import Path


CLASSES = ("SNP", "Indel", "INS", "DEL")


def read_vcf(path):
    if not Path(path).is_file():
        return
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8") as handle:
        samples = []
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError(f"Malformed VCF record in {path}")
            fmt = f[8].split(":") if len(f) > 8 else []
            genotypes = dict(zip(samples, (dict(zip(fmt, col.split(":"))) for col in f[9:])))
            info = {}
            for item in f[7].split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    info[key] = val
                elif item:
                    info[item] = True
            yield {
                "chrom": f[0], "pos": int(f[1]), "id": f[2], "ref": f[3],
                "alt": f[4], "info": info, "genotypes": genotypes,
            }


def gt(rec, person):
    raw = rec["genotypes"].get(person, {}).get("GT", ".")
    if "|" in raw:
        return raw.split("|")
    return raw.split("/")


def diploid_called(gt_value):
    return len(gt_value) == 2 and all(x not in (".", "") for x in gt_value)


def child_nonref(rec, person):
    value = gt(rec, person)
    return diploid_called(value) and any(int(x) > 0 for x in value)


def variants_key(rec):
    return (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])


def read_bed(path):
    if not path or not Path(path).is_file():
        return {}
    op = gzip.open if str(path).endswith(".gz") else open
    records = {}
    with op(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, stop, *_ = line.split()
            records.setdefault(chrom, []).append((int(start), int(stop)))
    for chrom, spans in records.items():
        spans.sort()
        merged = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        records[chrom] = merged
    return records


def confidence(rec, bed, sv):
    if not bed:
        return False
    start = rec["pos"] - 1
    end = start + len(rec["ref"])
    if sv:
        value = rec["info"].get("END")
        if isinstance(value, str) and value.isdigit():
            end = max(end, int(value))
    return any(start >= a and end <= b for a, b in bed.get(rec["chrom"], ()))


def happy_verdict(path):
    outcomes = {}
    for record in read_vcf(path):
        sample = record["genotypes"].get("QUERY", {})
        bd = sample.get("BD")
        if bd not in ("TP", "FP", "FN", "N"):
            continue
        status = "true" if bd == "TP" else "false" if bd == "FP" else "unknown"
        outcomes[variants_key(record)] = status
        if record["id"] not in (".", ""):
            outcomes[("id", record["id"])] = status
    return outcomes


def truvari_verdict(work):
    result = {}
    for file, source in (("tp-comp.vcf.gz", "match"), ("fp.vcf.gz", "false")):
        for r in read_vcf(work / file):
            if source == "match":
                value = r["info"].get("GTMatch")
                if value is None:
                    raise ValueError("Truvari tp-comp record lacks INFO/GTMatch")
                status = "true" if str(value) == "0" else "false"
            else:
                status = "false"
            result[variants_key(r)] = status
            if r["id"] not in (".", ""):
                result[("id", r["id"])] = status
    return result


def status_for(rec, outcomes):
    if ("id", rec["id"]) in outcomes and rec["id"] not in (".", ""):
        return outcomes[("id", rec["id"])]
    return outcomes.get(variants_key(rec), "unknown")


def calculate(directory, reference, trio, child, father, mother, small_bed, sv_bed):
    root = Path(directory)
    beds = {"small": read_bed(small_bed), "sv": read_bed(sv_bed)}
    results = []
    per_variant = []
    for cls in CLASSES:
        sv = cls in ("INS", "DEL")
        raw = list(read_vcf(root / f"{cls}.annot.vcf.gz"))
        observed = Counter()
        valid = []
        for r in raw:
            observed["total_variants"] += 1
            if not child_nonref(r, child):
                continue
            observed["child_alt_variants"] += 1
            gts = [gt(r, p) for p in (child, father, mother)]
            if not all(diploid_called(g) for g in gts) or "MCU" in r["info"]:
                continue
            observed["testable_child_variants"] += 1
            if "MCV" in r["info"]:
                observed["observed_inconsistent"] += 1
                valid.append(r)
            else:
                observed["observed_consistent"] += 1

        truth_available = (root / f"{cls}.truth.vcf.gz").exists() if not sv else (root / f"{cls}.truvari").is_dir()
        truth_available = truth_available and bool(beds["sv" if sv else "small"])
        if truth_available:
            comparison = (truvari_verdict(root / f"{cls}.truvari") if sv else
                          happy_verdict(root / f"{cls}.truth.vcf.gz"))
        else:
            comparison = {}
        adjudication = Counter()
        for r in valid:
            if not truth_available:
                status = "truth_unavailable"
            elif not confidence(r, beds["sv" if sv else "small"], sv):
                status = "outside_confidence"
            else:
                status = status_for(r, comparison)
            adjudication[status] += 1
            per_variant.append({
                "reference": reference, "trio": trio, "class": cls,
                "chrom": r["chrom"], "pos": r["pos"], "id": r["id"],
                "child_gt": "/".join(gt(r, child)), "status": status,
            })

        testable = observed["testable_child_variants"]
        inconsistent = observed["observed_inconsistent"]
        true_n = adjudication["true"] if truth_available else None
        false_n = adjudication["false"] if truth_available else None
        evaluated = true_n + false_n if truth_available else None
        results.append({
            "reference": reference, "trio": trio, "class": cls,
            "total_variants": observed["total_variants"],
            "child_alt_variants": observed["child_alt_variants"],
            "testable_child_variants": testable,
            "observed_consistent": observed["observed_consistent"],
            "observed_inconsistent": inconsistent,
            "observed_inconsistent_rate_pct": (100 * inconsistent / testable if testable else None),
            "intersected_with_truth_set": evaluated,
            "true_inconsistent": true_n, "false_inconsistent": false_n,
            "false_inconsistent_rate_pct": (100 * false_n / evaluated if evaluated else None),
            "outside_confidence": adjudication["outside_confidence"] if truth_available else None,
            "unresolved_inside_confidence": adjudication["unknown"] if truth_available else None,
            "truth_unavailable": adjudication["truth_unavailable"],
        })
    return results, per_variant


def write_tsv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8") as dest:
        writer = csv.DictWriter(dest, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("NA" if v is None else f"{v:.4f}" if isinstance(v, float) else v)
                             for k, v in row.items()})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--trio", required=True)
    p.add_argument("--child", required=True)
    p.add_argument("--father", required=True)
    p.add_argument("--mother", required=True)
    p.add_argument("--small-bed", required=True)
    p.add_argument("--sv-bed", default="")
    args = p.parse_args()
    rows, details = calculate(
        args.directory, args.reference, args.trio, args.child, args.father,
        args.mother, args.small_bed, args.sv_bed,
    )
    write_tsv(Path(args.directory) / "mendelian_summary.tsv", rows)
    write_tsv(Path(args.directory) / "mendelian_variant_adjudication.tsv", details)


if __name__ == "__main__":
    main()
