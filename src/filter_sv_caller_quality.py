#!/usr/bin/env python3
"""Filter caller-specific SV VCFs before Truvari matching.

For single-sample VCFs, exclude records missing required depth or read support and report them in the audit.
"""

import argparse
from collections import Counter
import gzip
import math
from pathlib import Path
import sys


def read_lines(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def numeric(value):
    if value is None or str(value) in {"", "."}:
        return None
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) else None


def info_fields(raw):
    out = {}
    if raw != ".":
        for item in raw.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                out[k] = v
            elif item:
                out[item] = True
    return out


def sample_fields(fields):
    if len(fields) < 10 or fields[8] in {"", "."}:
        return {}
    labels = fields[8].split(":")
    values = fields[9].split(":")
    return dict(zip(labels, values))


def allele_depth(genotype, info):
    for k in ("AD",):
        value = genotype.get(k, info.get(k))
        if value and value != ".":
            parts = [numeric(i) for i in str(value).split(",")]
            if len(parts) > 1 and all(i is not None for i in parts):
                return sum(parts), sum(parts[1:])
    for k in ("DV", "AO"):
        val = numeric(genotype.get(k, info.get(k)))
        if val is not None:
            ref = numeric(genotype.get("DR", info.get("DR")))
            return (val + ref if ref is not None else None), val
    return None, None


def depth_and_alt(caller, genotype, info):
    _, ad_alt = allele_depth(genotype, info)
    depth = numeric(genotype.get("DP", info.get("DP")))
    if caller == "svim":
        alt = numeric(info.get("SUPPORT"))
        if alt is None:
            alt = numeric(genotype.get("SUPPORT"))
        if alt is None:
            alt = ad_alt
    elif caller == "sniffles":
        alt = numeric(genotype.get("DV"))
        dr = numeric(genotype.get("DR"))
        depth = (alt + dr) if alt is not None and dr is not None else None
    else:
        alt = ad_alt
        if alt is None:
            alt = numeric(genotype.get("DV", info.get("SUPPORT")))
    return depth, alt


def sv_length(fields, info):
    val = info.get("SVLEN")
    if val:
        for part in str(val).split(","):
            num = numeric(part)
            if num is not None and num != 0:
                return abs(num)
    ref, alt = fields[3], fields[4]
    if not (alt.startswith("<") and alt.endswith(">")):
        return abs(len(alt) - len(ref))
    end = numeric(info.get("END"))
    if end is not None and info.get("SVTYPE") in {"DEL", "DUP", "INV"}:
        return abs(end - int(fields[1]))
    return None


def filter_record(fields, caller, min_dp, min_alt, sniffles_gq, svim_qual, min_sv_size=51):
    chrom = fields[0].removeprefix("chr")
    if not chrom.isdigit() or not 1 <= int(chrom) <= 22:
        return "non_autosomal"
    if fields[6] != "PASS":
        return "filter_not_pass"
    info = info_fields(fields[7])
    length = sv_length(fields, info)
    if length is None:
        return "missing_sv_length"
    if length < min_sv_size:
        return f"sv_length_below_{min_sv_size}"
    if caller == "svimasm":
        return None
    fmt = sample_fields(fields)
    if caller == "sniffles":
        gq = numeric(fmt.get("GQ"))
        if gq is None:
            return "missing_genotype_quality"
        if gq < sniffles_gq:
            return "genotype_quality_below_threshold"
    if caller == "svim":
        qual = numeric(fields[5])
        if qual is None:
            return "missing_svim_qual"
        if qual < svim_qual:
            return "svim_qual_below_threshold"
    depth, alt = depth_and_alt(caller, fmt, info)
    if depth is None:
        return "missing_depth"
    if alt is None:
        return "missing_alt_support"
    if depth < min_dp:
        return "depth_below_threshold"
    if alt < min_alt:
        return "alt_support_below_threshold"
    return None


def process(input_path, output_path, audit_path, caller, min_dp, min_alt,
            sniffles_gq, svim_qual, min_sv_size=51):
    stats = Counter()
    with read_lines(input_path) as source, open(output_path, "wt", encoding="utf-8") as target:
        n_samples = None
        for line in source:
            if line.startswith("#CHROM"):
                n_samples = max(0, len(line.rstrip("\n").split("\t")) - 9)
                if caller != "svimasm" and n_samples != 1:
                    raise ValueError(f"{caller}: requires one sample; got {n_samples}")
                target.write(line)
                continue
            if line.startswith("#"):
                target.write(line)
                continue
            if n_samples is None:
                raise ValueError("VCF lacks #CHROM header")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Malformed VCF record with fewer than 8 fields")
            if "," in fields[4]:
                raise ValueError("Multiallelic VCF: split first with bcftools norm -m -any")
            stats["input"] += 1
            failure = filter_record(fields, caller, min_dp, min_alt, sniffles_gq, svim_qual, min_sv_size)
            if failure:
                stats[failure] += 1
            else:
                stats["retained"] += 1
                target.write(line)
        if n_samples is None:
            raise ValueError("VCF lacks #CHROM header")
    with open(audit_path, "wt", encoding="utf-8") as out:
        out.write("reason\tcount\n")
        for name in ["input", "retained"] + sorted(set(stats) - {"input", "retained"}):
            out.write(f"{name}\t{stats[name]}\n")
    missing = sum(v for k, v in stats.items() if k.startswith("missing_"))
    if missing:
        print(
            f"WARNING: {caller} excluded {missing} SV records for missing required fields; "
            f"inspect {audit_path} and source FORMAT/INFO definitions.",
            file=sys.stderr,
        )
    return stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--caller", required=True, choices=("pbsv", "sniffles", "svim", "svimasm"))
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--audit", required=True)
    p.add_argument("--min-size", type=int, default=51)
    p.add_argument("--min-dp", type=int, default=7)
    p.add_argument("--min-alt", type=int, default=3)
    p.add_argument("--sniffles-gq", type=int, default=20)
    p.add_argument("--svim-qual", type=int, default=10)
    a = p.parse_args()
    try:
        stats = process(
            a.input, a.output, a.audit, a.caller, a.min_dp, a.min_alt,
            a.sniffles_gq, a.svim_qual, a.min_size,
        )
    except (OSError, ValueError) as e:
        p.exit(1, f"ERROR: {e}\n")
    print(f"{a.caller}: retained {stats['retained']}/{stats['input']} SV records", file=sys.stderr)


if __name__ == "__main__":
    main()
