#!/usr/bin/env python3
"""Select SVIM-asm VCF records supported by independent Truvari bench matches.

Each --matched-vcf is the tp-base.vcf.gz produced when SVIM-asm was the
Truvari baseline and one read-based caller was the comparison.
"""

import argparse
from collections import Counter
import gzip
from pathlib import Path
import sys


def lines(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8")


def key(fields):
    if len(fields) < 8:
        raise ValueError("VCF record has fewer than eight fields")
    return (fields[0], fields[1], fields[3], fields[4])


def load_matches(path):
    records = set()
    with lines(path) as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            rec = key(line.rstrip("\r\n").split("\t"))
            if rec in records:
                raise ValueError(f"Duplicated SVIM-asm base allele in Truvari result: {rec}")
            records.add(rec)
    return records


def select(base, matches, output, audit, min_support, strong=None):
    if min_support not in (1, 2, 3):
        raise ValueError("Minimum read-based callers must be 1, 2 or 3")
    if set(matches) != {"pbsv", "sniffles", "svim"}:
        raise ValueError("Exactly one matched VCF for each of pbsv, sniffles and svim is required")
    compared = {caller: load_matches(path) for caller, path in matches.items()}
    strong_set = None
    if strong is not None:
        if set(strong) != {"pbsv", "sniffles", "svim"}:
            raise ValueError("Require strong-match VCFs from all three read-based callers")
        strong_set = set().union(*(load_matches(path) for path in strong.values()))
    stats = Counter()
    seen = set()
    with lines(base) as inp, open(output, "wt", encoding="utf-8") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            variant = key(fields)
            if variant in seen:
                raise ValueError(f"Duplicate SVIM-asm baseline variant: {variant}")
            seen.add(variant)
            stats["SVIM_asm_baseline"] += 1
            supported = [caller for caller in ("pbsv", "sniffles", "svim") if variant in compared[caller]]
            stats[f"matched_read_callers_{len(supported)}"] += 1
            if len(supported) < min_support:
                stats["excluded_insufficient_read_callers"] += 1
                continue
            if strong_set is not None and variant not in strong_set:
                stats["excluded_no_strong_read_support"] += 1
                continue
            stats["retained"] += 1
            out.write(line)
    unexpected = set().union(*compared.values()) - seen
    if unexpected:
        raise ValueError(
            f"Truvari tp-base contains {len(unexpected)} variants absent from the exact SVIM-asm baseline; "
            "check input VCF version/normalization"
        )
    if strong_set is not None and strong_set - seen:
        raise ValueError("Strong-match calls do not exist in baseline SVIM-asm VCF")
    with open(audit, "wt", encoding="utf-8") as stream:
        stream.write("metric\tcount\n")
        for name in (
            "SVIM_asm_baseline",
            "matched_read_callers_0",
            "matched_read_callers_1",
            "matched_read_callers_2",
            "matched_read_callers_3",
            "excluded_insufficient_read_callers",
            "excluded_no_strong_read_support",
            "retained",
        ):
            stream.write(f"{name}\t{stats[name]}\n")
        for caller in ("pbsv", "sniffles", "svim"):
            stream.write(f"truvari_matched_{caller}\t{len(compared[caller])}\n")
    return stats["retained"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True, help="filtered, sorted SVIM-asm baseline VCF(.gz)")
    p.add_argument("--matched-vcf", action="append", required=True, metavar="CALLER=VCF")
    p.add_argument("--strong-matched-vcf", action="append", default=[], metavar="CALLER=VCF")
    p.add_argument("--min-read-callers", type=int, default=2, choices=(1, 2, 3))
    p.add_argument("--output", required=True)
    p.add_argument("--audit", required=True)
    a = p.parse_args()

    def parse_args(items):
        out = {}
        for item in items:
            if "=" not in item:
                raise ValueError("Each --matched-vcf requires CALLER=VCF format")
            caller, path = item.split("=", 1)
            if caller in out:
                raise ValueError(f"Duplicated matched caller: {caller}")
            out[caller] = path
        return out

    try:
        selected = select(
            a.base, parse_args(a.matched_vcf), a.output, a.audit,
            a.min_read_callers,
            parse_args(a.strong_matched_vcf) if a.strong_matched_vcf else None,
        )
    except (OSError, ValueError) as e:
        p.exit(1, f"ERROR: {e}\n")
    print(f"Selected {selected} SVIM-asm calls with >= {a.min_read_callers} independent Truvari read-caller matches", file=sys.stderr)


if __name__ == "__main__":
    main()
