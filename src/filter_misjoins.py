#!/usr/bin/env python3
"""Classify interchromosomal pairs reported by paftools.js misjoin."""
import argparse
from collections import defaultdict
from pathlib import Path

def run(args):
    arms = defaultdict(list)
    if args.bed:
        with open(args.bed) as stream:
            for line in stream:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end = line.split()[:3]
                start, end = int(start), int(end)
                if start < 0 or end <= start:
                    raise ValueError("Invalid BED interval")
                arms[chrom].append((start, end))
        required = {"chr13", "chr14", "chr15", "chr21", "chr22"}
        if not required.issubset(arms):
            raise ValueError("BED must contain all five CHM13 acrocentric short arms")
    def overlap(row):
        return row[5] in {"chr13", "chr14", "chr15", "chr21", "chr22"} and any(
            start < int(row[8]) and end > int(row[7]) for start, end in arms.get(row[5], []))
    excluded = set()
    pairs = []
    with open(args.input) as stream:
        for line in stream:
            row = line.split()
            if row and row[0] in ("J", "j"):
                if len(row) < 13:
                    raise ValueError("Invalid paftools misjoin record")
                pairs.append(row[1:13])
    if len(pairs) % 2:
        raise ValueError("Unpaired paftools misjoin record")
    for left, right in zip(pairs[::2], pairs[1::2]):
        if left[0] != right[0] or left[5] == right[5]:
            raise ValueError("Invalid interchromosomal pair")
        status = "UNCLASSIFIED" if not args.bed else (
            "ACROCENTRIC" if overlap(left) and overlap(right) else "EXCLUDE")
        if status == "EXCLUDE":
            excluded.add(left[0])
        print(args.assembly, left[0], left[1], left[5], *left[7:9],
              right[5], *right[7:9], status, sep="\t")
    Path(args.exclude).write_text("".join(x + "\n" for x in sorted(excluded)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--bed", default="")
    parser.add_argument("--exclude", required=True)
    run(parser.parse_args())
