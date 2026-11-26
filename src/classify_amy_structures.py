#!/usr/bin/env python3
"""Group identical gene-tag structures and apply a supplied haplotype nomenclature."""
import argparse
import csv
from collections import defaultdict

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures", required=True)
    p.add_argument("--mapping", default="")
    p.add_argument("--output", required=True)
    p.add_argument("--assignments", required=True)
    a = p.parse_args()
    groups = defaultdict(list)
    with open(a.structures) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            groups[row["structure"]].append(row["haplotype"])
    names = {}
    if a.mapping:
        with open(a.mapping) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["structure"] in names:
                    raise ValueError("Duplicate structure in haplotype mapping")
                names[row["structure"]] = row["haplotype_class"]
        if set(groups) - set(names):
            raise ValueError("Haplotype mapping does not cover all structures")
    else:
        names = {s: f"Structure{i:03d}" for i, s in enumerate(sorted(groups), 1)}
    if len(set(names[s] for s in groups)) != len(groups):
        raise ValueError("Haplotype classes must identify distinct structures")
    with open(a.output, "w") as f, open(a.assignments, "w") as g:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        assign = csv.writer(g, delimiter="\t", lineterminator="\n")
        out.writerow(["haplotype_class", "n_haplotypes", "example", "structure"])
        assign.writerow(["haplotype", "haplotype_class", "structure"])
        for s in sorted(groups, key=lambda x: (-len(groups[x]), names[x])):
            ids = sorted(groups[s])
            out.writerow([names[s], len(ids), ids[0], s])
            for hap in ids:
                assign.writerow([hap, names[s], s])

if __name__ == "__main__":
    main()
