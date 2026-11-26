#!/usr/bin/env python3
"""Collect unique node IDs traversed by each sample in GFA W lines."""

import re
import sys

SEPARATOR = re.compile(r"[<>]")


def collect_nodes(gfa_path):
    """Return a mapping of sample name to the set of nodes it traverses."""
    per_sample = {}

    with open(gfa_path) as handle:
        for line in handle:
            if not line.startswith("W"):
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                continue

            sample = fields[1]
            nodes = (node for node in SEPARATOR.split(fields[6]) if node)
            per_sample.setdefault(sample, set()).update(nodes)

    return per_sample


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python3 wline_sample_nodes.py <input.gfa> <output_prefix>")

    gfa_path, prefix = sys.argv[1], sys.argv[2]
    per_sample = collect_nodes(gfa_path)

    for sample, nodes in sorted(per_sample.items()):
        out_path = "{}.{}.nodes.txt".format(prefix, sample)
        with open(out_path, "w") as handle:
            for node in sorted(nodes, key=lambda n: (len(n), n)):
                handle.write(node + "\n")
        print("{}\t{}\t{}".format(sample, len(nodes), out_path))


if __name__ == "__main__":
    main()
