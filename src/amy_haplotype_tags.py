#!/usr/bin/env python3
"""Extract ordered gene tags and orientations from GFA walks."""

import re
import sys
from collections import OrderedDict

STEP = re.compile(r"([<>])([^<>]+)")


def read_node_genes(path):
    """Return a mapping of node identifier to gene name."""
    node_gene = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[1]:
                continue
            node_gene[fields[0]] = fields[1]
    return node_gene


def read_node_lengths(gfa_path):
    """Return a mapping of node identifier to sequence length."""
    lengths = {}
    with open(gfa_path) as handle:
        for line in handle:
            if not line.startswith("S"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                lengths[fields[1]] = len(fields[2])
    return lengths


def read_walks(gfa_path):
    """Return the ordered walk of every haplotype as (name, [(orient, node)])."""
    walks = OrderedDict()
    with open(gfa_path) as handle:
        for line in handle:
            if not line.startswith("W"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue

            sample, haplotype = fields[1], fields[2]
            name = "{}#{}".format(sample, haplotype)
            steps = [(m.group(1), m.group(2)) for m in STEP.finditer(fields[6])]
            walks.setdefault(name, []).extend(steps)
    return walks


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python3 amy_haplotype_tags.py "
                 "<subgraph.gfa> <node_gene.tsv> <output_prefix>")

    gfa_path, node_gene_path, prefix = sys.argv[1:4]

    node_gene = read_node_genes(node_gene_path)
    lengths = read_node_lengths(gfa_path)
    walks = read_walks(gfa_path)

    long_path = prefix + ".gene_tags.tsv"
    wide_path = prefix + ".haplotype_structure.tsv"

    with open(long_path, "w") as long_out, open(wide_path, "w") as wide_out:
        long_out.write("haplotype\torder\tnode\torientation\tgene\tlength\n")
        wide_out.write("haplotype\tn_genes\tstructure\n")

        for name, steps in walks.items():
            tags = []
            order = 0
            previous = None

            for orient, node in steps:
                gene = node_gene.get(node)
                if gene is None:
                    previous = None
                    continue

                strand = "+" if orient == ">" else "-"
                tag = gene + strand

                # Collapse consecutive identical gene tags.
                if tag == previous:
                    continue

                order += 1
                previous = tag
                tags.append(tag)
                long_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    name, order, node, strand, gene, lengths.get(node, 0)))

            wide_out.write("{}\t{}\t{}\n".format(name, len(tags), ",".join(tags)))

    print("wrote {}".format(long_path))
    print("wrote {}".format(wide_path))


if __name__ == "__main__":
    main()
