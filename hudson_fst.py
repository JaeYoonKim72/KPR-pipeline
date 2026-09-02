#!/usr/bin/env python3
"""Compute Hudson's FST between two groups of samples from a biallelic VCF.

Hudson's estimator is a ratio of averages: the numerators and the denominators
of the per-site estimates are summed separately and divided at the end. It is
robust to unequal sample sizes between the two groups, which matters when the
groups differ several fold in size.

For every biallelic site with allele frequencies p1 and p2 and sample sizes n1
and n2,

    numerator   = (p1 - p2)^2 - p1(1 - p1)/(n1 - 1) - p2(1 - p2)/(n2 - 1)
    denominator = p1(1 - p2) + p2(1 - p1)

and FST is the sum of the numerators divided by the sum of the denominators.

Usage:
    python3 hudson_fst.py <input.vcf[.gz]> <group1.list> <group2.list> [label1 label2]

Each group file holds one sample name per line. Haploid and diploid genotypes
are both accepted; every called allele counts as one observation.
"""

import gzip
import sys


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_list(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def allele_counts(genotype):
    """Return (alternative allele count, total called alleles)."""
    field = genotype.split(":")[0]
    if field in (".", "./.", ".|."):
        return 0, 0

    alleles = field.replace("|", "/").split("/")
    alt, total = 0, 0
    for allele in alleles:
        if allele == ".":
            continue
        total += 1
        if allele != "0":
            alt += 1
    return alt, total


def main():
    if len(sys.argv) not in (4, 6):
        sys.exit("Usage: python3 hudson_fst.py <input.vcf[.gz]> "
                 "<group1.list> <group2.list> [label1 label2]")

    vcf_path, list1, list2 = sys.argv[1:4]
    label1, label2 = (sys.argv[4], sys.argv[5]) if len(sys.argv) == 6 else ("pop1", "pop2")

    group1, group2 = set(read_list(list1)), set(read_list(list2))
    index1, index2 = [], []

    num_sum, den_sum, used = 0.0, 0.0, 0

    with open_vcf(vcf_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            fields = line.rstrip("\n").split("\t")

            if line.startswith("#CHROM"):
                samples = fields[9:]
                index1 = [i for i, s in enumerate(samples) if s in group1]
                index2 = [i for i, s in enumerate(samples) if s in group2]
                if not index1 or not index2:
                    sys.exit("one of the groups has no sample in the VCF")
                continue

            if "," in fields[4]:          # multiallelic, skipped
                continue

            calls = fields[9:]

            alt1 = tot1 = alt2 = tot2 = 0
            for i in index1:
                a, t = allele_counts(calls[i])
                alt1 += a
                tot1 += t
            for i in index2:
                a, t = allele_counts(calls[i])
                alt2 += a
                tot2 += t

            if tot1 < 2 or tot2 < 2:
                continue

            p1, p2 = alt1 / tot1, alt2 / tot2

            numerator = ((p1 - p2) ** 2
                         - p1 * (1 - p1) / (tot1 - 1)
                         - p2 * (1 - p2) / (tot2 - 1))
            denominator = p1 * (1 - p2) + p2 * (1 - p1)

            if denominator == 0:
                continue

            num_sum += numerator
            den_sum += denominator
            used += 1

    fst = num_sum / den_sum if den_sum else float("nan")
    print("pop1\tpop2\tsites\tHudson_FST")
    print("{}\t{}\t{}\t{:.6f}".format(label1, label2, used, fst))


if __name__ == "__main__":
    main()
