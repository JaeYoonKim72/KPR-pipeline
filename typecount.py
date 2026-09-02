#!/usr/bin/env python3
"""Classify every variant of a VCF by the number of alternative alleles it
carries across all samples.

Each record is assigned to exactly one of the following classes:

    singleton   exactly one alternative allele across all samples
    doubleton   exactly two alternative alleles across all samples
    common      more than two alternative alleles, frequency below 95%
    core        more than two alternative alleles, frequency of 95% or above

Polymorphic variants are all non-singleton variants, that is the sum of the
doubleton, common and core classes.

Usage:
    python3 typecount.py <input.vcf.gz>

The result is written next to the input file with the extension .typecount.
"""

import gzip
import os
import sys

CORE_FREQUENCY = 0.95


def parse_vcf(vcf_path):
    """Count the variants of each class in a bgzipped VCF."""
    counts = {
        "singleton": 0,
        "doubleton": 0,
        "common": 0,
        "core": 0,
        "polymorphism": 0,
    }
    total = 0
    samples = []

    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break

        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            format_keys = fields[8].split(":")
            if "GT" not in format_keys:
                continue
            gt_index = format_keys.index("GT")

            alleles = []
            for sample_field in fields[9:]:
                sample_values = sample_field.split(":")
                if len(sample_values) <= gt_index:
                    continue

                genotype = sample_values[gt_index]
                if "/" in genotype:
                    called = genotype.split("/")
                elif "|" in genotype:
                    called = genotype.split("|")
                else:
                    continue

                for allele in called:
                    if allele.isdigit():
                        alleles.append(int(allele))

            if not alleles:
                continue

            alt_count = sum(1 for allele in alleles if allele > 0)

            if alt_count == 1:
                variant_class = "singleton"
            elif alt_count == 2:
                variant_class = "doubleton"
            elif alt_count > 2:
                frequency = alt_count / len(alleles)
                variant_class = "core" if frequency >= CORE_FREQUENCY else "common"
            else:
                continue

            counts[variant_class] += 1
            if variant_class != "singleton":
                counts["polymorphism"] += 1
            total += 1

    return counts, total, samples


def write_results(counts, vcf_path):
    """Write the counts as a single-row tab-separated table."""
    base = os.path.basename(vcf_path)
    if base.endswith(".vcf.gz"):
        out_name = base[: -len(".vcf.gz")] + ".typecount"
    elif base.endswith(".gz"):
        out_name = base[: -len(".gz")] + ".typecount"
    else:
        out_name = base + ".typecount"

    out_path = os.path.join(os.path.dirname(vcf_path) or ".", out_name)

    header = ["file", "singleton", "doubleton", "common", "core", "polymorphism"]
    row = [
        base,
        counts["singleton"],
        counts["doubleton"],
        counts["common"],
        counts["core"],
        counts["polymorphism"],
    ]

    with open(out_path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        handle.write("\t".join(str(value) for value in row) + "\n")

    return out_path


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 typecount.py <input.vcf.gz>")

    vcf_path = sys.argv[1]
    if not os.path.isfile(vcf_path):
        sys.exit("File not found: {}".format(vcf_path))

    counts, total, samples = parse_vcf(vcf_path)
    out_path = write_results(counts, vcf_path)

    print("File: {}".format(vcf_path))
    print("Samples: {}".format(len(samples)))
    print("Variants: {}".format(total))
    for name in ("singleton", "doubleton", "common", "core", "polymorphism"):
        print("  {:<12} {}".format(name, counts[name]))
    print("Written to: {}".format(out_path))


if __name__ == "__main__":
    main()
