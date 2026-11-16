#!/usr/bin/env python3
"""Convert a multiple sequence alignment of gene copies into a VCF of
substitutions, in the coordinates of the reference sequence included in the
alignment.

Alignment columns where the reference carries a gap are skipped, so the
positions reported are those of the reference gene. Only substitutions are
emitted; insertions and deletions are ignored. Every gene copy is one
observation and is written as a haploid genotype, so that a haplotype carrying
several copies of a gene contributes one observation per copy.

Usage:
    python3 msa_to_vcf.py <alignment.fasta> <reference_name> <chromosome> <output.vcf>

    reference_name  identifier of the reference sequence inside the alignment
    chromosome      name written to the CHROM column, e.g. the gene name
"""

import sys

BASES = {"A", "C", "G", "T"}


def read_fasta(path):
    """Return the alignment as a list of (name, sequence)."""
    records = []
    name, chunks = None, []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line.strip())

    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: python3 msa_to_vcf.py <alignment.fasta> "
                 "<reference_name> <chromosome> <output.vcf>")

    aln_path, ref_name, chrom, out_path = sys.argv[1:5]
    records = read_fasta(aln_path)
    if not records:
        sys.exit("empty alignment: {}".format(aln_path))

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        sys.exit("sequences of unequal length; the input must be aligned")

    reference = None
    samples = []
    for name, seq in records:
        if name == ref_name and reference is None:
            reference = seq.upper()
        else:
            samples.append((name, seq.upper()))

    if reference is None:
        sys.exit("reference {} not found in the alignment".format(ref_name))
    if not samples:
        sys.exit("the alignment holds no sequence other than the reference")

    sample_names = [name for name, _ in samples]
    n_sites = 0
    position = 0

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("##contig=<ID={},length={}>\n".format(
            chrom, sum(1 for base in reference if base != "-")))
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(sample_names) + "\n")

        for column in range(len(reference)):
            ref_base = reference[column]
            if ref_base == "-":
                continue

            position += 1
            if ref_base not in BASES:
                continue

            observed = [seq[column] for _, seq in samples]
            alt_alleles = []
            for base in observed:
                if base in BASES and base != ref_base and base not in alt_alleles:
                    alt_alleles.append(base)

            if not alt_alleles:
                continue

            genotypes = []
            for base in observed:
                if base == ref_base:
                    genotypes.append("0")
                elif base in alt_alleles:
                    genotypes.append(str(alt_alleles.index(base) + 1))
                else:
                    genotypes.append(".")

            out.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n".format(
                chrom, position, ref_base, ",".join(alt_alleles),
                "\t".join(genotypes)))
            n_sites += 1

    print("{}\t{} copies\t{} sites\t{} bp".format(
        chrom, len(samples), n_sites, position))


if __name__ == "__main__":
    main()
