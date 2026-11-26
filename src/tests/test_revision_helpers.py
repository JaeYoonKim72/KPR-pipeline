#!/usr/bin/env python3
"""Regression tests for misjoin classification and AMY structure labels."""
import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

SRC = pathlib.Path(__file__).resolve().parents[1]

class RevisionHelpers(unittest.TestCase):
    def test_misjoins(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = pathlib.Path(tmp)
            bed = d / "arms.bed"
            bed.write_text("".join(f"chr{c}\t0\t3000000\n" for c in (13, 14, 15, 21, 22)))
            def alignment(contig, chrom, start):
                return f"J\t{contig}\t5000000\t0\t1000000\t+\t{chrom}\t10000000\t{start}\t{start+1000000}\t1000000\t1000000\t60\n"
            report = d / "report.txt"
            report.write_text(alignment("keep", "chr13", 0) + alignment("keep", "chr14", 0)
                              + alignment("remove", "chr13", 0) + alignment("remove", "chrX", 0)
                              + alignment("remove2", "chr13", 4000000) + alignment("remove2", "chr14", 0))
            excluded = d / "excluded.txt"
            result = subprocess.run([sys.executable, str(SRC / "filter_misjoins.py"),
                "--input", str(report), "--assembly", "sample", "--bed", str(bed),
                "--exclude", str(excluded)], capture_output=True, text=True, check=True)
            self.assertEqual(excluded.read_text().splitlines(), ["remove", "remove2"])
            self.assertEqual([l.split("\t")[-1] for l in result.stdout.splitlines()],
                             ["ACROCENTRIC", "EXCLUDE", "EXCLUDE"])
    def test_structure_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = pathlib.Path(tmp)
            structures = d / "structures.tsv"
            structures.write_text("haplotype\tn_genes\tstructure\na\t1\tAMY1A+\nb\t1\tAMY1A+\nc\t1\tAMY1B-\n")
            mapping = d / "map.tsv"
            mapping.write_text("structure\thaplotype_class\nAMY1A+\tH1A1-S1\nAMY1B-\tH1A1-S2\n")
            out = d / "classes.tsv"
            assignments = d / "assignments.tsv"
            subprocess.run([sys.executable, str(SRC / "classify_amy_structures.py"),
                "--structures", str(structures), "--mapping", str(mapping),
                "--output", str(out), "--assignments", str(assignments)], check=True)
            with out.open() as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual([(r["haplotype_class"], r["n_haplotypes"]) for r in rows],
                             [("H1A1-S1", "2"), ("H1A1-S2", "1")])

if __name__ == "__main__":
    unittest.main()
