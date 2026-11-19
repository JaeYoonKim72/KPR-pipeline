"""Regression tests for KPR-private variant count proportions."""

import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest


SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "summarize_variant_comparison.py"


class VariantComparisonSummaryTests(unittest.TestCase):
    def test_private_proportion_uses_all_three_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            counts = root / "counts"
            out = root / "results"
            counts.mkdir()
            out.mkdir()
            with (counts / "chr1.tsv").open("w") as handle:
                for category, value in (
                    ("target_only", 7),
                    ("shared", 45),
                    ("other_only", 50),
                ):
                    handle.write(f"chr1\tsv\tlocus\tINS\t{category}\t{value}\n")

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--counts", str(counts),
                    "--out", str(out),
                    "--target-label", "KPR",
                    "--other-label", "HPRC",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with (out / "variant_comparison_summary.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["KPR_private"], "7")
            self.assertEqual(row["shared"], "45")
            self.assertEqual(row["HPRC_private"], "50")
            self.assertEqual(row["proportion_KPR_private_pct"], "6.86")
            self.assertAlmostEqual(
                float(row["proportion_KPR_private_pct"]),
                100 * 7 / (7 + 45 + 50),
                delta=0.005,
            )

    def test_no_private_variants_reports_zero_percent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            counts = root / "counts"
            out = root / "results"
            counts.mkdir()
            out.mkdir()
            (counts / "chr2.tsv").write_text(
                "chr2\tsmall\tallele_exact\tSNV\tshared\t4\n"
            )
            subprocess.run(
                [sys.executable, str(SCRIPT), "--counts", str(counts), "--out", str(out)],
                check=True,
                capture_output=True,
                text=True,
            )
            with (out / "variant_comparison_summary.tsv").open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["proportion_KPR_private_pct"], "0.00")


if __name__ == "__main__":
    unittest.main()
