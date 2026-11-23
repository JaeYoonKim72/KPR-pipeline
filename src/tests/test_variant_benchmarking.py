"""Regression tests for variant benchmarking and trio adjudication."""

import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from summarize_benchmark import truvari_row, hap_row
from summarize_mendelian import calculate, confidence, read_bed
from classify_benchmark_fn import classify_sv, classify_small

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000000>\n"
    "##INFO=<ID=MCV,Number=1,Type=Integer,Description=\"Violation\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG003\tHG004\n"
)


def write_gzip(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        t = tempfile.TemporaryDirectory()
        self.addCleanup(t.cleanup)
        self.d = Path(t.name)

    def test_sv_gt_aware_metrics(self):
        p = self.d / "summary.json"
        p.write_text(json.dumps({
            "base cnt": 10, "comp cnt": 12, "TP-base": 7, "TP-comp": 7,
            "TP-base_TP-gt": 5, "TP-comp_TP-gt": 5,
            "TP-base_FP-gt": 2, "TP-comp_FP-gt": 2,
            "FP": 5, "FN": 3, "gt_concordance": 5/7,
        }))
        out = truvari_row(p)
        self.assertEqual((out["TP"], out["FN"], out["FP"]), (5, 5, 7))
        self.assertAlmostEqual(out["Precision"], 5/12)
        self.assertAlmostEqual(out["Recall"], 5/10)
        self.assertAlmostEqual(out["Presence"], 7/10)
        self.assertAlmostEqual(out["GT_concordance"], 5/7)

    def test_small_presence_and_failure_cause(self):
        p = self.d / "a.SNP.summary.csv"
        p.write_text(
            "Type,Filter,TRUTH.TOTAL,QUERY.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.FP,"
            "METRIC.Precision,METRIC.Recall,METRIC.F1_Score\n"
            "SNP,ALL,3,3,2,1,1,0.666667,0.666667,0.666667\n"
        )
        annotated = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n"
            "chr1\t10\ta\tA\tG\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm\n"
            "chr1\t20\tb\tA\tG\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm\n"
            "chr1\t30\tc\tA\tG\t.\tPASS\t.\tGT:BD:BK\t0/1:FN:am\t1/1:FP:am\n"
        )
        write_gzip(self.d / "a.SNP.vcf.gz", annotated)
        out = hap_row(p, "SNP")
        self.assertEqual(out["FN_matching_genotype"], 1)
        self.assertAlmostEqual(out["Presence"], 1)
        self.assertAlmostEqual(out["GT_concordance"], 2/3)

    def test_mendelian_confidence_and_true_false(self):
        bed = self.d / "truth.bed"
        bed.write_text("chr1\t0\t1000\n")
        annotated = self.d / "SNP.annot.vcf.gz"
        write_gzip(annotated, HEADER +
                   "chr1\t10\tv1\tA\tG\t.\tPASS\tMCV=1\tGT\t0/1\t0/0\t0/0\n" +
                   "chr1\t20\tv2\tA\tG\t.\tPASS\tMCV=1\tGT\t1/1\t0/0\t0/0\n" +
                   "chr1\t30\tv3\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0\n" +
                   "chr1\t40\tv4\tA\tG\t.\tPASS\tMCV=1\tGT\t0/1\t./.\t0/0\n")
        hap = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n"
            "chr1\t10\tv1\tA\tG\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm\n"
            "chr1\t20\tv2\tA\tG\t.\tPASS\t.\tGT:BD:BK\t0/1:FN:am\t1/1:FP:am\n"
        )
        write_gzip(self.d / "SNP.truth.vcf.gz", hap)
        for cls in ("Indel", "INS", "DEL"):
            write_gzip(self.d / f"{cls}.annot.vcf.gz", HEADER)
        rows, details = calculate(self.d, "graph", "AJ", "HG002", "HG003", "HG004", bed, "")
        r = rows[0]
        self.assertEqual(r["observed_inconsistent"], 2)
        self.assertEqual(r["testable_child_variants"], 3)
        self.assertEqual((r["true_inconsistent"], r["false_inconsistent"]), (1, 1))
        self.assertEqual(r["false_inconsistent_rate_pct"], 50)
        self.assertIsNone(rows[2]["false_inconsistent_rate_pct"])

    def test_sv_confidence_requires_complete_interval(self):
        bed = self.d / "sv.bed"
        bed.write_text("chr1\t100\t200\n")
        region = read_bed(str(bed))
        v = {"chrom": "chr1", "pos": 150, "ref": "A", "info": {"END": "210"}}
        self.assertFalse(confidence(v, region, True))
        v["info"]["END"] = "190"
        self.assertTrue(confidence(v, region, True))

    def test_fn_sv_and_gbz_missing(self):
        fn = self.d / "fn.vcf.gz"
        tpb = self.d / "tp-base.vcf.gz"
        query = self.d / "q.vcf.gz"
        base = ("##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        write_gzip(fn, base +
                   "chr1\t100\t.\tA\t<INS>\t.\tPASS\tSVLEN=100\n" +
                   "chr1\t1000\t.\tA\t<INS>\t.\tPASS\tSVLEN=100\n")
        write_gzip(tpb, base)
        write_gzip(query, base + "chr1\t120\t.\tA\t<INS>\t.\tPASS\tSVLEN=100\n")
        calls = list(classify_sv(fn, tpb, query, "kpr"))
        self.assertEqual([x[1] for x in calls], ["unresolved", "detection_failure"])
        linear = list(classify_sv(fn, tpb, query, "manta"))
        self.assertEqual(linear[0][1], "matching_failure")

    def test_bash_syntax(self):
        for p in ("5.Variant_benchmark.sh", "src/trio_mendelian_inconsistency.sh"):
            subprocess.run(["bash", "-n", str(ROOT / p)], check=True)

    def test_mendelian_script_preserves_reference_calls(self):
        text = (SRC / "trio_mendelian_inconsistency.sh").read_text()
        self.assertIn("bcftools merge -Oz", text)
        self.assertNotIn("bcftools merge -0", text)
        self.assertIn("rtg mendelian", text)
        self.assertIn('INFO/MCV', text)


if __name__ == "__main__":
    unittest.main()
