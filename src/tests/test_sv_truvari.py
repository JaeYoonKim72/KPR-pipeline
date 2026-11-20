"""Tests for SV caller filtering and Truvari matched-base intersections."""

import sys
import tempfile
import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SRC))
from filter_sv_caller_quality import filter_record, process
from select_svim_asm_truvari import select

HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def row(pos, length=100):
    return f"chr1\t{pos}\t.\tA\t<DEL>\t30\tPASS\tSVLEN=-{length};SVTYPE=DEL\n"


def record(caller, fmt="GT:DP:AD", sample="0/1:10:5,5", qual="30", info="SVLEN=-100;SVTYPE=DEL"):
    data = ["chr1", "100", ".", "A", "<DEL>", qual, "PASS", info]
    if caller != "svimasm":
        data.extend([fmt, sample])
    return data


class QualityTests(unittest.TestCase):
    def test_pbsv_depth_alt(self):
        self.assertIsNone(filter_record(record("pbsv"), "pbsv", 7, 3, 20, 10))
        self.assertEqual(filter_record(record("pbsv", sample="0/1:6:3,3"), "pbsv", 7, 3, 20, 10), "depth_below_threshold")
        self.assertEqual(filter_record(record("pbsv", sample="0/1:10:8,2"), "pbsv", 7, 3, 20, 10), "alt_support_below_threshold")

    def test_sniffles_gq_dr_dv(self):
        v = record("sniffles", fmt="GT:GQ:DR:DV", sample="0/1:20:4:3")
        self.assertIsNone(filter_record(v, "sniffles", 7, 3, 20, 10))
        v[-1] = "0/1:19:4:3"
        self.assertEqual(filter_record(v, "sniffles", 7, 3, 20, 10), "genotype_quality_below_threshold")
        v[-1] = "0/1:20:3:3"
        self.assertEqual(filter_record(v, "sniffles", 7, 3, 20, 10), "depth_below_threshold")

    def test_svim_qual_dp_support(self):
        v = record("svim", fmt="GT:DP", sample="0/1:7", qual="10",
                   info="SVLEN=-100;SVTYPE=DEL;SUPPORT=3")
        self.assertIsNone(filter_record(v, "svim", 7, 3, 20, 10))
        v[5] = "9"
        self.assertEqual(filter_record(v, "svim", 7, 3, 20, 10), "svim_qual_below_threshold")
        v[5] = "10"
        v[7] = "SVLEN=-100;SVTYPE=DEL"
        self.assertEqual(filter_record(v, "svim", 7, 3, 20, 10), "missing_alt_support")

    def test_svimasm_pass_autosome_sv_length(self):
        v = record("svimasm", info="SVLEN=-51;SVTYPE=DEL")
        self.assertIsNone(filter_record(v, "svimasm", 7, 3, 20, 10))
        v[7] = "SVLEN=-50;SVTYPE=DEL"
        self.assertEqual(filter_record(v, "svimasm", 7, 3, 20, 10), "sv_length_below_51")
        v[7] = "SVLEN=-51;SVTYPE=DEL"
        v[0] = "chrX"
        self.assertEqual(filter_record(v, "svimasm", 7, 3, 20, 10), "non_autosomal")
        v[0] = "chr1"
        v[6] = "."
        self.assertEqual(filter_record(v, "svimasm", 7, 3, 20, 10), "filter_not_pass")

    def test_missing_fields_excluded_and_audited(self):
        with tempfile.TemporaryDirectory() as td:
            t = Path(td)
            (t / "raw.vcf").write_text(HEAD.rstrip() + "\tFORMAT\tSAMPLE\n" + "\t".join(record("pbsv", "GT:DP", "0/1:10")) + "\n")
            counts = process(str(t / "raw.vcf"), str(t / "pass.vcf"), str(t / "audit.tsv"), "pbsv", 7, 3, 20, 10)
            self.assertEqual(counts["missing_alt_support"], 1)
            self.assertEqual(counts["retained"], 0)
            self.assertIn("missing_alt_support\t1", (t / "audit.tsv").read_text())


class AggregateTests(unittest.TestCase):
    def setUp(self):
        t = tempfile.TemporaryDirectory()
        self.addCleanup(t.cleanup)
        self.t = Path(t.name)
        self.base = self.t / "svimasm.vcf"
        self.base.write_text(HEAD + "".join(row(p) for p in (100, 200, 300, 400)))
        self.calls = {}
        for caller, positions in [("pbsv", (100, 200)), ("sniffles", (100, 300)), ("svim", (100, 200, 400))]:
            fn = self.t / f"{caller}.vcf"
            fn.write_text(HEAD + "".join(row(p) for p in positions))
            self.calls[caller] = str(fn)

    def run_select(self, n, strong=None):
        out = self.t / f"final_{n}.vcf"
        audit = self.t / f"audit_{n}.tsv"
        selected = select(str(self.base), self.calls, str(out), str(audit), n, strong)
        positions = [int(x.split("\t")[1]) for x in out.read_text().splitlines() if not x.startswith("#")]
        return selected, positions, audit.read_text()

    def test_primary_two_match(self):
        count, positions, audit = self.run_select(2)
        self.assertEqual((count, positions), (2, [100, 200]))
        self.assertIn("matched_read_callers_3\t1", audit)

    def test_truthset_one_two_three(self):
        self.assertEqual(self.run_select(1)[1], [100, 200, 300, 400])
        self.assertEqual(self.run_select(2)[1], [100, 200])
        self.assertEqual(self.run_select(3)[1], [100])

    def test_strong_support_option(self):
        strong = {}
        for caller, positions in [("pbsv", (100,)), ("sniffles", ()), ("svim", ())]:
            fn = self.t / f"{caller}.high.vcf"
            fn.write_text(HEAD + "".join(row(p) for p in positions))
            strong[caller] = str(fn)
        self.assertEqual(self.run_select(2, strong)[1], [100])

    def test_missing_required_caller_rejected(self):
        with self.assertRaises(ValueError):
            select(str(self.base), {"pbsv": self.calls["pbsv"]}, str(self.t / "out"), str(self.t / "audit"), 2)

    def test_stage1_truvari_then_cohort_survivor(self):
        sh = (SRC.parent / "1.genome_assembly_and_statistics.sh").read_text()
        self.assertIn("truvari bench --base", sh)
        self.assertIn("--min-read-callers 2", sh)
        self.assertIn("cohort.filelist", sh)
        self.assertNotIn("survivor.min3", sh)
        self.assertNotIn("filter_svim_asm_support.py", sh)


if __name__ == "__main__":
    unittest.main()
