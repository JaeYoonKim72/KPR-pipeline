"""Synthetic tests for SVIM-asm-required structural-variant consensus."""

import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[2]
FILTER = ROOT / "src" / "filter_svim_asm_support.py"
PIPELINE = ROOT / "1.genome_assembly_and_statistics.sh"

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Caller support vector">\n'
    '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Supporting callers">\n'
    "##contig=<ID=chr1,length=100000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def record(pos, support):
    return (
        f"chr1\t{pos}\t.\tA\t<DEL>\t60\tPASS\t"
        f"SVTYPE=DEL;SUPP_VEC={support};SUPP={support.count('1')}\n"
    )


class SvimAsmConsensusTests(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.dir = pathlib.Path(self.work.name)

    def run_filter(self, vectors, caller_order=None):
        names = caller_order or ["pbsv.vcf", "sniffles.vcf", "svim.vcf", "svimasm.vcf"]
        filelist = self.dir / "filelist"
        filelist.write_text("".join(str(self.dir / name) + "\n" for name in names))
        vcf = self.dir / "candidate.vcf"
        vcf.write_text(VCF_HEADER + "".join(
            record(100 + 100 * i, vector) for i, vector in enumerate(vectors)
        ))
        output = self.dir / "final.vcf"
        audit = self.dir / "audit.tsv"
        result = subprocess.run(
            [
                sys.executable, str(FILTER),
                "--input", str(vcf),
                "--output", str(output),
                "--filelist", str(filelist),
                "--required-vcf", str(self.dir / "svimasm.vcf"),
                "--min-other-callers", "2",
                "--audit", str(audit),
            ],
            text=True,
            capture_output=True,
        )
        return result, output, audit

    def test_svim_asm_and_two_other_callers(self):
        vectors = ["1111", "1101", "1011", "0111", "1110", "1001", "0101", "0001", "0000"]
        result, output, audit = self.run_filter(vectors)
        self.assertEqual(result.returncode, 0, result.stderr)
        text = output.read_text()
        observed = [x.split("\t")[1] for x in text.splitlines() if not x.startswith("#")]
        self.assertEqual(observed, ["100", "200", "300", "400"])
        self.assertIn("SVTYPE=DEL;SUPP_VEC=1101;SUPP=3", text)
        with audit.open(newline="") as fh:
            counts = {row["metric"]: int(row["count"]) for row in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(counts["input_records"], 9)
        self.assertEqual(counts["retained_records"], 4)
        self.assertEqual(counts["excluded_no_svim_asm"], 2)
        self.assertEqual(counts["excluded_insufficient_other_callers"], 3)

    def test_required_caller_is_located_from_filelist_not_assumed_last(self):
        order = ["pbsv.vcf", "svimasm.vcf", "sniffles.vcf", "svim.vcf"]
        result, output, _ = self.run_filter(["1110", "1011", "0101"], order)
        self.assertEqual(result.returncode, 0, result.stderr)
        positions = [
            line.split("\t")[1]
            for line in output.read_text().splitlines() if not line.startswith("#")
        ]
        self.assertEqual(positions, ["100"])

    def test_invalid_support_vector_fails_closed(self):
        result, output, _ = self.run_filter(["110"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid SUPP_VEC", result.stderr)
        self.assertFalse(output.exists())

    def test_header_only_vcf_is_valid_with_zero_passes(self):
        result, output, audit = self.run_filter([])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(output.read_text().endswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
        self.assertIn("retained_records\t0\n", audit.read_text())

    def test_missing_required_caller_is_rejected(self):
        result, output, _ = self.run_filter(
            ["1111"], ["pbsv.vcf", "sniffles.vcf", "svim.vcf", "other.vcf"]
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("absent from filelist", result.stderr)
        self.assertFalse(output.exists())

    def test_shell_pipeline_uses_three_caller_survivor_and_support_filter(self):
        script = PIPELINE.read_text()
        self.assertIn("SURVIVOR merge", script)
        self.assertIn("100 3 1 0 0 51", script)
        self.assertIn("filter_svim_asm_support.py", script)
        self.assertIn("--required-vcf", script)
        self.assertIn("--min-other-callers 2", script)
        self.assertIn("cohort.filelist", script)


if __name__ == "__main__":
    unittest.main()
