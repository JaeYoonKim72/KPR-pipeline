"""Static and help-mode tests for ensemble annotation stage 2."""

import ast
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[2]
SHELL = ROOT / "2.Ensemble_gene_annotation.sh"
DOCUMENT = ROOT / "src" / "Ensemble_gene_annotation.md"
FIGURE = ROOT / "gene_annotation_workflow.png"


def embedded_python(source):
    start = '  "${PYTHON_BIN}" - "$@" "${ENGINE_ARGS[@]}" <<' + "'PY'\n"
    begin = source.find(start)
    if begin < 0:
        raise ValueError("Embedded Python delimiter missing")
    begin += len(start)
    end = source.find("\nPY\n}", begin)
    if end < 0:
        raise ValueError("Embedded Python closing marker missing")
    return source[begin:end]


class GeneAnnotationWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.script = SHELL.read_text(encoding="utf-8")
        cls.doc = DOCUMENT.read_text(encoding="utf-8")

    def test_help(self):
        p = subprocess.run(
            ["bash", str(SHELL), "--help"],
            capture_output=True, text=True, check=True
        )
        for step in (
            "softmask", "predict_orfs", "protein_support", "gene_annotation",
            "compare_annotations", "search_unannotated", "filter_flagger",
            "tables", "check", "all"
        ):
            self.assertIn(step, p.stdout)
        self.assertIn("2.Ensemble_gene_annotation.sh", p.stdout)

    def test_shell_syntax(self):
        subprocess.run(["bash", "-n", str(SHELL)], check=True)

    def test_inline_python_syntax(self):
        ast.parse(embedded_python(self.script), filename="gene_annotation_engine.py")

    def test_separate_protein_thresholds(self):
        self.assertIn("MIN_PROTEIN_IDENTITY=${MIN_PROTEIN_IDENTITY:-95}", self.script)
        self.assertIn("MIN_BLASTP_IDENTITY=${MIN_BLASTP_IDENTITY:-90}", self.script)
        self.assertIn("if aa['identity']<self.a.protein_identity", self.script)
        self.assertIn("if best['identity']<self.a.blastp_identity", self.script)
        self.assertIn("p.add_argument('--blastp-identity',type=float,default=90", self.script)

    def test_manuscript_commands(self):
        for snippet in (
            "--softmasking=1", "'--percent','80'", "'-copies','-sc','0.95'",
            "'-polish','-exclude_partial'", "ic=chain(parts,t['cds_parts'])",
            "FULLY_HAP",
        ):
            self.assertIn(snippet, self.script)

    def test_docs_and_figure(self):
        for heading in (
            "## 1. Set up", "## 2. Run", "## 3. Steps",
            "## 4. Configuration", "## 5. Output", "## 6. Analysis procedures",
        ):
            self.assertIn(heading, self.doc)
        self.assertIn("../gene_annotation_workflow.png", self.doc)
        self.assertTrue(FIGURE.is_file())
        with FIGURE.open("rb") as f:
            self.assertEqual(f.read(8), b"\x89PNG\r\n\x1a\n")
        self.assertIn("≥90%", self.doc)
        self.assertIn("≥95%", self.doc)

    def test_placeholder_removed(self):
        self.assertNotIn("NOT IMPLEMENTED", self.script)
        self.assertNotIn("run_test.sh", self.script)
        self.assertNotIn("RUNNING_gene_annotation.md", self.script)
        self.assertNotIn("Status: placeholder only", self.doc)


if __name__ == "__main__":
    unittest.main()
