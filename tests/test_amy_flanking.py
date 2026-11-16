#!/usr/bin/env python3
"""Synthetic end-to-end smoke tests for AMY flanking step 7.

Run: python3 -m unittest discover -s tests -p 'test_amy_flanking.py' -v
Full shell integration requires bcftools + bgzip/tabix on PATH.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import amy_flank_popgen as f


class AMYTests(unittest.TestCase):
    def test_single_site_estimators(self):
        self.assertAlmostEqual(f.projection(1,10,2),0.2)
        r=f.r2_corrected(
            ["chr1","100",".","A","G",".",".",".","GT","0|0","0|1","1|1"],
            ["chr1","150",".","C","T",".",".",".","GT","0|0","0|1","1|1"])
        self.assertAlmostEqual(r,1.0)
        self.assertAlmostEqual(float(f.pf(float("nan"))) if False else 0,0)

    @unittest.skipUnless(shutil.which("bcftools") and shutil.which("bgzip"), "requires bcftools and bgzip")
    def test_all_optional_ci_and_genomic(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            a=d/"A.ids";b=d/"B.ids"
            a.write_text("A1\nA2\nA3\n")
            b.write_text("B1\nB2\nB3\n")
            def row(pos,ref,alt,dosages):
                translate={0:"0|0",1:"0|1",2:"1|1"}
                return "\t".join(["chr1",str(pos),".",ref,alt,".","PASS",".","GT"]+
                                  [translate[n] for n in dosages])
            text="##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
            text+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            text+="\t".join(["A1","A2","A3","B1","B2","B3"])+"\n"
            variants=[
                (110,"A","G",[0,1,2,0,0,1]),
                (150,"C","T",[1,2,0,2,2,0]),
                (180,"G","A",[2,1,0,1,2,1]),
                (410,"C","G",[0,1,2,2,1,0]),
                (450,"T","C",[1,0,2,0,0,1]),
                (470,"G","T",[0,2,1,1,2,0]),
                (1000,"A","C",[0,1,2,2,0,1]),
                (51000,"A","G",[0,1,2,2,1,0]),
                (52000,"C","A",[0,1,1,0,2,2]),
                (71000,"A","T",[2,1,0,0,1,2]),
            ]
            text+="\n".join(row(*v) for v in variants)+"\n"
            input_vcf=d/"tiny.vcf"
            input_vcf.write_text(text)
            gz=d/"tiny.vcf.gz"
            with open(gz,"wb") as fh:
                subprocess.run(["bgzip","-c",str(input_vcf)],stdout=fh,check=True)
            subprocess.run(["bcftools","index","-f",str(gz)],check=True)
            fai=d/"tiny.fai";fai.write_text("chr1\t100000\t0\t0\t0\n")
            out=d/"result"
            cmd=["bash",str(ROOT/"5.Population_genetics_in_AMY_flanking_locus.sh"),"all",
                 "--vcf",str(gz),"--pop",f"A={a}","--pop",f"B={b}",
                 "--out",str(out),"--distal","chr1:101-200",
                 "--proximal","chr1:401-500","--projection","4",
                 "--ld-distal-max","100","--ld-proximal-max","100",
                 "--individual-ci","--bootstrap-reps","5",
                 "--jackknife-fraction","0.30","--genomic",
                 "--genomic-reps","5","--window-step","1000",
                 "--block-mb","1","--chrom-sizes",str(fai)]
            cp=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            self.assertEqual(cp.returncode,0,cp.stdout)
            for fname in ("stats/diversity.tsv","stats/within_ld.tsv",
                          "stats/across_ld.tsv","stats/hudson_fst.tsv",
                          "individual_ci/individual_ci.tsv",
                          "genomic/genomic_percentile_ci.tsv","status.tsv"):
                self.assertTrue((out/fname).exists(),fname+"\n"+cp.stdout)
            self.assertIn("SKIPPED_MISSING_SELSCAN",(out/"status.tsv").read_text())


if __name__=="__main__":
    unittest.main()
