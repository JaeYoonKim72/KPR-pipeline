#!/usr/bin/env python3
"""Helpers for building the Korean SV truth sets (pangenome_variant_benchmarking.sh,
steps sv_truth_filter, sv_truth_support and sv_truth_readsupport).

Sub-commands
  filter     Keep PASS INS/DEL >= 50 bp on chr1-22/X/Y with a complete non-reference GT.
             Read-based callers additionally need read support:
               pbsv     : DP >= 7, ALT reads (AD[1]) >= 3
               SVIM     : DP >= 7, ALT reads (AD[1]) >= 3, QUAL >= 10
               Sniffles : GQ >= 20, DR + DV >= 7, DV >= 3
             SVIM-asm records get stable IDs (<sample>_ASM_00000001, ...).
  collect    For each SVIM-asm SV, record which read-based callers matched it
             (Truvari tp-base of bench_<caller>), and write <sample>.SV.support.tsv.
  levels     Write <sample>.SV.truth.ge1/ge2/ge3.vcf.gz: SVIM-asm SVs supported by
             >= 1, >= 2 or all 3 read-based callers (assembly genotypes are kept).
  fixheader  Set the INFO/SVLEN header to Number=A (or add it) without touching records.
"""
import argparse
import collections
import csv
import gzip
import re
import sys

import pysam

CHROM_RE = re.compile(r'^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$')
READ_CALLERS = ('pbsv', 'sniffles', 'svim')
LEVELS = {'ge1': 1, 'ge2': 2, 'ge3': 3}


def is_alt_complete(gt):
    return gt is not None and len(gt) > 0 and all(a is not None for a in gt) and any(a > 0 for a in gt)


def size(rec):
    val = rec.info.get('SVLEN')
    if isinstance(val, (tuple, list)):
        val = val[0] if val else None
    try:
        return abs(int(val))
    except (TypeError, ValueError):
        if rec.info.get('SVTYPE') == 'DEL' and rec.stop > rec.start:
            return rec.stop - rec.start
        return 0


def fmt_number(sample, key, index=None):
    if key not in sample:
        return None
    v = sample.get(key)
    if index is not None:
        if not isinstance(v, (tuple, list)) or len(v) <= index:
            return None
        v = v[index]
    elif isinstance(v, (tuple, list)):
        v = v[0] if v else None
    try:
        return float(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def passes(rec, caller, sample):
    if not CHROM_RE.fullmatch(rec.contig):
        return False
    if rec.filter.keys() and 'PASS' not in rec.filter.keys():
        return False
    if rec.info.get('SVTYPE') not in ('INS', 'DEL') or size(rec) < 50:
        return False
    if not is_alt_complete(sample.get('GT')):
        return False
    if caller == 'svim_asm':
        return True
    if caller in ('pbsv', 'svim'):
        dp, ad_alt = fmt_number(sample, 'DP'), fmt_number(sample, 'AD', 1)
        if dp is None or ad_alt is None or dp < 7 or ad_alt < 3:
            return False
        if caller == 'svim' and (rec.qual is None or rec.qual < 10):
            return False
    if caller == 'sniffles':
        gq, dr, dv = (fmt_number(sample, k) for k in ('GQ', 'DR', 'DV'))
        if any(x is None for x in (gq, dr, dv)) or gq < 20 or dr + dv < 7 or dv < 3:
            return False
    return True


def cmd_filter(a):
    inv = pysam.VariantFile(a.input)
    if len(inv.header.samples) != 1:
        raise ValueError(f'Expected a single-sample VCF: {a.input}')
    name = next(iter(inv.header.samples))
    n = 0
    with pysam.VariantFile(a.output, 'wz', header=inv.header.copy()) as out:
        for rec in inv:
            if not passes(rec, a.caller, rec.samples[name]):
                continue
            n += 1
            if a.caller == 'svim_asm':
                rec.id = f'{a.sample}_ASM_{n:08d}'
            out.write(rec)
    print(f'[FILTER] {a.sample} {a.caller}: {n} SVs')


def cmd_collect(a):
    d, s = a.dir, a.sample
    with pysam.VariantFile(f'{d}/{s}.svim_asm.alt.vcf.gz') as f:
        assembly = list(f)
    ids = {r.id for r in assembly}
    if len(ids) != len(assembly):
        raise ValueError('SVIM-asm IDs are not unique')
    support = collections.defaultdict(set)
    for c in READ_CALLERS:
        with pysam.VariantFile(f'{d}/bench_{c}/tp-base.vcf.gz') as f:
            tp = {r.id for r in f}
        unknown = tp - ids
        if unknown:
            raise ValueError(f'{c}: {len(unknown)} tp-base IDs not in the assembly VCF')
        for i in tp:
            support[i].add(c)
    with open(f'{d}/{s}.SV.support.tsv', 'w') as fh:
        fh.write('ASM_ID\tCHROM\tPOS\tSVTYPE\tSVLEN\tSupport_count\tSupporting_callers\n')
        for r in assembly:
            c = sorted(support.get(r.id, ()))
            fh.write(f'{r.id}\t{r.contig}\t{r.pos}\t{r.info.get("SVTYPE")}\t{size(r)}\t{len(c)}\t{",".join(c) or "."}\n')
    print(f'[COLLECT] {s}: {len(assembly)} SVIM-asm SVs, {len(support)} supported by >= 1 caller')


def cmd_levels(a):
    d, s = a.dir, a.sample
    with open(f'{d}/{s}.SV.support.tsv') as fh:
        n_sup = {row['ASM_ID']: int(row['Support_count']) for row in csv.DictReader(fh, delimiter='\t')}
    with pysam.VariantFile(f'{d}/{s}.svim_asm.alt.vcf.gz') as src:
        paths = {k: f'{d}/{s}.SV.truth.{k}.vcf.gz' for k in LEVELS}
        outs = {k: pysam.VariantFile(p, 'wz', header=src.header.copy()) for k, p in paths.items()}
        counts = collections.Counter()
        for r in src:
            for k, m in LEVELS.items():
                if n_sup[r.id] >= m:
                    outs[k].write(r)
                    counts[k] += 1
        for o in outs.values():
            o.close()
    for p in paths.values():
        pysam.tabix_index(p, preset='vcf', force=True)
    print(f'[LEVELS] {s}: ' + ', '.join(f'{k}={counts[k]}' for k in LEVELS))


def cmd_fixheader(a):
    op = gzip.open if a.input.endswith('.gz') else open
    seen = False
    with op(a.input, 'rt') as fin, open(a.output, 'wt') as fout:
        for line in fin:
            if line.startswith('##INFO=<ID=SVLEN,'):
                line = re.sub(r'Number=[^,>]+', 'Number=A', line, count=1)
                seen = True
            elif line.startswith('#CHROM') and not seen:
                fout.write('##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length for each ALT allele">\n')
                seen = True
            fout.write(line)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    x = sub.add_parser('filter')
    for k in ('input', 'output', 'caller', 'sample'):
        x.add_argument('--' + k, required=True)
    for name in ('collect', 'levels'):
        x = sub.add_parser(name)
        x.add_argument('--sample', required=True)
        x.add_argument('--dir', required=True)
    x = sub.add_parser('fixheader')
    x.add_argument('--input', required=True)
    x.add_argument('--output', required=True)
    a = p.parse_args()
    {'filter': cmd_filter, 'collect': cmd_collect, 'levels': cmd_levels, 'fixheader': cmd_fixheader}[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
