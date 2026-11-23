#!/usr/bin/env python3
"""Read-support filter of the Korean SV truth sets (5.Variant_benchmark.sh, step sv_truth_readsupport).

Each truth SV is matched to the raw read-based calls (Sniffles, pbsv, SVIM) of the same
sample: same SVTYPE, start position within --dist bp, size similarity >= --pctsize and,
for insertions, inserted-sequence similarity >= --pctseq. For every caller the best match
is kept and its depth / ALT-read support is recorded:
    Sniffles : DP = DR + DV, ALT = DV (or INFO/SUPPORT)
    pbsv     : DP = FORMAT/DP, ALT = AD[1]
    SVIM     : DP = FORMAT/DP, ALT = AD[1] (or INFO/SUPPORT)
A truth SV is kept only when at least one matched caller has ALT >= --min-alt AND
DP >= --min-dp (both from the same caller). SVs without such a call, including SVs for
which no matched caller reports these values, are removed.

Outputs (in --outdir):
    truth.readsupport.vcf.gz  filtered truth (INFO: EVCALL, EVMAXALT, EVMAXDP, EVSTATUS)
    evidence.tsv              per-SV evidence
    filter_counts.tsv         kept / total / unmatched
"""
import argparse
import collections
import csv
import difflib
import os
import subprocess

import pysam

def svtype(r):
    try:
        return str(r.info['SVTYPE'])
    except (KeyError, ValueError):
        if r.alts and not r.alts[0].startswith('<'):
            if len(r.alts[0]) > len(r.ref):
                return 'INS'
            if len(r.alts[0]) < len(r.ref):
                return 'DEL'
        return 'OTHER'

def length(r):
    try:
        v = r.info['SVLEN']
        v = v[0] if isinstance(v, tuple) else v
        if v is not None:
            return abs(int(v))
    except (KeyError, ValueError, TypeError):
        pass
    if r.alts and len(r.alts) == 1 and not r.alts[0].startswith('<'):
        return abs(len(r.alts[0]) - len(r.ref))
    return None

def insseq(r):
    if svtype(r) != 'INS' or not r.alts or len(r.alts) != 1 or r.alts[0].startswith('<'):
        return None
    a, b = r.alts[0], r.ref
    return a[len(b):] if a.startswith(b) else a[1:] if len(a) > 1 else None

def sample_value(r, k):
    try:
        return next(iter(r.samples.values())).get(k)
    except (KeyError, ValueError, StopIteration):
        return None

def integer(v):
    try:
        return int(v) if v is not None else None
    except (TypeError, ValueError):
        return None

def evidence(r, caller):
    """Return (depth, alt_reads) of a read-based call."""
    if caller == 'sniffles':
        dv, dr = integer(sample_value(r, 'DV')), integer(sample_value(r, 'DR'))
        if dv is None:
            try:
                dv = integer(r.info['SUPPORT'])
            except (KeyError, ValueError):
                pass
        return (dr + dv if dr is not None and dv is not None else None), dv
    dp, ad = integer(sample_value(r, 'DP')), sample_value(r, 'AD')
    alt = integer(ad[1]) if ad is not None and len(ad) > 1 else None
    if caller == 'svim' and alt is None:
        try:
            alt = integer(r.info['SUPPORT'])
        except (KeyError, ValueError):
            pass
    return dp, alt

def candidates(truth, vcf, caller, dist, pctsize, pctseq):
    lt = length(truth)
    if not lt:
        return []
    try:
        it = vcf.fetch(truth.contig, max(0, truth.pos - 1 - dist), truth.pos + dist)
    except (ValueError, KeyError):
        return []
    seq, hits = insseq(truth), []
    for r in it:
        if svtype(r) != svtype(truth):
            continue
        lr = length(r)
        if not lr or min(lt, lr) / max(lt, lr) < pctsize or abs(r.pos - truth.pos) > dist:
            continue
        sim = None
        if seq is not None:
            other = insseq(r)
            if other is not None and max(len(seq), len(other)) <= 2000:
                sim = difflib.SequenceMatcher(None, seq, other, autojunk=False).ratio()
                if sim < pctseq:
                    continue
        dp, alt = evidence(r, caller)
        rank = (1 if sim is not None else 0, sim or 0, min(lt, lr) / max(lt, lr), -abs(r.pos - truth.pos), alt or -1)
        hits.append((rank, dp, alt))
    return hits

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ('truth', 'sniffles', 'pbsv', 'svim', 'outdir'):
        ap.add_argument('--' + k, required=True)
    ap.add_argument('--dist', type=int, default=500)
    ap.add_argument('--pctsize', type=float, default=0.7)
    ap.add_argument('--pctseq', type=float, default=0.7)
    ap.add_argument('--support-tsv', required=True)
    ap.add_argument('--min-alt', type=int, default=5)
    ap.add_argument('--min-dp', type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    with open(a.support_tsv, newline='') as support_file:
        support = {r['ASM_ID']: set() if r['Supporting_callers'] == '.' else set(r['Supporting_callers'].split(',')) for r in csv.DictReader(support_file, delimiter='\t')}
    vf = pysam.VariantFile(a.truth)
    callers = {c: pysam.VariantFile(getattr(a, c)) for c in ('sniffles', 'pbsv', 'svim')}
    hdr = vf.header.copy()
    hdr.add_line('##INFO=<ID=EVCALL,Number=.,Type=String,Description="Matched read-based callers">')
    hdr.add_line('##INFO=<ID=EVMAXALT,Number=1,Type=Integer,Description="Maximum ALT-read support among matched callers">')
    hdr.add_line('##INFO=<ID=EVMAXDP,Number=1,Type=Integer,Description="Maximum depth among matched callers">')
    hdr.add_line('##INFO=<ID=EVSTATUS,Number=1,Type=String,Description="MATCHED or UNMATCHED">')
    out_vcf = os.path.join(a.outdir, 'truth.readsupport.vcf.gz')
    out = pysam.VariantFile(out_vcf, 'wz', header=hdr)
    counts = collections.Counter()
    with open(os.path.join(a.outdir, 'evidence.tsv'), 'w') as qc:
        qc.write('chrom\tpos\tid\ttype\tlength\tmatched_callers\tmax_alt\tmax_dp\t'
                 'sniffles_alt\tsniffles_dp\tpbsv_alt\tpbsv_dp\tsvim_alt\tsvim_dp\tkept\n')
        for r in vf:
            typ = svtype(r)
            if typ not in ('INS', 'DEL'):
                continue
            rr = r.copy()
            rr.translate(hdr)
            info, names = {}, []
            for c, v in callers.items():
                hits = candidates(r, v, c, a.dist, a.pctsize, a.pctseq) if c in support.get(r.id, set()) else []
                if hits:
                    best = max(hits, key=lambda h: h[0])
                    info[c] = (best[1], best[2])
                    names.append(c)
                else:
                    info[c] = (None, None)
            alts = [x[1] for x in info.values() if x[1] is not None]
            dps = [x[0] for x in info.values() if x[0] is not None]
            ma, md = (max(alts) if alts else None), (max(dps) if dps else None)
            rr.info['EVSTATUS'] = 'MATCHED' if names else 'UNMATCHED'
            if names:
                rr.info['EVCALL'] = tuple(names)
            if ma is not None:
                rr.info['EVMAXALT'] = ma
            if md is not None:
                rr.info['EVMAXDP'] = md
            observed = [(dp, alt) for dp, alt in info.values() if dp is not None and alt is not None]
            keep = any(alt >= a.min_alt and dp >= a.min_dp for dp, alt in observed)
            counts['total'] += 1
            counts['unmatched'] += not names
            if keep:
                out.write(rr)
                counts['kept'] += 1
            fields = [r.contig, r.pos, r.id or '.', typ, length(r) or '.', ','.join(names) or '.',
                      ma if ma is not None else '.', md if md is not None else '.']
            fields += [v if v is not None else '.' for c in ('sniffles', 'pbsv', 'svim') for v in (info[c][1], info[c][0])]
            qc.write('\t'.join(map(str, fields + [int(keep)])) + '\n')
    out.close()
    subprocess.run(['bcftools', 'index', '-f', '-t', out_vcf], check=True)
    with open(os.path.join(a.outdir, 'filter_counts.tsv'), 'w') as f:
        f.write('kept\ttotal\tunmatched\n')
        f.write(f'{counts["kept"]}\t{counts["total"]}\t{counts["unmatched"]}\n')
    print(f'total {counts["total"]} unmatched {counts["unmatched"]} kept {counts["kept"]}')

if __name__ == '__main__':
    main()
