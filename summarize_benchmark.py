#!/usr/bin/env python3
"""Summarize the hap.py and Truvari results (pangenome_variant_benchmarking.sh, step summary).

Reads <outdir>/bench/small/<sample>/<dataset>.<type>.summary.csv and
<outdir>/bench/sv/<sample>/<dataset>/<level>/<type>/summary.json, and writes to <outdir>/summary/:
  kor_small_detail.tsv   per sample x dataset x type (SNP, DEL, INS), hap.py Filter=ALL
  kor_sv_detail.tsv      per sample x dataset x support level x type (DEL, INS), Truvari
  kor_small_mean.tsv     mean / SD over samples
  kor_sv_mean.tsv
  kor_wilcoxon.tsv       paired two-sided exact Wilcoxon signed-rank tests between graphs
"""
import argparse
import csv
import itertools
import json
import os
import statistics

METRICS = ('Precision', 'Recall', 'F1')

def happy_rows(path, filt='ALL'):
    with open(path) as f:
        return {r['Type']: r for r in csv.DictReader(f) if r['Filter'] == filt}

def happy_metrics(r):
    return dict(Truth=int(r['TRUTH.TOTAL']), Query=int(r['QUERY.TOTAL']), TP=int(r['TRUTH.TP']),
                FN=int(r['TRUTH.FN']), FP=int(r['QUERY.FP']), Precision=float(r['METRIC.Precision']),
                Recall=float(r['METRIC.Recall']), F1=float(r['METRIC.F1_Score']))

def truvari_metrics(path):
    d = json.load(open(path))
    return dict(Truth=d['TP-base'] + d['FN'], Query=d['TP-comp'] + d['FP'], TP=d['TP-base'],
                TP_comp=d['TP-comp'], FN=d['FN'], FP=d['FP'], Precision=d['precision'] or 0.0,
                Recall=d['recall'] or 0.0, F1=d['f1'] or 0.0)

def write(path, rows):
    if not rows:
        return
    cols = list(rows[0])
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.6f}' if isinstance(v, float) else v) for k, v in r.items()})

def mean_table(rows, keys):
    out = []
    groups = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in keys), []).append(r)
    for g, rs in groups.items():
        d = dict(zip(keys, g), N=len(rs))
        for k in ('Truth', 'Query', 'TP', 'FN', 'FP'):
            d[k] = statistics.mean(r[k] for r in rs)
        for m in METRICS:
            v = [r[m] for r in rs]
            d[m + '_mean'] = statistics.mean(v)
            d[m + '_sd'] = statistics.stdev(v) if len(v) > 1 else 0.0
        out.append(d)
    return out

def signed_rank_exact(diff):
    """Two-sided exact Wilcoxon signed-rank p-value (zero differences dropped, mid-ranks for ties)."""
    d = [x for x in diff if x != 0]
    n = len(d)
    if n == 0:
        return 1.0
    a = sorted(abs(x) for x in d)
    rank = {v: statistics.mean(i + 1 for i, x in enumerate(a) if x == v) for v in set(a)}
    r = [rank[abs(x)] for x in d]
    mu = sum(r) / 2
    obs = abs(sum(ri for ri, x in zip(r, d) if x > 0) - mu)
    hits = sum(abs(sum(ri for ri, s in zip(r, sg) if s) - mu) >= obs - 1e-9
               for sg in itertools.product((0, 1), repeat=n))
    return hits / 2 ** n

def wilcoxon(rows, keys, graphs, samples):
    out = []
    idx = {tuple(r[k] for k in keys) + (r['Dataset'], r['Sample']): r for r in rows}
    groups = sorted({tuple(r[k] for k in keys) for r in rows})
    for g in groups:
        for a, b in itertools.combinations(graphs, 2):
            for m in METRICS:
                try:
                    x = [idx[g + (a, s)][m] for s in samples]
                    y = [idx[g + (b, s)][m] for s in samples]
                except KeyError:
                    continue
                d = [i - j for i, j in zip(x, y)]
                out.append(dict(zip(keys, g), A=a, B=b, Metric=m, Median_A=statistics.median(x),
                                Median_B=statistics.median(y), A_higher=sum(v > 0 for v in d),
                                B_higher=sum(v < 0 for v in d), P_value=signed_rank_exact(d)))
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--graphs', nargs='+', default=['kpr', 'hprc_kpr', 'hprc'])
    a = ap.parse_args()
    bench, summ = os.path.join(a.outdir, 'bench'), os.path.join(a.outdir, 'summary')
    os.makedirs(summ, exist_ok=True)

    small = []
    d = os.path.join(bench, 'small')
    for s in a.samples:
        sd = os.path.join(d, s)
        if not os.path.isdir(sd):
            continue
        for f in sorted(os.listdir(sd)):
            if not f.endswith('.summary.csv'):
                continue
            ds, typ = f[:-len('.summary.csv')].rsplit('.', 1)
            r = happy_rows(os.path.join(sd, f))['SNP' if typ == 'SNP' else 'INDEL']
            small.append(dict(Sample=s, Dataset=ds, Type=typ, **happy_metrics(r)))
    write(os.path.join(summ, 'kor_small_detail.tsv'), small)
    write(os.path.join(summ, 'kor_small_mean.tsv'), mean_table(small, ('Type', 'Dataset')))

    sv = []
    d = os.path.join(bench, 'sv')
    for s in a.samples:
        for root, _, files in os.walk(os.path.join(d, s)):
            if 'summary.json' in files:
                ds, lv, typ = os.path.relpath(root, os.path.join(d, s)).split(os.sep)
                sv.append(dict(Sample=s, Dataset=ds, Support=lv, Type=typ,
                               **truvari_metrics(os.path.join(root, 'summary.json'))))
    write(os.path.join(summ, 'kor_sv_detail.tsv'), sv)
    write(os.path.join(summ, 'kor_sv_mean.tsv'), mean_table(sv, ('Support', 'Type', 'Dataset')))

    tests = [dict(Variant='small', Support='-', **t) for t in wilcoxon(small, ('Type',), a.graphs, a.samples)]
    tests += [dict(Variant='SV', **t) for t in wilcoxon(sv, ('Support', 'Type'), a.graphs, a.samples)]
    write(os.path.join(summ, 'kor_wilcoxon.tsv'), tests)
    print('[SUMMARY]', summ)

if __name__ == '__main__':
    main()
