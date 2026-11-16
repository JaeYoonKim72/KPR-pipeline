#!/usr/bin/env python3
"""Merge duplicate records of one sample (same CHROM/POS/REF/ALT) into one record.

After per-sample vcfwave + norm -m-any + vcfdecompose, an allele carried by both haplotypes can be
written as two records (e.g. 1|0 and 0|1). Without merging, a homozygous variant would be counted
as two heterozygous ones. Haplotype i of the merged GT = 1 if any record has the ALT there,
0 if any record has REF there, '.' otherwise. Other fields are taken from the first record.

Input: sorted single-sample VCF (text) on stdin. Output: VCF on stdout.
--report FILE: number of input records, output records and merged duplicates.
"""
import argparse
import re
import sys


def merge_gt(gts):
    """gts: list of GT strings of the same allele -> merged GT string"""
    parsed = [re.split(r'([|/])', g) for g in gts]
    n = max(len(p) for p in parsed)
    sep = '|' if all('/' not in g for g in gts) else '/'
    out = []
    for i in range(0, n, 2):
        vals = [p[i] for p in parsed if i < len(p)]
        out.append('1' if '1' in vals else ('0' if '0' in vals else '.'))
    return sep.join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report')
    a = p.parse_args()
    n_in = n_out = 0
    group = {}          # (ref, alt) -> [fields, [GTs]] for the current CHROM:POS
    order = []
    cur = None

    def flush():
        nonlocal n_out
        for key in order:
            f, gts = group[key]
            if len(gts) > 1:
                fmt = f[8].split(':')
                gi = fmt.index('GT')
                sd = f[9].split(':')
                sd[gi] = merge_gt(gts)
                f[9] = ':'.join(sd)
            sys.stdout.write('\t'.join(f) + '\n')
            n_out += 1
        group.clear()
        order.clear()

    for line in sys.stdin:
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
        n_in += 1
        f = line.rstrip('\n').split('\t')
        if len(f) != 10:
            sys.exit(f'ERROR: expected 1 sample column, got {len(f) - 9}')
        pos = (f[0], f[1])
        if pos != cur:
            flush()
            cur = pos
        gt = f[9].split(':')[f[8].split(':').index('GT')]
        key = (f[3], f[4])
        if key in group:
            group[key][1].append(gt)
        else:
            group[key] = [f, [gt]]
            order.append(key)
    flush()
    if a.report:
        with open(a.report, 'w') as r:
            r.write(f'records_in\t{n_in}\nrecords_out\t{n_out}\nduplicates_merged\t{n_in - n_out}\n')


if __name__ == '__main__':
    main()
