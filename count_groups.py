#!/usr/bin/env python3
"""Count target-only / shared / other-only variants between two sample groups.

Reads `bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'` from stdin.
A group carries an allele if any haplotype of any of its samples has it
(missing '.' ignored).

--mode locus   : one record = one locus (multi-allelic VCF; SV after vcfcreatemulti
                 + PanSVMerger, so one ALT = one allele cluster). Reports
    locus                  : locus carried by target / other / both
    allele_in_record       : each ALT of each record
    allele_in_shared_locus : ALTs inside loci carried by both groups (APGp1 Fig. 6c text)
    shared_locus_with_target_private_allele
--mode overlap : sorted biallelic VCF (small indels); records whose REF spans overlap
                 form one locus (same locus definition as vcfcreatemulti for SVs),
                 alleles inside a locus are exact CHROM/POS/REF/ALT. Same units as locus.
--mode exact   : biallelic VCF; identical CHROM/POS/REF/ALT = same allele
    allele_exact
Output: TSV  varset unit class category count
"""
import argparse
import re
import sys
from collections import Counter

TARGET, OTHER = 1, 2
CATEGORY = {TARGET: 'target_only', OTHER: 'other_only', TARGET | OTHER: 'shared'}


def allele_class(ref, alt):
    lr, la = len(ref), len(alt)
    if la == lr:
        return 'SNV' if lr == 1 else 'OTHER'   # SNV only occurs in the small varset
    return 'INS' if la > lr else 'DEL'


def read_groups(path, target_pop):
    groups = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            sample, pop = line.split()[:2]
            groups[sample] = TARGET if pop == target_pop else OTHER
    return groups


def parse_header(line, groups):
    cols = line.rstrip('\n').lstrip('#').split('\t')
    flags = []
    for c in cols[4:]:
        name = re.sub(r'^\[\d+\]', '', c).rsplit(':GT', 1)[0]
        flags.append(groups.get(name, 0))
    unknown = [re.sub(r'^\[\d+\]', '', c) for c, g in zip(cols[4:], flags) if g == 0]
    if unknown:
        sys.stderr.write(f'WARN samples not in population table (ignored): {",".join(unknown)}\n')
    return flags


def carriers(gts, flags):
    """allele index -> group bit flags"""
    out = {}
    for gt, g in zip(gts, flags):
        if not g:
            continue
        for a in re.split(r'[|/]', gt):
            if a in ('.', '0', ''):
                continue
            a = int(a)
            out[a] = out.get(a, 0) | g
    return out


def count_locus(cnt, alleles):
    """alleles: list of (allele_class, group flags) of the carried alleles of one locus"""
    locus_g = 0
    for _, g in alleles:
        locus_g |= g
    if not locus_g:
        return
    classes = {ac for ac, _ in alleles}
    lc = classes.pop() if len(classes) == 1 else 'COMPLEX'
    cnt[('locus', lc, CATEGORY[locus_g])] += 1
    has_target_private = False
    for ac, g in alleles:
        cnt[('allele_in_record', ac, CATEGORY[g])] += 1
        if locus_g == TARGET | OTHER:
            cnt[('allele_in_shared_locus', ac, CATEGORY[g])] += 1
            has_target_private |= g == TARGET
    if has_target_private:
        cnt[('shared_locus_with_target_private_allele', lc, 'shared')] += 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mode', choices=['locus', 'overlap', 'exact'], required=True)
    p.add_argument('--varset', choices=['small', 'sv'], required=True)
    p.add_argument('--pop', required=True)
    p.add_argument('--target', required=True)
    a = p.parse_args()

    groups = read_groups(a.pop, a.target)
    cnt = Counter()
    exact = {}
    flags = None
    # overlap mode: current locus
    cur_chrom, cur_end, cur = None, -1, {}

    def flush():
        count_locus(cnt, [(allele_class(ref, al), g) for (_, ref, al), g in cur.items()])
        cur.clear()

    for line in sys.stdin:
        if line.startswith('#'):
            flags = parse_header(line, groups)
            continue
        f = line.rstrip('\n').split('\t')
        chrom, pos, ref, alt = f[:4]
        alts = alt.split(',')
        car = carriers(f[4:], flags)

        if a.mode == 'exact':
            for i, al in enumerate(alts, 1):
                g = car.get(i, 0)
                if g:
                    key = (chrom, pos, ref, al)
                    exact[key] = exact.get(key, 0) | g
            continue

        if a.mode == 'locus':
            count_locus(cnt, [(allele_class(ref, al), car[i]) for i, al in enumerate(alts, 1) if car.get(i)])
            continue

        # overlap: REF span [pos, pos+len(REF)-1]; uncarried records do not extend a locus
        if not car:
            continue
        start = int(pos)
        if chrom != cur_chrom or start > cur_end:
            flush()
            cur_chrom, cur_end = chrom, -1
        cur_end = max(cur_end, start + len(ref) - 1)
        for i, al in enumerate(alts, 1):
            g = car.get(i, 0)
            if g:
                key = (pos, ref, al)
                cur[key] = cur.get(key, 0) | g

    if a.mode == 'overlap':
        flush()

    for (chrom, pos, ref, al), g in exact.items():
        cnt[('allele_exact', allele_class(ref, al), CATEGORY[g])] += 1

    for (unit, cls, cat), n in sorted(cnt.items()):
        print(f'{a.varset}\t{unit}\t{cls}\t{cat}\t{n}')


if __name__ == '__main__':
    main()
