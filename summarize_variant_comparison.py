#!/usr/bin/env python3
"""Sum per-chromosome counts (count_groups.py output) into the summary tables.

  variant_comparison_summary.tsv : SV / small variant x locus / allele,
                                   target-specific, shared, other-specific, novel %
  variant_comparison_by_class.tsv: the same units split by variant class (SNV / INS / DEL / COMPLEX)
  SV new alleles at shared loci  : printed and written to variant_comparison_shared_loci.tsv

Units used (see compare_pangenome_variant_counts.sh, step count):
  SV    locus  -> unit 'locus'            (vcfcreatemulti locus, mode locus)
  SV    allele -> unit 'allele_in_record' (PanSVMerger allele cluster)
  small locus  -> unit 'locus'            (overlapping REF spans, mode overlap)
  small allele -> unit 'allele_exact'     (exact CHROM/POS/REF/ALT, mode exact)
Novel % = target-specific / (shared + other-specific) x 100.
Locus class OTHER (equal-length multi-base substitution) is reported under COMPLEX.
"""
import argparse
import glob
import os
from collections import defaultdict

UNITS = [('SV (>50 bp)', 'sv', 'locus', 'locus'),
         ('SV (>50 bp)', 'sv', 'allele', 'allele_in_record'),
         ('Small variant (<=50 bp)', 'small', 'locus', 'locus'),
         ('Small variant (<=50 bp)', 'small', 'allele', 'allele_exact')]
CATS = ['target_only', 'shared', 'other_only']
CLASS_ORDER = ['ALL', 'SNV', 'INS', 'DEL', 'COMPLEX']


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--counts', required=True, help='directory with <chr>.tsv')
    p.add_argument('--out', required=True, help='output directory')
    p.add_argument('--target-label', default='KPR')
    p.add_argument('--other-label', default='HPRC')
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.counts, '*.tsv')))
    if not files:
        raise SystemExit(f'no <chr>.tsv in {a.counts}')
    tab = defaultdict(lambda: defaultdict(int))      # (varset, unit, class) -> category -> n
    for path in files:
        with open(path) as f:
            for line in f:
                _, varset, unit, cls, cat, n = line.rstrip('\n').split('\t')
                if cls == 'OTHER' and unit == 'locus':
                    cls = 'COMPLEX'
                tab[(varset, unit, cls)][cat] += int(n)
                tab[(varset, unit, 'ALL')][cat] += int(n)

    t, o = a.target_label, a.other_label
    head = ['variant', 'unit', 'class', f'{t}_specific', 'shared', f'{o}_specific', 'novel_pct']

    def row(label, varset, unit_label, unit, cls):
        c = tab.get((varset, unit, cls))
        if not c:
            return None
        n = [c.get(x, 0) for x in CATS]
        pct = 100 * n[0] / (n[1] + n[2]) if n[1] + n[2] else 0
        return [label, unit_label, cls] + n + [f'{pct:.2f}']

    summary, by_class = [], []
    for label, varset, unit_label, unit in UNITS:
        summary.append(row(label, varset, unit_label, unit, 'ALL'))
        for cls in CLASS_ORDER:
            r = row(label, varset, unit_label, unit, cls)
            if r:
                by_class.append(r)

    def write(name, rows, cols):
        with open(os.path.join(a.out, name), 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for r in rows:
                if r:
                    f.write('\t'.join(map(str, r)) + '\n')

    write('variant_comparison_summary.tsv', summary, head)
    write('variant_comparison_by_class.tsv', by_class, head)

    # SV loci shared by both groups that carry a target-specific allele (new alleles at known SV sites)
    n_loci = tab.get(('sv', 'shared_locus_with_target_private_allele', 'ALL'), {}).get('shared', 0)
    n_alleles = tab.get(('sv', 'allele_in_shared_locus', 'ALL'), {}).get('target_only', 0)
    write('variant_comparison_shared_loci.tsv',
          [[f'shared SV loci with a {t}-specific allele', n_loci],
           [f'{t}-specific SV alleles in shared loci', n_alleles]], ['item', 'count'])

    print('\t'.join(head))
    for r in summary:
        if r:
            print('\t'.join(map(str, r)))
    print(f'{n_alleles} {t}-specific SV alleles in {n_loci} shared SV loci')


if __name__ == '__main__':
    main()
