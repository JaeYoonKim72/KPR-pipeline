#!/usr/bin/env python3
"""Keep autosomal (chr1-22) SVs of one type (INS or DEL) and >= --min-size bp.

Works for sequence-resolved and symbolic (<DEL>, <INS>) records. INFO/SVTYPE and a
signed INFO/SVLEN are written for every kept record. Input must be biallelic
(bcftools norm -m -any first). Used by pangenome_variant_benchmarking.sh, step sv_truth_readsupport.
"""
import argparse
import sys

import pysam

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--input', required=True)
p.add_argument('--output', required=True)
p.add_argument('--type', required=True, choices=['INS', 'DEL'])
p.add_argument('--min-size', type=int, default=51)
a = p.parse_args()

chroms = {f'chr{i}' for i in range(1, 23)}
counts = dict(total=0, autosome=0, type=0, size=0, kept=0, unresolved=0)

with pysam.VariantFile(a.input) as inv:
    hdr = inv.header.copy()
    if 'SVTYPE' not in hdr.info:
        hdr.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">')
    if 'SVLEN' not in hdr.info:
        hdr.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Signed structural variant length">')
    with pysam.VariantFile(a.output, 'w', header=hdr) as out:
        for rec in inv:
            counts['total'] += 1
            if rec.chrom not in chroms:
                continue
            counts['autosome'] += 1
            ref = rec.ref or ''
            for alt in rec.alts or ():
                if not alt or alt == '*':
                    continue
                typ = length = None
                if alt.startswith('<') and alt.endswith('>'):
                    typ = alt[1:-1].split(':')[0]
                    val = rec.info.get('SVLEN')
                    if isinstance(val, tuple):
                        val = val[0] if val else None
                    if val is not None:
                        length = abs(int(val))
                    elif rec.stop is not None and rec.start is not None and typ == 'DEL':
                        length = abs(rec.stop - rec.start)
                elif all(c in 'ACGTNacgtn' for c in alt + ref):
                    diff = len(alt) - len(ref)
                    typ = 'INS' if diff > 0 else 'DEL' if diff < 0 else None
                    length = abs(diff)
                if typ != a.type:
                    continue
                counts['type'] += 1
                if length is None:
                    counts['unresolved'] += 1
                    continue
                if length < a.min_size:
                    continue
                counts['size'] += 1
                if len(rec.alts or ()) != 1:
                    raise ValueError(f'Multi-allelic record at {rec.chrom}:{rec.pos}; run bcftools norm -m -any first')
                rec.translate(hdr)
                rec.info['SVTYPE'] = typ
                rec.info['SVLEN'] = length if typ == 'INS' else -length
                out.write(rec)
                counts['kept'] += 1
                break

print('\t'.join(f'{k}={v}' for k, v in counts.items()), file=sys.stderr)
