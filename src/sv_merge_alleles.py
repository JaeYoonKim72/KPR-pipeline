#!/usr/bin/env python3
"""Cluster SV alleles with PanSVMerger and VSEARCH; remap genotypes and allele counts."""
import argparse
import importlib.util
import multiprocessing as mp
import os
import random
import string
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

import pysam


def load_pansvmerger(path):
    spec = importlib.util.spec_from_file_location('pansvmerger', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- worker side ----
W = {}


def worker_init(psm_path, tmpdir, identity, maxseqlength, minseqlength, vsearch, threads):
    W.update(psm=load_pansvmerger(psm_path), tmpdir=tmpdir, identity=identity,
             maxseqlength=maxseqlength, minseqlength=minseqlength, vsearch=vsearch, threads=threads)


def cluster_job(job):
    """job: (chrom, pos, ref, alts) -> PanSVMerger allele->cluster mapping (or None)"""
    chrom, pos, ref, alts = job
    psm = W['psm']
    tag = ''.join(random.choices(string.ascii_letters + string.digits, k=12))
    fa = os.path.join(W['tmpdir'], f'{tag}.fa')
    try:
        psm.write_tmp_fasta(chrom, pos, ref, alts, fa)
        # Same command as PanSVMerger run_vsearch(); output captured instead of printed.
        cmd = [W['vsearch'], '--cluster_fast', fa, '--id', str(W['identity']), '--strand', 'both',
               '--uc', f'{fa}.uc', '--threads', str(W['threads']),
               '--maxseqlength', str(W['maxseqlength']), '--minseqlength', str(W['minseqlength'])]
        r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f'vsearch failed at {chrom}:{pos}: {r.stderr[-500:]}')
        clusters, _ = psm.parse_vsearch_clusters(f'{fa}.uc')
        if (chrom, pos) not in clusters:
            return None
        clusters = psm.ensure_ref_and_alts_in_clusters(SimpleNamespace(chrom=chrom, pos=pos, alts=alts), clusters)
        return clusters[(chrom, pos)]
    except Exception as e:
        raise RuntimeError(f'{chrom}:{pos} (n_alt={len(alts)}): {type(e).__name__}: {e}') from e
    finally:
        for f in (fa, f'{fa}.uc'):
            if os.path.exists(f):
                os.remove(f)


# ---- main side ----
DROP = object()  # merged_record() result: record removed


def merged_record(rec, mapping, header, keep_info):
    ref_cluster = mapping[0]
    cluster_to_gt = {ref_cluster: 0}
    for c in sorted(set(mapping.values())):
        if c not in cluster_to_gt:
            cluster_to_gt[c] = len(cluster_to_gt)
    n_alt = len(cluster_to_gt) - 1
    if n_alt == 0:
        return DROP  # every ALT clustered with REF -> no SV left at this locus
    if n_alt == len(rec.alts):
        return None  # nothing merged

    # Longest member sequence represents each cluster
    alts = [None] * n_alt
    for allele_idx, c in mapping.items():
        g = cluster_to_gt[c]
        if g == 0 or allele_idx > len(rec.alts):
            continue
        seq = rec.alts[allele_idx - 1]
        if alts[g - 1] is None or len(seq) > len(alts[g - 1]):
            alts[g - 1] = seq
    alts = [a if a is not None else rec.ref for a in alts]

    old_to_new = {i: cluster_to_gt[mapping[i]] for i in range(len(rec.alts) + 1)}
    info = {}
    for k, v in rec.info.items():
        if k not in keep_info:
            continue
        # Normalize merged scalar INFO and flag fields.
        if header.info[k].type == 'Flag':
            v = True
        elif header.info[k].number == 1:
            if isinstance(v, tuple):
                if header.info[k].type != 'String':
                    continue
                v = '|'.join(map(str, v))
            elif isinstance(v, str):
                v = v.replace(',', '|')
        info[k] = v
    new = header.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop,
                            alleles=[rec.ref] + alts, id=rec.id, qual=rec.qual,
                            filter=list(rec.filter), info=info)
    ac = [0] * n_alt
    an = 0
    for s in rec.samples:
        gt = rec.samples[s]['GT']
        new_gt = tuple(None if a is None else old_to_new[a] for a in gt)
        new.samples[s]['GT'] = new_gt
        new.samples[s].phased = rec.samples[s].phased
        for a in new_gt:
            if a is not None:
                an += 1
                if a > 0:
                    ac[a - 1] += 1
    new.info['AC'] = tuple(ac)
    new.info['AN'] = an
    if 'AF' in header.info:
        new.info['AF'] = tuple(x / an if an else 0.0 for x in ac)
    return new


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--in', dest='inp', required=True, help='SV multi-allelic VCF (.vcf.gz, indexed)')
    p.add_argument('--out', required=True, help='output .vcf.gz')
    p.add_argument('--pansvmerger', required=True, help='path to SV.merge_variants.py')
    p.add_argument('--identity', type=float, default=0.9)
    p.add_argument('--alt-length', type=int, default=50000, help='same as PanSVMerger -al')
    p.add_argument('--alt-number', type=int, default=100, help='same as PanSVMerger -an')
    p.add_argument('--maxseqlength', type=int, default=100000)
    p.add_argument('--minseqlength', type=int, default=2)
    p.add_argument('--workers', type=int, default=1, help='parallel vsearch processes')
    p.add_argument('--vsearch-threads', type=int, default=1, help='threads per vsearch')
    p.add_argument('--batch', type=int, default=2000, help='multi-ALT records per batch')
    p.add_argument('--tmpdir', default=None)
    a = p.parse_args()

    vsearch = __import__('shutil').which('vsearch')
    if not vsearch:
        sys.exit('ERROR: vsearch not in PATH')
    vin = pysam.VariantFile(a.inp)
    header = vin.header
    # Per-allele INFO fields (Number=A/R/G) cannot be carried over after merging
    keep_info = {k for k, v in header.info.items() if v.number not in ('A', 'R', 'G')}
    vout = pysam.VariantFile(a.out, 'wz', header=header)
    tmpdir = tempfile.mkdtemp(dir=a.tmpdir or os.path.dirname(os.path.abspath(a.out)))

    pool = mp.Pool(a.workers, initializer=worker_init,
                   initargs=(a.pansvmerger, tmpdir, a.identity, a.maxseqlength, a.minseqlength,
                             vsearch, a.vsearch_threads))
    n = {'in': 0, 'merged': 0, 'dropped': 0, 'clustered': 0}
    t0 = time.time()

    def flush(buf):
        jobs = [(r.chrom, r.pos, r.ref, tuple(r.alts)) for r, need in buf if need]
        maps = iter(pool.imap(cluster_job, jobs, chunksize=1))
        for rec, need in buf:
            if not need:
                vout.write(rec)
                continue
            mapping = next(maps)
            new = merged_record(rec, mapping, header, keep_info) if mapping else None
            if new is DROP:
                n['dropped'] += 1
            elif new is None:
                vout.write(rec)
            else:
                vout.write(new)
                n['merged'] += 1
        n['clustered'] += len(jobs)
        print(f'{time.strftime("%T")} records={n["in"]} clustered={n["clustered"]} '
              f'({time.time() - t0:.0f}s)', file=sys.stderr, flush=True)

    buf, n_jobs = [], 0
    for rec in vin.fetch():
        n['in'] += 1
        alts = rec.alts or ()
        need = len(alts) > 1 and sum(len(x) >= a.alt_length for x in alts) <= a.alt_number
        buf.append((rec, need))
        n_jobs += need
        if n_jobs >= a.batch:
            flush(buf)
            buf, n_jobs = [], 0
    flush(buf)
    pool.close()
    pool.join()
    vout.close()
    os.rmdir(tmpdir)
    pysam.tabix_index(a.out, preset='vcf', force=True)
    print(f'records={n["in"]}\tclustered={n["clustered"]}\tmerged_records={n["merged"]}'
          f'\tdropped_all_ref={n["dropped"]}')


if __name__ == '__main__':
    main()
