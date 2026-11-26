#!/usr/bin/env python3
"""Individual-level bootstrap and delete-one jackknife for focal AMY SNP statistics.

Uses the same allele-count, Hudson FST and dosage-r2 estimators as
amy_flank_popgen.py. No paralog-level genotypes are used.
"""
import argparse
import csv
import math
import random
import statistics
from pathlib import Path
import amy_flank_popgen as focal


def numeric(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (ValueError, TypeError):
        return float("nan")


def percentile(vals, q):
    values = sorted(v for v in vals if math.isfinite(v))
    if not values:
        return float("nan")
    i = (len(values)-1) * q
    a = int(math.floor(i))
    b = int(math.ceil(i))
    return values[a]*(b-i)+values[b]*(i-a) if a!=b else values[a]


def update_rows(rows, indices):
    if not indices:
        return rows
    result = []
    for pos, data in rows:
        new = data.copy()
        for pop, picked in indices.items():
            gs = data[pop][2]
            vals = [gs[i] for i in picked]
            observed = [x for x in vals if x is not None]
            new[pop] = (sum(observed), 2*len(observed), vals)
        result.append((pos, new))
    return result


def value(stat, region, p, ref, data, lengths, args, fst_rate=None):
    if stat in ("pi_per_bp", "thetaW_per_bp", "Tajima_D"):
        result = focal.diversity(data[region], p, lengths[region], args.projection, args.callrate)
        return numeric(result[{"pi_per_bp":2,"thetaW_per_bp":3,"Tajima_D":4}[stat]])
    if stat == "mean_corrected_r2":
        maxlen = args.distal_ld_max if region == "distal" else args.proximal_ld_max
        return numeric(focal.ld(data[region], p, None, maxlen, args.callrate, args.ld_maf)[1])
    if stat == "across_corrected_r2":
        return numeric(focal.cross_ld(data["distal"], data["proximal"], p, args.callrate, args.ld_maf)[1])
    if stat == "Hudson_FST":
        return numeric(focal.fst(data[region], p, ref, args.callrate if fst_rate is None else fst_rate)[1])
    raise ValueError("unknown statistic: "+stat)


def sample_variance(values):
    v = [x for x in values if math.isfinite(x)]
    return statistics.variance(v) if len(v) >= 2 else float("nan")


def run(args):
    chosen = {}
    for item in args.pop:
        p, path = item.split("=",1)
        if p in chosen or not p:
            raise ValueError("population repeated/empty: "+p)
        chosen[p] = focal.read_lines(path)
    if len(set(x for s in chosen.values() for x in s)) != sum(map(len,chosen.values())):
        raise ValueError("population lists must be disjoint and contain no duplicate IDs")
    regs = {"distal":focal.parse_region(args.distal),
            "proximal":focal.parse_region(args.proximal)}
    mask = focal.mask_load(args.callable_bed, regs["distal"][0])
    data = focal.read_vcf(args.vcf, chosen, regs, args.callrate, mask)
    lengths = {k:focal.callable_len(x[1],x[2],mask) for k,x in regs.items()}
    if any(v<=0 for v in lengths.values()):
        raise ValueError("zero callable length in flank")
    pops = list(chosen)
    definitions = []
    for region in ("distal","proximal"):
        for p in pops:
            definitions.extend((metric,region,p,"") for metric in ("pi_per_bp","thetaW_per_bp","Tajima_D"))
            definitions.append(("mean_corrected_r2",region,p,""))
        if len(pops)>1:
            definitions.extend(("Hudson_FST",region,pops[0],ref) for ref in pops[1:])
    definitions.extend(("across_corrected_r2","across_flanks",p,"") for p in pops)
    out = Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    rnd = random.Random(args.seed)
    summaries, replicates = [], []
    for item_idx, (metric, region, pop, ref) in enumerate(definitions):
        metric_data = data
        fst_rate = None
        if metric == "Hudson_FST":
            metric_data = {region: focal.fst_panel(data[region], pop, ref, args.callrate)}
            fst_rate = 0.0
        point = value(metric,region,pop,ref,metric_data,lengths,args,fst_rate)
        involved = [pop]+([ref] if ref else [])
        boot=[]
        for i in range(args.bootstrap_reps):
            picks = {p:[rnd.randrange(len(chosen[p])) for _ in chosen[p]] for p in involved}
            affected = {r:update_rows(rows,picks) for r,rows in metric_data.items()
                        if r==region or region=="across_flanks"}
            if region!="across_flanks":
                affected={region:affected[region]}
            score=value(metric,region,pop,ref,affected,lengths,args,fst_rate)
            replicates.append((metric,region,pop,ref,"bootstrap",i+1,focal.pf(score)))
            if math.isfinite(score):
                boot.append(score)
        bootmean=statistics.fmean(boot) if boot else float("nan")
        bootse=statistics.stdev(boot) if len(boot)>=2 else float("nan")
        low,high=percentile(boot,.025),percentile(boot,.975)
        bias=bootmean-point
        jack_var=0.
        jack_valid=True
        jack_counts=[]
        for p in involved:
            n=len(chosen[p])
            requested=max(2, math.ceil(args.jackknife_fraction*n))
            used=min(n,requested)
            chosen_i=sorted(rnd.sample(range(n),used))
            jack_counts.append("%s:%d/%d"%(p,used,n))
            samples=[]
            for j in chosen_i:
                picks={p:[k for k in range(n) if k!=j]}
                affected={r:update_rows(rows,picks) for r,rows in metric_data.items()
                          if r==region or region=="across_flanks"}
                if region!="across_flanks":
                    affected={region:affected[region]}
                score=value(metric,region,pop,ref,affected,lengths,args,fst_rate)
                replicates.append((metric,region,pop,ref,"jackknife",j+1,focal.pf(score)))
                if math.isfinite(score):
                    samples.append(score)
            variance=sample_variance(samples)
            if not math.isfinite(variance):
                jack_valid=False
            else:
                # Complete delete-one JK: (n-1)/n * sum((theta_-i - mean)^2).
                # For m<n randomly selected delete-one estimates, their
                # sample variance estimates the finite-population variance.
                # Hence this is explicitly an APPROXIMATE subset jackknife.
                jack_var+=((n-1)**2/n)*variance
        jackse=math.sqrt(jack_var) if jack_valid else float("nan")
        status="FULL_DELETE_ONE" if args.jackknife_fraction>=1 else "SUBSET_JACKKNIFE_APPROX"
        if not math.isfinite(point):
            status="INVALID_POINT"
        summaries.append((metric,region,pop,ref,focal.pf(point),len(boot),
                          focal.pf(bootmean),focal.pf(bootse),focal.pf(bias),
                          focal.pf(low),focal.pf(high),
                          focal.pf(low-bias),focal.pf(high-bias),
                          ",".join(jack_counts),focal.pf(jackse),
                          focal.pf(point-1.96*jackse),
                          focal.pf(point+1.96*jackse),status))
        print("[CI] %s %s %s %s: bootstrap %d/%d, JK %s" %
              (metric,region,pop,ref,len(boot),args.bootstrap_reps,status),flush=True)
    focal.write(out/"individual_ci.tsv",
       ["statistic","region","population","reference","point_estimate","bootstrap_n_valid",
        "bootstrap_mean","bootstrap_se","bootstrap_bias",
        "bootstrap_raw_025","bootstrap_raw_975",
        "bootstrap_shifted_025","bootstrap_shifted_975",
        "jackknife_used_per_group","jackknife_se",
        "jackknife_CI_025","jackknife_CI_975","jackknife_status"],
       summaries)
    focal.write(out/"individual_replicates.tsv",
       ["statistic","region","population","reference","method","replicate_or_deleted_index","value"],
       replicates)
    print("[OK] output:",out/"individual_ci.tsv")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf",required=True)
    ap.add_argument("--pop",action="append",required=True)
    ap.add_argument("--distal",required=True)
    ap.add_argument("--proximal",required=True)
    ap.add_argument("--callable-bed",default="")
    ap.add_argument("--out",required=True)
    ap.add_argument("--projection",type=int,default=168)
    ap.add_argument("--callrate",type=float,default=.90)
    ap.add_argument("--bootstrap-reps",type=int,default=100)
    ap.add_argument("--jackknife-fraction",type=float,default=.30)
    ap.add_argument("--seed",type=int,default=20261117)
    ap.add_argument("--ld-maf",type=float,default=.05)
    ap.add_argument("--distal-ld-max",type=int,default=100000)
    ap.add_argument("--proximal-ld-max",type=int,default=60000)
    args=ap.parse_args()
    if args.bootstrap_reps<2:
        ap.error("bootstrap-reps must be >=2")
    if not 0<args.jackknife_fraction<=1:
        ap.error("jackknife-fraction must be (0,1]")
    if not 0<=args.ld_maf<.5:
        ap.error("ld-maf must be between 0 and 0.5")
    run(args)

if __name__=="__main__":
    main()
