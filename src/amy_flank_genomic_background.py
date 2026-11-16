#!/usr/bin/env python3
"""Sliding-window genomic backgrounds and chromosome-stratified block-bootstrap
upper-tail genomic percentiles for AMY focal diversity and Hudson FST.

Main target: chr1 pilot (50-Mb physical blocks). Publication configuration:
GRCh38 autosomes 1-22, 5-Mb blocks, 10000 bootstrap replicates.
Requires an indexed genomewide joint VCF and bcftools on PATH.
"""
import argparse
import bisect
import csv
import gzip
import math
import random
import subprocess
import sys
from array import array
from collections import defaultdict
from pathlib import Path
import amy_flank_popgen as focal

GRCH38_LENGTHS = {
"1":248956422,"2":242193529,"3":198295559,"4":190214555,
"5":181538259,"6":170805979,"7":159345973,"8":145138636,
"9":138394717,"10":133797422,"11":135086622,"12":133275309,
"13":114364328,"14":107043718,"15":101991189,"16":90338345,
"17":83257441,"18":80373285,"19":58617616,"20":64444167,
"21":46709983,"22":50818468
}


class Coverage:
    def __init__(self, bed, chrom, size):
        raw=focal.mask_load(bed,chrom)
        if raw is None:
            raw=[(1,size)]
        merged=[]
        for a,b in sorted(raw):
            a=max(1,a); b=min(size,b)
            if b<a:
                continue
            if merged and a<=merged[-1][1]+1:
                merged[-1]=(merged[-1][0],max(merged[-1][1],b))
            else:
                merged.append((a,b))
        self.intervals=merged
        self.starts=[v[0] for v in merged]
        self.totals=[]
        n=0
        for a,b in merged:
            self.totals.append(n)
            n+=b-a+1

    def before(self,x):
        i=bisect.bisect_right(self.starts,x)-1
        if i<0:
            return 0
        a,b=self.intervals[i]
        return self.totals[i]+max(0,min(x,b)-a+1)

    def length(self,start,end):
        return self.before(end)-self.before(start-1)

    def contains(self,position):
        i=bisect.bisect_right(self.starts,position)-1
        return i>=0 and position<=self.intervals[i][1]


def load_lengths(path):
    d={}
    if path:
        with open(path) as f:
            for line in f:
                a=line.split()
                if len(a)>=2 and a[1].isdigit():
                    d[a[0]]=int(a[1])
    else:
        for k,v in GRCH38_LENGTHS.items():
            d["chr"+k]=v
            d[k]=v
    return d


def read_observed(root):
    points={}
    with open(Path(root)/"diversity.tsv") as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            for metric in ("pi_per_bp","thetaW_per_bp","Tajima_D"):
                z=focal.pf(float(row[metric])) if row[metric]!="NA" else "NA"
                if z!="NA":
                    points[(row["region"],row["population"],"",metric)]=float(z)
    fst_path=Path(root)/"hudson_fst.tsv"
    if fst_path.exists():
        with open(fst_path) as fh:
            for row in csv.DictReader(fh,delimiter="\t"):
                if row["Hudson_FST"]!="NA":
                    points[(row["region"],row["query"],row["reference"],"Hudson_FST")]=float(row["Hudson_FST"])
    return points


def series(chrom,populations,lists,args,mask):
    """Stream one chromosome once, retaining only per-site aggregate contributions."""
    metrics={}
    for pop in populations:
        for stat in ("pi","S"):
            metrics[(pop,"",stat)]=array("d",[0.])
    query=populations[0]
    for ref in populations[1:]:
        for stat in ("fn","fd"):
            metrics[(query,ref,stat)]=array("d",[0.])
    positions=array("I")
    command=["bcftools","view","-r",chrom,"-m2","-M2","-v","snps","-Ov",args.vcf]
    proc=subprocess.Popen(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                          bufsize=64*1024)
    groups=None
    seen=0; eligible=0
    for line in proc.stdout:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            names=line.rstrip("\n").split("\t")[9:]
            groups={}
            for p,samples in lists.items():
                if not set(samples).issubset(names):
                    raise ValueError("VCF missing sample(s) of "+p+" on "+chrom)
                groups[p]=[names.index(x) for x in samples]
            continue
        if not line.strip() or groups is None:
            continue
        r=line.rstrip("\n").split("\t")
        if len(r)<10 or r[3] not in "ACGT" or r[4] not in "ACGT":
            continue
        try:
            pos=int(r[1])
        except ValueError:
            continue
        if not mask.contains(pos):
            continue
        seen+=1
        genotypes=r[9:]
        vals={}
        for pop,ids in groups.items():
            ds=[focal.genotype(genotypes[i]) for i in ids]
            valid=[x for x in ds if x is not None]
            an=2*len(valid)
            min_an=math.ceil(2*len(ids)*args.callrate-1e-9)
            if an>=max(2,min_an):
                ac=sum(valid)
                vals[pop]=(ac,an)
            else:
                vals[pop]=None
        positions.append(pos)
        eligible+=1
        for pop in populations:
            q=vals[pop]
            pi=sp=0.
            if q is not None:
                ac,an=q
                if 0<ac<an:
                    pi=2*ac*(an-ac)/(an*(an-1))
                    projected=focal.projection(ac,an,args.projection)
                    if projected is not None:sp=projected
            A=metrics[(pop,"","pi")];A.append(A[-1]+pi)
            B=metrics[(pop,"","S")];B.append(B[-1]+sp)
        for ref in populations[1:]:
            q,rval=vals[query],vals[ref]
            fn=fd=0.
            if q is not None and rval is not None:
                ac,an=q;bc,bn=rval
                p1=ac/an;p2=bc/bn
                fd=p1*(1-p2)+p2*(1-p1)
                if fd>0 and an>1 and bn>1:
                    fn=(p1-p2)**2-p1*(1-p1)/(an-1)-p2*(1-p2)/(bn-1)
                else:
                    fn=fd=0.
            A=metrics[(query,ref,"fn")];A.append(A[-1]+fn)
            B=metrics[(query,ref,"fd")];B.append(B[-1]+fd)
        if eligible%100000==0:
            print("[GENOME] %s: %d callable SNPs processed"%(chrom,eligible),flush=True)
    rc=proc.wait()
    if rc:
        raise RuntimeError("bcftools failed on "+chrom+": "+proc.stderr.read()[:2000])
    if groups is None:
        raise ValueError("VCF header not found for "+chrom)
    if seen==0:
        print("[WARNING] no callable SNPs on "+chrom,file=sys.stderr)
    return positions,metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",required=True)
    p.add_argument("--pop",action="append",required=True)
    p.add_argument("--focal-stats",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--callable-bed",default="")
    p.add_argument("--distal",required=True)
    p.add_argument("--proximal",required=True)
    p.add_argument("--chromosomes",default="chr1")
    p.add_argument("--chrom-sizes",default="")
    p.add_argument("--block-mb",type=int,default=50)
    p.add_argument("--window-step",type=int,default=1000)
    p.add_argument("--bootstrap-reps",type=int,default=100)
    p.add_argument("--projection",type=int,default=168)
    p.add_argument("--callrate",type=float,default=.90)
    p.add_argument("--seed",type=int,default=20261117)
    p.add_argument("--fst-callable-match",type=float,default=.10,
                   help="FST background callable-bp match fraction; -1 disables")
    args=p.parse_args()
    if args.block_mb<1 or args.bootstrap_reps<2 or args.window_step<1:
        p.error("invalid block/window/replicate options")
    pops={}; chosen=[]
    for x in args.pop:
        name,path=x.split("=",1)
        if name in pops:
            p.error("duplicate population name")
        pops[name]=focal.read_lines(path);chosen.append(name)
    if len(set(x for group in pops.values() for x in group))!=sum(map(len,pops.values())):
        p.error("sample IDs must be unique and population lists disjoint")
    lengths=load_lengths(args.chrom_sizes)
    chroms=[x.strip() for x in args.chromosomes.split(",") if x.strip()]
    for chrom in chroms:
        if chrom not in lengths:
            p.error("chromosome length unknown: "+chrom+"; provide --chrom-sizes GRCh38.fai")
    regs={"distal":focal.parse_region(args.distal),"proximal":focal.parse_region(args.proximal)}
    if regs["distal"][0]!=regs["proximal"][0]:
        p.error("both flanks must be on same chromosome")
    span_start=min(x[1] for x in regs.values())
    span_end=max(x[2] for x in regs.values())
    amychrom=regs["distal"][0]
    points=read_observed(args.focal_stats)
    if not points:
        raise ValueError("no focal point estimates; run diversity and fst first")
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    block_bp=args.block_mb*1000000
    physical={}
    for chrom in chroms:
        physical[chrom]=math.ceil(lengths[chrom]/block_bp)
    total_phys=sum(physical.values())
    observations={key:{"windows":0,"tail":0,"blocks":defaultdict(lambda:[0,0])} for key in points}
    n0=args.projection
    a1=sum(1/i for i in range(1,n0))
    a2=sum(1/(i*i) for i in range(1,n0))
    b1=(n0+1)/(3*(n0-1));b2=2*(n0*n0+n0+3)/(9*n0*(n0-1))
    e1=(b1-1/a1)/a1
    e2=(b2-(n0+2)/(a1*n0)+a2/(a1*a1))/(a1*a1+a2)
    focal_callable={}
    for reg,(_,start,end) in regs.items():
        cover=Coverage(args.callable_bed,amychrom,lengths[amychrom])
        focal_callable[reg]=cover.length(start,end)
    for chrom in chroms:
        cover=Coverage(args.callable_bed,chrom,lengths[chrom])
        positions,prefix=series(chrom,chosen,pops,args,cover)
        size=lengths[chrom]
        for reg,(_,start,end) in regs.items():
            winlen=end-start+1
            key_subset=[k for k in points if k[0]==reg]
            for left in range(1,size-winlen+2,args.window_step):
                right=left+winlen-1
                if chrom==amychrom and left<=span_end and right>=span_start:
                    continue
                length=cover.length(left,right)
                if length==0:
                    continue
                lo=bisect.bisect_left(positions,left)
                hi=bisect.bisect_right(positions,right)
                block=(chrom,(left+winlen//2-1)//block_bp)
                by_value={}
                for name in chosen:
                    s=prefix[(name,"","S")][hi]-prefix[(name,"","S")][lo]
                    pi_sum=prefix[(name,"","pi")][hi]-prefix[(name,"","pi")][lo]
                    by_value[(name,"","pi_per_bp")]=pi_sum/length
                    by_value[(name,"","thetaW_per_bp")]=s/(a1*length)
                    var=e1*s+e2*s*(s-1)
                    by_value[(name,"","Tajima_D")]=(pi_sum-s/a1)/math.sqrt(var) if var>0 else float("nan")
                for other in chosen[1:]:
                    n=prefix[(chosen[0],other,"fn")][hi]-prefix[(chosen[0],other,"fn")][lo]
                    d=prefix[(chosen[0],other,"fd")][hi]-prefix[(chosen[0],other,"fd")][lo]
                    by_value[(chosen[0],other,"Hudson_FST")]=n/d if d>0 else float("nan")
                for key in key_subset:
                    _,p1,p2,metric=key
                    value=by_value[(p1,p2,metric)]
                    if not math.isfinite(value):
                        continue
                    if metric=="Hudson_FST" and args.fst_callable_match>=0:
                        focalL=focal_callable[reg]
                        if focalL==0 or abs(length-focalL)>args.fst_callable_match*focalL:
                            continue
                    dest=observations[key]
                    tail=int(value>=points[key])
                    dest["windows"]+=1;dest["tail"]+=tail
                    b=dest["blocks"][block]
                    b[0]+=1;b[1]+=tail
        print("[GENOME] finished "+chrom,flush=True)
    rng=random.Random(args.seed)
    rows=[]
    for key,d in observations.items():
        reg,p1,p2,metric=key
        n=d["windows"]
        point=100*d["tail"]/n if n else float("nan")
        boot=[]
        if n:
            for iteration in range(args.bootstrap_reps):
                tot=hit=0
                for chrom in chroms:
                    for _ in range(physical[chrom]):
                        index=rng.randrange(physical[chrom])
                        count,extreme=d["blocks"].get((chrom,index),(0,0))
                        tot+=count;hit+=extreme
                if tot:
                    boot.append(100*hit/tot)
        qlow=percentile(boot,.025) if boot else float("nan")
        qhigh=percentile(boot,.975) if boot else float("nan")
        rows.append((reg,p1,p2,metric,focal.pf(points[key]),focal.pf(point),n,
                     len(d["blocks"]),total_phys,args.block_mb,args.window_step,
                     args.bootstrap_reps,len(boot),
                     focal.pf(sum(boot)/len(boot)) if boot else "NA",
                     focal.pf(statistics.stdev(boot)) if len(boot)>1 else "NA",
                     focal.pf(qlow),focal.pf(qhigh),args.chromosomes,
                     args.fst_callable_match if metric=="Hudson_FST" else "unmatched"))
    focal.write(out/"genomic_percentile_ci.tsv",
       ["region","population","reference","statistic","focal_value",
        "upper_tail_genomic_percentile","background_windows",
        "occupied_physical_blocks","all_physical_blocks",
        "block_size_mb","window_step_bp","bootstrap_reps","valid_reps",
        "bootstrap_mean","bootstrap_sd","CI_025","CI_975","chromosomes",
        "fst_callable_match_fraction"],rows)
    print("[OK] "+str(out/"genomic_percentile_ci.tsv"))


def percentile(values,q):
    data=sorted(x for x in values if math.isfinite(x))
    if not data:return float("nan")
    z=(len(data)-1)*q
    a=int(z);b=math.ceil(z)
    return data[a]*(b-z)+data[b]*(z-a) if a!=b else data[a]


if __name__=="__main__":
    import statistics
    main()
