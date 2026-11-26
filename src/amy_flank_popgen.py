#!/usr/bin/env python3
"""Diversity, LD and Hudson FST in AMY flanking regions."""
import argparse
import csv
import gzip
import math
from pathlib import Path

DEFAULT_DISTAL = (103456163, 103561526)
DEFAULT_PROXIMAL = (103760698, 103826698)

def read_lines(path):
    with open(path) as f:
        return [s.strip().split()[0] for s in f if s.strip() and not s.startswith("#")]

def vcf_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def parse_region(s):
    x = s.replace(",", "")
    ch, rest = x.split(":")
    a, b = map(int, rest.split("-"))
    if a < 1 or b < a: raise ValueError("invalid interval: "+s)
    return ch, a, b

def mask_load(bed, ch):
    if not bed: return None
    segments = []
    opener = gzip.open if bed.endswith(".gz") else open
    with opener(bed, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            s=line.split()
            if s[0] in (ch, ch.replace("chr",""), "chr"+ch):
                segments.append((int(s[1])+1, int(s[2])))  # BED 0-based -> 1-based
    return sorted(segments)

def in_mask(pos, mask):
    if mask is None: return True
    return any(a <= pos <= b for a,b in mask)

def callable_len(a, b, mask):
    if mask is None: return b-a+1
    pieces=sorted((max(a,x),min(b,y)) for x,y in mask if min(b,y)>=max(a,x))
    n=0
    prev=-1
    for x,y in pieces:
        x=max(x,prev+1)
        if y>=x: n+=y-x+1; prev=y
    return n

def genotype(call):
    gt=call.split(":",1)[0]
    if "|" in gt: sep="|"
    elif "/" in gt: sep="/"
    else: return None
    x=gt.split(sep)
    if len(x)!=2 or any(v not in ("0","1") for v in x): return None
    return int(x[0])+int(x[1])

def read_vcf(path, selected, regions, callrate, mask):
    result={k:[] for k in regions}
    with vcf_open(path) as f:
        cols=None
        for line in f:
            if line.startswith("##"):continue
            if line.startswith("#CHROM"):
                names=line.rstrip().split("\t")[9:]
                missing={p:set(s)-set(names) for p,s in selected.items()}
                if any(missing.values()):
                    raise ValueError("samples absent from VCF: "+str({p:sorted(v)[:5] for p,v in missing.items() if v}))
                cols={p:[names.index(s) for s in samples] for p,samples in selected.items()}
                if any(not i for i in cols.values()):raise ValueError("empty population")
                continue
            if cols is None:raise ValueError("VCF missing #CHROM header")
            a=line.rstrip().split("\t")
            if len(a)<10 or a[0] not in (regions["distal"][0], regions["distal"][0].replace("chr",""),"chr"+regions["distal"][0]):continue
            try:pos=int(a[1])
            except ValueError:continue
            where=[k for k,(_,x,y) in regions.items() if x<=pos<=y]
            if not where or not in_mask(pos,mask):continue
            if len(a[3])!=1 or len(a[4])!=1 or a[3] not in "ACGT" or a[4] not in "ACGT" or a[3]==a[4]:continue
            cells=a[9:]
            counts={}
            for pop,ids in cols.items():
                ds=[genotype(cells[i]) for i in ids]
                ok=[v for v in ds if v is not None]
                counts[pop]=(sum(ok),2*len(ok),ds)
            # Apply call-rate filters separately to each population comparison.
            record=(pos,counts)
            for k in where:result[k].append(record)
    return result

def logchoose(n,k):
    if k<0 or k>n:return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def projection(a,n,n0):
    if n<n0 or not (0<a<n):return None
    # 1 - P(all reference) - P(all alternative)
    return max(0.,min(1., 1.-math.exp(logchoose(n-a,n0)-logchoose(n,n0))-
                             math.exp(logchoose(a,n0)-logchoose(n,n0))))

def tajima_d(pi_sum,s,n):
    if n<3 or s<=0:return float("nan")
    a1=sum(1/i for i in range(1,n));a2=sum(1/i**2 for i in range(1,n))
    b1=(n+1)/(3*(n-1));b2=2*(n*n+n+3)/(9*n*(n-1))
    c1=b1-1/a1;c2=b2-(n+2)/(a1*n)+a2/(a1*a1)
    e1=c1/a1;e2=c2/(a1*a1+a2)
    v=e1*s+e2*s*(s-1)
    return (pi_sum-s/a1)/math.sqrt(v) if v>0 else float("nan")

def pf(x):
    return "NA" if not isinstance(x,(float,int)) or not math.isfinite(x) else "%.9g"%x

def diversity(rows,pop,l,n0,rate):
    pi=0.;sp=0.;used=0;raw=0
    for _,d in rows:
        ac,an,_=d[pop]
        if an < 2 or an < math.ceil(2*len(d[pop][2])*rate-1e-9):continue
        if 0<ac<an:
            raw+=1
            pi+=2*ac*(an-ac)/(an*(an-1))
            q=projection(ac,an,n0)
            if q is not None:sp+=q
            used+=1
    a1=sum(1/i for i in range(1,n0))
    theta=sp/a1/l if l and a1 else float("nan")
    return (used,raw,pf(pi/l if l else float("nan")),pf(theta),pf(tajima_d(pi,sp,n0)),pf(sp))

def r2_corrected(x,y):
    pairs=[(a,b) for a,b in zip(x,y) if a is not None and b is not None]
    n=len(pairs)
    if n<=2:return None
    mx=sum(a for a,b in pairs)/n;my=sum(b for a,b in pairs)/n
    xx=sum((a-mx)**2 for a,b in pairs); yy=sum((b-my)**2 for a,b in pairs)
    if xx==0 or yy==0:return None
    xy=sum((a-mx)*(b-my) for a,b in pairs)
    return ((n-1)*(xy*xy/(xx*yy))-1)/(n-2)

def ld(rows,pop,dmin,dmax,rate,maf=.05):
    use=[]
    for pos,d in rows:
        ac,an,g=d[pop]
        if an<2*len(g)*rate-1e-9 or not (0<ac<an):continue
        site_maf=min(ac/an,1-ac/an)
        if site_maf>=maf:use.append((pos,g))
    vals=[]
    for i,(p,x) in enumerate(use):
        for q,y in use[i+1:]:
            dist=abs(q-p)
            if dmin is not None and dist<dmin:continue
            if dmax is not None and dist>=dmax:continue
            v=r2_corrected(x,y)
            if v is not None:vals.append(v)
    return len(vals), pf(sum(vals)/len(vals) if vals else float("nan"))

def cross_ld(a,b,pop,rate,maf=.05):
    left=[(p,d[pop]) for p,d in a]; right=[(p,d[pop]) for p,d in b]
    def good(c):
        ac,an,g=c
        return an>=2*len(g)*rate-1e-9 and an>0 and 0<ac<an and min(ac/an,1-ac/an)>=maf
    x=[c[2] for p,c in left if good(c)]
    y=[c[2] for p,c in right if good(c)]
    s=0.;n=0
    for i in x:
        for j in y:
            v=r2_corrected(i,j)
            if v is not None:s+=v;n+=1
    return n,pf(s/n if n else float("nan"))

def fst_panel(rows, p, q, rate):
    return [(pos, d) for pos, d in rows
            if all(d[g][1] >= max(2, 2*len(d[g][2])*rate-1e-9) for g in (p, q))]

def fst(rows,p,q,rate):
    num=den=0.;used=0
    for _,d in rows:
        a,n,_=d[p];b,m,_=d[q]
        if n<2 or m<2:continue
        if n<2*len(d[p][2])*rate-1e-9 or m<2*len(d[q][2])*rate-1e-9:continue
        pa=a/n;pb=b/m
        dd=pa*(1-pb)+pb*(1-pa)
        if dd<=0:continue
        nn=(pa-pb)**2-pa*(1-pa)/(n-1)-pb*(1-pb)/(m-1)
        num+=nn;den+=dd;used+=1
    return used,pf(num/den if den else float("nan"))

def write(path,header,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(header);w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--vcf",required=True)
    ap.add_argument("--pop",action="append",required=True,help="NAME=sample-list.txt; first population is query")
    ap.add_argument("--distal",default="chr1:103456163-103561526")
    ap.add_argument("--proximal",default="chr1:103760698-103826698")
    ap.add_argument("--callable-bed",default="")
    ap.add_argument("--out",required=True)
    ap.add_argument("--projection",type=int,default=168)
    ap.add_argument("--callrate",type=float,default=.9)
    ap.add_argument("--ld-maf",type=float,default=.05)
    ap.add_argument("--distal-ld-max",type=int,default=100000)
    ap.add_argument("--proximal-ld-max",type=int,default=60000)
    ap.add_argument("--only",choices=("all","diversity","ld","across_ld","fst"),default="all")
    args=ap.parse_args()
    regs={"distal":parse_region(args.distal),"proximal":parse_region(args.proximal)}
    if regs["distal"][0]!=regs["proximal"][0]:ap.error("regions must use same chromosome")
    if args.projection<3:ap.error("projection must be >=3")
    if not (0<args.callrate<=1):ap.error("invalid call rate")
    selected={}
    for entry in args.pop:
        if "=" not in entry:ap.error("--pop requires NAME=sample-list.txt")
        name,path=entry.split("=",1)
        if not name or name in selected:ap.error("duplicate/empty population")
        selected[name]=read_lines(path)
        if not selected[name]:ap.error("empty population list: "+name)
    all_samples=[x for xs in selected.values() for x in xs]
    if len(set(all_samples))!=len(all_samples):ap.error("sample assigned to more than one population")
    mask=mask_load(args.callable_bed,regs["distal"][0])
    if not args.callable_bed:
        print("[WARNING] no callable BED: region physical length used as denominator (not StrictMask callable bp)",flush=True)
    data=read_vcf(args.vcf,selected,regs,args.callrate,mask)
    out=Path(args.out);rate=args.callrate
    if args.only in ("all","diversity"):
        res=[]
        for region,rows in data.items():
            l=callable_len(regs[region][1],regs[region][2],mask)
            for pop in selected:
                n,raw,pi,theta,taj,s=diversity(rows,pop,l,args.projection,rate)
                res.append((region,pop,len(selected[pop]),l,n,raw,args.projection,s,pi,theta,taj))
        write(out/"diversity.tsv",["region","population","sample_n","callable_bp","usable_polymorphic_snps","raw_segregating_snps","projection_n_alleles","projected_S","pi_per_bp","thetaW_per_bp","Tajima_D"],res)
    if args.only in ("all","ld"):
        res=[]
        for region,rows in data.items():
            cap=args.distal_ld_max if region=="distal" else args.proximal_ld_max
            for pop in selected:
                n,v=ld(rows,pop,None,cap,rate,args.ld_maf);res.append((region,pop,cap,n,v))
        write(out/"within_ld.tsv",["region","population","distance_lt_bp","n_pairs","mean_corrected_r2"],res)
    if args.only in ("all","across_ld"):
        res=[]
        for pop in selected:
            n,v=cross_ld(data["distal"],data["proximal"],pop,rate,args.ld_maf);res.append((pop,n,v))
        write(out/"across_ld.tsv",["population","n_pairs","mean_corrected_r2"],res)
    if args.only in ("all","fst"):
        res=[]
        pops=list(selected)
        for region,rows in data.items():
            for ref in pops[1:]:
                n,val=fst(rows,pops[0],ref,rate);res.append((region,pops[0],ref,n,val))
        write(out/"hudson_fst.tsv",["region","query","reference","shared_valid_snps","Hudson_FST"],res)
    print("[OK] focal summaries written to",out)

if __name__=="__main__":
    main()
