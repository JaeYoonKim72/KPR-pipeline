#!/usr/bin/env python3
"""Create aligned, GT-only phased VCFs at the same POS/REF/ALT SNP sites.

Used for the query and reference inputs of selscan --xpnsl.
Reads only short flank VCFs; missing genotypes remain missing.
"""
import argparse
import gzip

def op(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def load(path):
    meta=[]; samples=[]; rows={}
    with op(path) as fh:
        for line in fh:
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#CHROM"):
                samples=line.rstrip("\n").split("\t")[9:]
            elif line.strip() and not line.startswith("#"):
                f=line.rstrip("\n").split("\t")
                if len(f)<10: continue
                key=(f[0],int(f[1]),f[3],f[4])
                if key in rows:
                    raise ValueError("duplicate variant: "+str(key))
                ix=f[8].split(":").index("GT")
                gt=[c.split(":")[ix] for c in f[9:]]
                rows[key]=f[:8]+["GT"]+gt
    if not samples:raise ValueError("no samples in VCF: "+path)
    return meta,samples,rows

def output(path,meta,names,records):
    with open(path,"w") as fh:
        fh.writelines(meta)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(names)+"\n")
        for r in records: fh.write("\t".join(r)+"\n")

def main():
    p=argparse.ArgumentParser()
    for a in ("query","reference-vcf","query-label","reference-label","output","query-list","reference-list","query-matched","reference-matched"):
        p.add_argument("--"+a,required=True)
    a=p.parse_args()
    m1,n1,r1=load(a.query)
    m2,n2,r2=load(a.reference_vcf)
    keys=sorted(set(r1)&set(r2))
    q=[r1[k] for k in keys]
    r=[r2[k] for k in keys]
    output(a.query_matched,m1,n1,q)
    output(a.reference_matched,m2,n2,r)
    output(a.output,m1,[a.query_label+"__"+n for n in n1]+[a.reference_label+"__"+n for n in n2],
           [r1[k]+r2[k][9:] for k in keys])
    with open(a.query_list,"w") as f:
        f.writelines(a.query_label+"__"+n+"\n" for n in n1)
    with open(a.reference_list,"w") as f:
        f.writelines(a.reference_label+"__"+n+"\n" for n in n2)
    print("MATCHED\t"+str(len(keys)))

if __name__=="__main__":
    main()
