#!/usr/bin/env bash

set -euo pipefail


THREADS=${THREADS:-8}

SCRIPT_DIR=${SCRIPT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}
REF_DIR=${REF_DIR:-$(pwd)/reference}
LIST_DIR=${LIST_DIR:-$(pwd)/lists}
OUT_DIR=${OUT_DIR:-$(pwd)/results/gene_annotation}
TMP_DIR=${TMP_DIR:-}
LOG_OUT=${LOG_OUT:-}

RUN_MODE=${RUN_MODE:-tools}
PYTHON_BIN=${PYTHON_BIN:-${PYTHON:-python3}}
AUGUSTUS_SPECIES=${AUGUSTUS_SPECIES:-human}

GENCODE_GTF=${GENCODE_GTF:-${REF_DIR}/gencode.v46.annotation.gtf.gz}

SWISSPROT_FA=${SWISSPROT_FA:-${REF_DIR}/uniprot_sprot.human.fasta}

REFERENCE_FA=${REFERENCE_FA:-}

CHM13_ANNOTATION=${CHM13_ANNOTATION:-}
SATELLITE_LIBRARY=${SATELLITE_LIBRARY:-}
PROTEIN_GENE_MAP=${PROTEIN_GENE_MAP:-}

SAMPLE_TABLE=${SAMPLE_TABLE:-${LIST_DIR}/samples.tsv}

MIN_PROTEIN_IDENTITY=${MIN_PROTEIN_IDENTITY:-95}
MIN_BLASTP_IDENTITY=${MIN_BLASTP_IDENTITY:-90}
MIN_EXONERATE_IDENTITY=${MIN_EXONERATE_IDENTITY:-90}
MIN_PROTEIN_COVERAGE=${MIN_PROTEIN_COVERAGE:-80}
MIN_LIFTOFF_OVERLAP=${MIN_LIFTOFF_OVERLAP:-80}
MIN_EXON_SUPPORT=${MIN_EXON_SUPPORT:-80}
MIN_INTRON_SUPPORT=${MIN_INTRON_SUPPORT:-80}
MIN_INTRON_MATCH=${MIN_INTRON_MATCH:-80}
MIN_UNANNOTATED_OVERLAP=${MIN_UNANNOTATED_OVERLAP:-95}
EXONERATE_BATCH_BP=${EXONERATE_BATCH_BP:-4000000}


ENGINE_ARGS=()

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${LOG_OUT}/pipeline.log"
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 2
}

parse_options() {
  while [[ $# -gt 0 ]]; do
    if [[ "$1" == --*=* ]]; then
      set -- "${1%%=*}" "${1#*=}" "${@:2}"
    fi

    case "$1" in
      -h|--help) usage; exit 0 ;;
      --samples|--reference-gtf|--proteins|--reference-fasta|--chm13-annotation|\
      --satellite-library|--protein-gene-map|--out|--mode|--threads|--species|\
      --protein-identity|--blastp-identity|--exonerate-identity|--support-coverage|--liftoff-overlap|\
      --min-exon-support|--min-intron-support|--intron-match|--unannotated-overlap|\
      --exonerate-batch-bp)
        [[ $# -ge 2 && "$2" != --* ]] || die "missing value for $1"
        ;;
      *) die "unknown option: $1 (use --help)" ;;
    esac

    case "$1" in
      --samples)              SAMPLE_TABLE=$2 ;;
      --reference-gtf)        GENCODE_GTF=$2 ;;
      --proteins)             SWISSPROT_FA=$2 ;;
      --reference-fasta)      REFERENCE_FA=$2 ;;
      --chm13-annotation)     CHM13_ANNOTATION=$2 ;;
      --satellite-library)    SATELLITE_LIBRARY=$2 ;;
      --protein-gene-map)     PROTEIN_GENE_MAP=$2 ;;
      --out)                  OUT_DIR=$2 ;;
      --mode)                 RUN_MODE=$2 ;;
      --threads)              THREADS=$2 ;;
      --species)              AUGUSTUS_SPECIES=$2 ;;
      --protein-identity)     MIN_PROTEIN_IDENTITY=$2 ;;
      --blastp-identity)      MIN_BLASTP_IDENTITY=$2 ;;
      --exonerate-identity)   MIN_EXONERATE_IDENTITY=$2 ;;
      --support-coverage)     MIN_PROTEIN_COVERAGE=$2 ;;
      --liftoff-overlap)      MIN_LIFTOFF_OVERLAP=$2 ;;
      --min-exon-support)     MIN_EXON_SUPPORT=$2 ;;
      --min-intron-support)   MIN_INTRON_SUPPORT=$2 ;;
      --intron-match)         MIN_INTRON_MATCH=$2 ;;
      --unannotated-overlap)  MIN_UNANNOTATED_OVERLAP=$2 ;;
      --exonerate-batch-bp)    EXONERATE_BATCH_BP=$2 ;;
    esac
    shift 2
  done
}

configure_runtime() {
  command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python not found: ${PYTHON_BIN}"
  [[ "${RUN_MODE}" == tools || "${RUN_MODE}" == outputs ]] || \
    die "RUN_MODE must be tools or outputs"

  TMP_DIR=${TMP_DIR:-${OUT_DIR}/tmp}
  LOG_OUT=${LOG_OUT:-${OUT_DIR}/logs}
  mkdir -p "${OUT_DIR}" "${LOG_OUT}" "${TMP_DIR}"

  export TMPDIR="$(cd -- "${TMP_DIR}" && pwd)"
  export CNG_SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"

  if [[ -z "${REFERENCE_FA}" && -f "${REF_DIR}/GRCh38.fa" ]]; then
    REFERENCE_FA=${REF_DIR}/GRCh38.fa
  fi

  ENGINE_ARGS=(
    --samples "${SAMPLE_TABLE}"
    --reference-gtf "${GENCODE_GTF}"
    --proteins "${SWISSPROT_FA}"
    --out "${OUT_DIR}"
    --mode "${RUN_MODE}"
    --threads "${THREADS}"
    --species "${AUGUSTUS_SPECIES}"
    --protein-identity "${MIN_PROTEIN_IDENTITY}"
    --blastp-identity "${MIN_BLASTP_IDENTITY}"
    --exonerate-identity "${MIN_EXONERATE_IDENTITY}"
    --support-coverage "${MIN_PROTEIN_COVERAGE}"
    --liftoff-overlap "${MIN_LIFTOFF_OVERLAP}"
    --min-exon-support "${MIN_EXON_SUPPORT}"
    --min-intron-support "${MIN_INTRON_SUPPORT}"
    --intron-match "${MIN_INTRON_MATCH}"
    --unannotated-overlap "${MIN_UNANNOTATED_OVERLAP}"
    --exonerate-batch-bp "${EXONERATE_BATCH_BP}"
  )
  [[ -z "${REFERENCE_FA}" ]] || ENGINE_ARGS+=(--reference-fasta "${REFERENCE_FA}")
  [[ -z "${CHM13_ANNOTATION}" ]] || ENGINE_ARGS+=(--chm13-annotation "${CHM13_ANNOTATION}")
  [[ -z "${SATELLITE_LIBRARY}" ]] || ENGINE_ARGS+=(--satellite-library "${SATELLITE_LIBRARY}")
  [[ -z "${PROTEIN_GENE_MAP}" ]] || ENGINE_ARGS+=(--protein-gene-map "${PROTEIN_GENE_MAP}")
  return 0
}

run_stage() {
  run_engine "$1" 2> >(tee -a "${LOG_OUT}/pipeline.log" >&2) | \
    tee -a "${LOG_OUT}/pipeline.log"
}


softmask() {
  log "step 1: repeat soft masking"

  run_stage mask

  log "step 1 finished"
}


predict_orfs() {
  log "step 2: AUGUSTUS coding transcript prediction"

  run_stage augustus

  log "step 2 finished"
}


protein_support() {
  log "step 3: Exonerate protein support"

  run_stage exonerate

  log "step 3 finished"
}


gene_annotation() {
  log "step 4: Liftoff annotation and extra-copy detection"

  run_stage liftoff

  log "step 4 finished"
}


compare_annotations() {
  log "step 5: comparison with Liftoff extra-copy transcripts"

  run_stage compare

  log "step 5 finished"
}


search_unannotated() {
  log "step 6: BLASTp review of Liftoff-unannotated candidates"

  run_stage blastp

  log "step 6 finished"
}


filter_flagger() {
  log "step 7: assembly-error filtering and final annotation"

  run_stage flagger

  log "step 7 finished"
}


make_tables() {
  log "tables: annotation ratios and CNG summaries"

  run_stage tables

  log "tables finished"
}

check_inputs() {
  log "checking inputs and required software"
  run_stage check
}

run_all() {
  log "running steps 1-7 and table generation"

  run_stage all

  log "all steps finished"
}


run_engine() {
  "${PYTHON_BIN}" - "$@" "${ENGINE_ARGS[@]}" <<'PY'
import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote, unquote

if sys.version_info < (3,9):
    raise SystemExit('Python >=3.9 is required')
STEPS = ['mask','augustus','exonerate','liftoff','compare','blastp','flagger','tables']
GROUPS = ['KPR','HPRC-EAS','HPRC-nEAS']
BLAST_FIELDS = 'qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send evalue bitscore stitle'.split()
LABELS = {'hap':'Hap','haploid':'Hap','dup':'Dup','duplicated':'Dup','col':'Col','collapsed':'Col','err':'Err','erroneous':'Err','unk':'Unk','unknown':'Unk'}
EXCLUDED = {'Q9UN81':'TRANSPOSABLE_ELEMENT_LINE1_ORF1'}
REVIEW = {'Q8IX06':'UNCERTAIN_REFERENCE_REXO1L1P','Q96PT4':'REPEAT_ASSOCIATED_DUX_FAMILY','O43812':'REPEAT_ASSOCIATED_DUX_FAMILY'}
HASHES = {}

def need(ok,msg):
    if not ok: raise ValueError(msg)
def log(s): print(s,flush=True)
def sha(p):
    p=Path(p); key=(str(p.resolve()),p.stat().st_size,p.stat().st_mtime_ns)
    if key not in HASHES:
        h=hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
        HASHES[key]=h.hexdigest()
    return HASHES[key]
def jread(p): return json.loads(Path(p).read_text())
def jsave(p,obj):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.part');q.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n');q.replace(p)
def table(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    for ext,sep,enc in [('tsv','\t','utf-8'),('csv',',','utf-8-sig')]:
        q=Path(str(p)+'.'+ext)
        with q.open('w',encoding=enc,newline='') as f:
            w=csv.DictWriter(f,fields,delimiter=sep,extrasaction='ignore');w.writeheader();w.writerows(rows)
def attr(text):
    if re.match(r'^\s*[^\s=;]+=',text):
        return {k.strip():unquote(v.strip()) for x in text.split(';') if '=' in x for k,v in [x.split('=',1)]}
    return {a:b or c.strip() for a,b,c in re.findall(r'([^\s;]+)\s+(?:"([^"]*)"|([^;]*))\s*;?',text)}
def merged(parts):
    out=[]
    for s,e in sorted(set(map(tuple,parts))):
        need(0<=s<e,'Invalid interval')
        if out and s<=out[-1][1]: out[-1][1]=max(out[-1][1],e)
        else: out.append([s,e])
    return out
def length(p): return sum(e-s for s,e in merged(p))
def shared(a,b): return sum(max(0,min(e,y)-max(s,x)) for s,e in merged(a) for x,y in merged(b))
def introns(parts):
    x=merged(parts);return {(a[1],b[0]) for a,b in zip(x,x[1:]) if a[1]<b[0]}
def ratios(a,b):
    n=shared(a,b);return (100*n/length(a),100*n/length(b)) if a and b else (0.,0.)
def chain(a,b):
    x,y=introns(a),introns(b)
    if not x and not y: return 100.
    if not x or not y: return 0.
    return min(100*len(x&y)/len(x),100*len(x&y)/len(y))
def canonical(identifier,copy):
    suf='_'+str(copy)
    return identifier[:-len(suf)] if copy>0 and identifier.endswith(suf) else identifier

def gff(path,liftoff=False):
    genes={}; tx={}; features=[]
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for no,line in enumerate(f,1):
            if line.startswith('##FASTA'):break
            if not line.strip() or line.startswith('#'):continue
            a=line.rstrip('\r\n').split('\t')
            need(len(a)==9,f'{path}:{no}: expected GFF/GTF9')
            at=attr(a[8]); kind=a[2].lower();s,e=int(a[3])-1,int(a[4]);need(0<=s<e,f'{path}:{no}: coordinates')
            base=dict(contig=unquote(a[0]),start0=s,end0=e,strand=a[6])
            if kind=='gene':
                gid=at.get('ID') or at.get('gene_id');need(gid,f'{path}:{no}: missing gene ID')
                copy=at.get('extra_copy_number','0');need(copy.isdigit(),f'{path}:{no}: bad copy index')
                if liftoff:need('extra_copy_number' in at,f'{path}:{no}: Liftoff extra_copy_number required')
                row=dict(base,id=gid,gene=at.get('gene_name') or at.get('gene') or at.get('Name') or '',
                         biotype=at.get('gene_type') or at.get('gene_biotype') or at.get('biotype') or '',copy=int(copy),attrs=at)
                need(gid not in genes or genes[gid]==row,f'{path}: conflicting duplicate gene {gid}')
                genes[gid]=row
            elif kind in ('mrna','transcript','primary_transcript') or kind.endswith('rna'):
                tid=at.get('ID') or at.get('transcript_id');gid=at.get('Parent') or at.get('gene_id')
                need(tid and gid and ',' not in gid,f'{path}:{no}: transcript parent required')
                row=dict(base,id=tid,gene_id=gid,biotype=at.get('transcript_type') or at.get('transcript_biotype') or '',cds=[],starts=[],stops=[],exons=[])
                need(tid not in tx or tx[tid]==row,f'{path}: duplicate transcript {tid}')
                tx[tid]=row
            elif kind in ('cds','exon','start_codon','stop_codon'):
                parents=at.get('Parent') or at.get('transcript_id') or at.get('gene_id')
                need(parents,f'{path}:{no}: missing feature parent')
                features.append((base,kind,a[7],parents.split(','),at))
    for base,kind,phase,parents,at in features:
        for tid in parents:
            if tid not in tx:
                gid=at.get('gene_id') or (tid if tid in genes else None)
                if gid and gid in genes:
                    tx[tid]=dict(genes[gid],id=tid,gene_id=gid,cds=[],starts=[],stops=[],exons=[])
                else: raise ValueError(f'{path}: feature parent not resolved: {tid}')
            t=tx[tid];need(t['contig']==base['contig'] and t['strand']==base['strand'],f'{path}: transcript/feature coordinates disagree')
            t[{'cds':'cds','exon':'exons','start_codon':'starts','stop_codon':'stops'}[kind]].append([base['start0'],base['end0'],phase] if kind=='cds' else [base['start0'],base['end0']])
    for t in tx.values():
        need(t['gene_id'] in genes,f'{path}: missing parent gene {t["gene_id"]}')
        gene=genes[t['gene_id']]
        need(gene['contig']==t['contig'] and gene['strand']==t['strand'] and gene['start0']<=t['start0']<t['end0']<=gene['end0'],f'{path}: transcript outside gene')
        t['cds']=sorted({tuple(c) for c in t['cds']});t['starts']=merged(t['starts']);t['stops']=merged(t['stops']);t['exons']=merged(t['exons'])
        need(all(gene['start0']<=x[0]<x[1]<=gene['end0'] for x in t['cds']),f'{path}: CDS outside gene')
    return genes,tx

class Fasta:
    """Seekable plain FASTA with a small in-memory index, not a whole-genome dict."""
    def __init__(self,path):
        self.path=Path(path);self.idx={};self.f=self.path.open('rb');name=None;last=0
        while True:
            off=self.f.tell();line=self.f.readline()
            if not line:break
            if line.startswith(b'>'):
                name=line[1:].split()[0].decode();need(name not in self.idx,f'Duplicate FASTA ID {name}')
                self.idx[name]=[0,self.f.tell(),0,0];last=0
            else:
                need(name is not None,'Sequence before FASTA header')
                seq=line.rstrip(b'\r\n');need(seq and not re.search(rb'\s',seq),'Blank/whitespace FASTA sequence line')
                row=self.idx[name]
                if row[2]==0:row[2:]=[len(seq),len(line)]
                elif last:need(last==row[2],f'Irregular FASTA wrapping: {name}')
                row[0]+=len(seq);last=len(seq)
        need(self.idx and all(x[0]>0 for x in self.idx.values()),f'Empty FASTA {path}')
    def get(self,seq,s,e):
        need(seq in self.idx,f'Contig absent from FASTA: {seq}')
        n,off,bases,width=self.idx[seq];need(0<=s<e<=n,f'FASTA interval out of range: {seq}:{s}-{e}')
        self.f.seek(off+(s//bases)*width+s%bases)
        raw=self.f.read((e-s)+((e//bases)-(s//bases)+2)*(width-bases))
        val=raw.replace(b'\n',b'').replace(b'\r',b'')[:e-s].decode()
        need(len(val)==e-s,'Short FASTA interval');return val
    def close(self):self.f.close()

def plain(source,destination):
    source=Path(source);destination=Path(destination)
    if source.suffix!='.gz':return source
    if not destination.exists():
        tmp=destination.with_suffix(destination.suffix+'.part')
        with gzip.open(source,'rb') as f,tmp.open('wb') as out:shutil.copyfileobj(f,out)
        tmp.replace(destination)
    return destination

def read_proteins(path,gene_map=None):
    refs={};seq=[];head=None
    def add():
        if head is None:return
        pid=head.split()[0];aa=''.join(seq).rstrip('*').upper();need(pid not in refs and aa,f'Invalid protein {pid}')
        m=re.search(r'\bGN=([^\s]+)',head);accession=pid.split('|')[1] if pid.startswith(('sp|','tr|')) else pid
        refs[pid]=dict(id=pid,accession=accession,gene=m.group(1) if m else '',aa=aa)
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        for line in f:
            if line.startswith('>'):add();head=line[1:].strip();seq=[]
            elif line.strip():seq.append(line.strip())
    add();need(refs,'Empty reference protein FASTA')
    if gene_map:
        for r in table(gene_map):
            need(r['protein_id'] in refs,'Protein-map ID absent from FASTA');refs[r['protein_id']]['gene']=r['gene']
    aliases={}
    for pid,r in refs.items():
        for key in (pid,r['accession']):
            need(key not in aliases or aliases[key]==pid,'Ambiguous protein alias');aliases[key]=pid
    return refs,aliases

def translation(t,fa):
    from Bio.Seq import Seq
    blocks=sorted(t['cds'],reverse=t['strand']=='-');parts=[[x[0],x[1]] for x in blocks]
    need(blocks,f'No CDS for {t["id"]}')
    dna='';phaseok=True
    for i,(s,e,p) in enumerate(blocks):
        p=0 if p=='.' else int(p)
        if i==0:phaseok=(p==0)
        else:phaseok=phaseok and p==((3-len(dna)%3)%3)
        seg=fa.get(t['contig'],s,e)
        dna+=str(Seq(seg).reverse_complement()) if t['strand']=='-' else seg
    if t['stops'] and not dna[-3:].upper() in ('TAA','TAG','TGA'):
        remaining=[]
        for s,e in t['stops']:
            for pos in range(s,e):
                if not any(a<=pos<b for a,b in parts):remaining.append(pos)
        expected=[max(e for s,e in parts)+i for i in range(3)] if t['strand']=='+' else [min(s for s,e in parts)-3+i for i in range(3)]
        if sorted(remaining)==expected:
            s,e=min(remaining),max(remaining)+1;seg=fa.get(t['contig'],s,e)
            dna+=str(Seq(seg).reverse_complement()) if t['strand']=='-' else seg
    triplet=len(dna)%3==0
    aa=str(Seq(dna[:len(dna)//3*3]).translate()) if dna else ''
    stop=aa.endswith('*');protein=aa[:-1] if stop else aa
    complete=bool(triplet and phaseok and dna[:3].upper()=='ATG' and stop)
    internal=protein.count('*');ambiguous=any(x not in 'ACDEFGHIKLMNPQRSTVWY*' for x in protein)
    coding=merged(parts)
    if stop and dna[-3:].upper() in ('TAA','TAG','TGA'):
        n=3 if len(dna)==sum(e-s for s,e in parts) else 0
        order=list(reversed(coding)) if t['strand']=='+' else coding[:]
        for piece in order:
            cut=min(n,piece[1]-piece[0]);n-=cut
            if t['strand']=='+':piece[1]-=cut
            else:piece[0]+=cut
            if n==0:break
        coding=[x for x in coding if x[0]<x[1]]
    lower=100*sum(x in 'acgt' for x in dna)/len(dna) if dna else 0.
    return dict(aa=protein,complete=complete,internal_stops=internal,ambiguous=ambiguous,
                lowercase_cds_pct=lower,cds_parts=merged(parts),protein_parts=coding,cds_nt=len(dna))

def aa_compare(query,pred):
    from Bio import Align
    from Bio.Align import substitution_matrices
    if not query or not pred:return dict(identity=0.,query_coverage=0.,pred_coverage=0.,score=0.)
    aligner=Align.PairwiseAligner();aligner.mode='local';aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
    alns=aligner.align(query,pred)
    if not alns:return dict(identity=0.,query_coverage=0.,pred_coverage=0.,score=0.)
    aln=alns[0];c=aln.coordinates;ident=0;columns=0;qn=0;pn=0
    for i in range(c.shape[1]-1):
        a,b=int(c[0,i]),int(c[0,i+1]);x,y=int(c[1,i]),int(c[1,i+1]);dq,dp=b-a,y-x
        columns+=max(dq,dp);qn+=dq;pn+=dp
        if dq and dp:ident+=sum(u==v for u,v in zip(query[a:b],pred[x:y]))
    return dict(identity=100*ident/columns if columns else 0.,query_coverage=100*qn/len(query),pred_coverage=100*pn/len(pred),score=float(aln.score))

def exo_alignments(path,aliases):
    rows=[];current=None
    with Path(path).open() as f:
        for line in f:
            if 'END OF GFF DUMP' in line:
                if current:rows.append(current);current=None
                continue
            if not line.strip() or line.startswith('#'):continue
            a=line.rstrip('\n').split('\t')
            if len(a)!=9 or not a[1].startswith('exonerate:protein2genome'):continue
            if a[2].lower()=='gene':
                if current:rows.append(current)
                m=re.search(r'\bsequence\s+(\S+)\s*(?:;|$)',a[8]);pi=re.search(r'\bidentity\s+([\d.]+)',a[8])
                need(m and pi,'Exonerate gene lacks sequence/identity attributes')
                q=m.group(1);need(q in aliases,f'Exonerate query absent from Swiss-Prot: {q}')
                current=dict(id='exo'+str(len(rows)+1),target=a[0],strand=a[6],query=aliases[q],identity=float(pi.group(1)),score=float(a[5]),cds=[])
            elif a[2].lower()=='cds':
                need(current is not None and current['target']==a[0] and current['strand']==a[6],'Exonerate CDS/gene mismatch')
                current['cds'].append([int(a[3])-1,int(a[4])])
    if current:rows.append(current)
    for r in rows:
        r['cds']=merged(r['cds']);need(r['cds'],'Exonerate alignment without CDS')
    return rows

def command(argv,folder,stdout=None):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    log('[COMMAND] '+shlex.join(list(map(str,argv))))
    with (folder/'commands.jsonl').open('a') as f:f.write(json.dumps(list(map(str,argv)))+'\n')
    with (folder/'stderr.log').open('ab') as err:
        if stdout:
            tmp=Path(str(stdout)+'.part')
            with tmp.open('wb') as out:r=subprocess.run(list(map(str,argv)),stdout=out,stderr=err)
            need(r.returncode==0,f'Tool failed ({r.returncode}): see {folder}/stderr.log')
            tmp.replace(stdout)
        else:
            with (folder/'stdout.log').open('ab') as out:r=subprocess.run(list(map(str,argv)),stdout=out,stderr=err)
            need(r.returncode==0,f'Tool failed ({r.returncode}): see {folder}/stderr.log')

def inputs(a):
    manifest=a.samples.resolve();rows=table(manifest);required={'sample','group','fasta','flagger_bed'}
    need(rows and required<=rows[0].keys(),'samples.tsv needs sample, group, fasta, flagger_bed')
    need(len({r['sample'] for r in rows})==len(rows),'Duplicate sample ID')
    paths={'fasta','flagger_bed','masked_fasta','augustus_gff','exonerate_gff','liftoff_gff','blastp_tsv','contig_map','augustus_hints','extrinsic_cfg'}
    for r in rows:
        need(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',r['sample']) and r['group'] in GROUPS,'Unsafe sample ID or unknown group')
        for k in paths:
            val=r.get(k,'').strip()
            if val and val!='.':
                p=Path(val);r[k]=str((manifest.parent/p).resolve() if not p.is_absolute() else p.resolve())
                need(Path(r[k]).is_file(),f'{r["sample"]}: input not found: {r[k]}')
            else:r[k]=''
        need(r['fasta'] and r['flagger_bed'],f'{r["sample"]}: fasta and flagger_bed paths are required')
        if a.mode=='outputs':need(all(r[k] for k in ('masked_fasta','augustus_gff','exonerate_gff','liftoff_gff','blastp_tsv')),f'{r["sample"]}: outputs mode requires masked_fasta, augustus_gff, exonerate_gff, liftoff_gff and blastp_tsv (empty file allowed for no hits)')
    return rows

class Pipeline:
    def __init__(self,a):
        self.a=a;self.rows=inputs(a);self.root=a.out.resolve();self.root.mkdir(parents=True,exist_ok=True)
        self.ref,self.alias=read_proteins(a.proteins,a.protein_gene_map)
        self.reference_genes,self.reference_tx=gff(a.reference_gtf)
        self.refcounts=defaultdict(set)
        for g in self.reference_genes.values():
            if g['gene']:self.refcounts[g['gene']].add(g['id'])
        need(self.refcounts,'Reference annotation contains no named gene features')
        checks={str(p.resolve()):sha(p) for p in [a.samples,a.reference_gtf,a.proteins]}
        for name in ('reference_fasta','protein_gene_map','satellite_library','chm13_annotation'):
            val=getattr(a,name)
            if val:checks[str(val.resolve())]=sha(val)
        for r in self.rows:
            for key,val in r.items():
                if val and key not in ('sample','group') and Path(val).is_file():checks[val]=sha(val)
        spec=dict(schema=1,script_sha256=sha(os.environ['CNG_SCRIPT_PATH']),parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items() if k not in ('step','out')},inputs=checks,samples=self.rows)
        existing=self.root/'run.json'
        if existing.exists():need(jread(existing)==spec,'Output directory belongs to different inputs/options/code. Choose a NEW --out directory.')
        else:jsave(existing,spec)
        self.signature=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
        shared=Counter(r['flagger_bed'] for r in self.rows)
        self.shared_beds={path for path,n in shared.items() if n>1}
        (self.root/'references').mkdir(exist_ok=True)
        self.protein_fasta=plain(a.proteins,self.root/'references/human_swissprot.fa')
    def folder(self,r,step):return self.root/'samples'/r['sample']/step
    def data(self,r,step):return jread(self.folder(r,step)/'result.json')
    def verify(self,r,step):
        d=self.folder(r,step);p=d/'completed.json'
        need(p.is_file(),f'{r["sample"]}: run {step} first')
        m=jread(p);need(m['signature']==self.signature and all((d/k).is_file() and sha(d/k)==v for k,v in m['outputs'].items()),f'Changed/incomplete prerequisite: {d}')
    def run(self,step):
        if step=='tables':return self.tables()
        previous=STEPS[:STEPS.index(step)]
        for r in self.rows:
            for pre in previous:self.verify(r,pre)
            d=self.folder(r,step);marker=d/'completed.json'
            if marker.exists():
                m=jread(marker);need(m['signature']==self.signature and all((d/p).is_file() and sha(d/p)==v for p,v in m['outputs'].items()),f'Changed/incomplete checkpoint: {d}')
                log('[REUSE] '+r['sample']+' '+step);continue
            d.mkdir(parents=True,exist_ok=True);log('[STEP] '+r['sample']+' '+step)
            result=getattr(self,'do_'+step)(r,d)
            jsave(d/'result.json',result)
            hashes={str(p.relative_to(d)):sha(p) for p in d.rglob('*') if p.is_file() and p.name!='completed.json' and not p.name.endswith('.part')}
            jsave(marker,dict(signature=self.signature,outputs=hashes))
    def do_mask(self,r,d):
        original=plain(r['fasta'],d/'original.fa')
        if r['masked_fasta']:masked=plain(r['masked_fasta'],d/'cached_masked.fa')
        else:
            local=d/'assembly.fa';shutil.copyfile(original,local)
            command(['RepeatMasker','-species','human','-xsmall','-pa',str(self.a.threads),'-dir',d/'human',local],d/'human')
            masked=d/'human'/'assembly.fa.masked';need(masked.is_file(),'RepeatMasker did not produce .masked FASTA')
            if self.a.satellite_library:
                command(['RepeatMasker','-lib',self.a.satellite_library,'-xsmall','-pa',str(self.a.threads),'-dir',d/'satellite',masked],d/'satellite')
                masked=d/'satellite'/'assembly.fa.masked.masked';need(masked.is_file(),'Satellite RepeatMasker output missing')
        fa=Fasta(original);fm=Fasta(masked)
        need(set(fa.idx)==set(fm.idx),'Masking changed FASTA contig IDs')
        for seq in fa.idx:
            n=fa.idx[seq][0];need(n==fm.idx[seq][0],'Masking changed sequence length')
            for s in range(0,n,1000000):need(fa.get(seq,s,min(n,s+1000000)).upper()==fm.get(seq,s,min(n,s+1000000)).upper(),'Masked FASTA is not the same assembly sequence')
        fa.close();fm.close();return dict(fasta=str(masked.resolve()))
    def do_augustus(self,r,d):
        path=d/'augustus.gff3'
        if r['augustus_gff']:shutil.copyfile(r['augustus_gff'],path)
        else:
            cmd=['augustus','--species='+self.a.species,'--softmasking=1','--gff3=on']
            if r['augustus_hints']:cmd+=['--hintsfile='+r['augustus_hints']]
            if r['extrinsic_cfg']:cmd+=['--extrinsicCfgFile='+r['extrinsic_cfg']]
            command(cmd+[self.data(r,'mask')['fasta']],d,path)
        genes,tx=gff(path);fa=Fasta(self.data(r,'mask')['fasta']);models={};audit=[]
        for tid,t in tx.items():
            reasons=[]
            if not t['cds']:reasons=['NO_CDS'];x={}
            else:
                x=translation(t,fa)
                if not x['complete']:reasons.append('INCOMPLETE_OR_NONCANONICAL_ORF')
                if x['internal_stops']:reasons.append('INTERNAL_STOP')
                if x['ambiguous']:reasons.append('AMBIGUOUS_TRANSLATION')
            audit.append(dict(transcript=tid,gene=t['gene_id'],complete=x.get('complete',False),internal_stops=x.get('internal_stops',''),decision='PASS' if not reasons else 'EXCLUDE',reason=';'.join(reasons)))
            if not reasons:models[tid]=dict(t,**x)
        fa.close();write(d/'orf_review',audit,['transcript','gene','complete','internal_stops','decision','reason'])
        return dict(genes=genes,models=models,prediction_mode='hints' if r['augustus_hints'] else ('source_output' if r['augustus_gff'] else 'de_novo'))
    def do_exonerate(self,r,d):
        aug=self.data(r,'augustus');models=aug['models'];path=d/'exonerate.gff';windows={};native=not r['exonerate_gff']
        if not native:
            shutil.copyfile(r['exonerate_gff'],path);align=exo_alignments(path,self.alias)
            for x in align:x['contig']=x.pop('target')
        else:
            fa=Fasta(self.data(r,'mask')['fasta']);targets=[]
            for gid,g in sorted(aug['genes'].items()):
                if not any(t['gene_id']==gid for t in models.values()):continue
                wid='w'+str(len(windows)+1).zfill(7);windows[wid]=dict(contig=g['contig'],start0=g['start0'],end0=g['end0']);targets.append((wid,fa.get(g['contig'],g['start0'],g['end0'])))
            fa.close();jsave(d/'window_coordinates.json',windows);align=[];batch=[];total=0
            def flush():
                nonlocal batch,total
                if not batch:return
                ix=len(list(d.glob('part_*.fa')))+1;target=d/f'part_{ix:04d}.fa';raw=d/f'part_{ix:04d}.gff'
                with target.open('w') as f:
                    for wid,seq in batch:f.write('>'+wid+'\n'+seq+'\n')
                command(['exonerate','--model','protein2genome','--percent','80','--softmasktarget','no','--showvulgar','no','--showalignment','no','--showquerygff','no','--showtargetgff','yes','-q',self.protein_fasta,'-t',target],d/f'part_{ix:04d}',raw)
                for x in exo_alignments(raw,self.alias):
                    need(x['target'] in windows,'Unrecognized Exonerate target');w=windows[x.pop('target')];x['contig']=w['contig'];x['cds']=[[a+w['start0'],b+w['start0']] for a,b in x['cds']];x['id']=f'part{ix}.'+x['id'];align.append(x)
                batch=[];total=0
            for wid,seq in targets:
                if batch and total+len(seq)>self.a.exonerate_batch_bp:flush()
                batch.append((wid,seq));total+=len(seq)
            flush()
        byctg=defaultdict(list);fa_check=Fasta(self.data(r,'mask')['fasta'])
        for x in align:
            need(x['contig'] in fa_check.idx and all(0<=a<b<=fa_check.idx[x['contig']][0] for a,b in x['cds']),'Exonerate coordinates are not in assembly space; project fragment coordinates before replay')
            byctg[x['contig']].append(x)
        fa_check.close()
        supported={};audit=[];cache={}
        for tid,t in models.items():
            passing=[]
            for x in byctg[t['contig']]:
                if t['strand']!=x['strand'] or not shared(t['protein_parts'],x['cds']):continue
                rec=ratios(t['protein_parts'],x['cds']);key=(tid,x['query'])
                if key not in cache:cache[key]=aa_compare(self.ref[x['query']]['aa'],t['aa'])
                aa=cache[key];reasons=[]
                if x['identity']<self.a.exonerate_identity:reasons.append('EXONERATE_IDENTITY')
                if aa['identity']<self.a.protein_identity:reasons.append('AUGUSTUS_PROTEIN_IDENTITY')
                if min(aa['query_coverage'],aa['pred_coverage'])<self.a.support_coverage:reasons.append('PROTEIN_COVERAGE')
                if min(rec)<self.a.liftoff_overlap:reasons.append('AUGUSTUS_EXONERATE_CDS_OVERLAP')
                coding_exons=t['protein_parts']
                exon_support=100*sum(shared([part],x['cds'])/(part[1]-part[0])>=.8 for part in coding_exons)/len(coding_exons) if coding_exons else 0.
                ai,ei=introns(t['cds_parts']),introns(x['cds'])
                intron_support=100*len(ai&ei)/len(ai) if ai else 100.
                if exon_support<self.a.min_exon_support:reasons.append('EXON_ALIGNMENT_SUPPORT')
                if intron_support<self.a.min_intron_support:reasons.append('INTRON_ALIGNMENT_SUPPORT')
                strict=x['identity']>=self.a.protein_identity and min(rec)>=self.a.unannotated_overlap
                item=dict(transcript=tid,alignment=x['id'],reference=x['query'],exonerate_identity=x['identity'],protein_identity=aa['identity'],query_coverage=aa['query_coverage'],prediction_coverage=aa['pred_coverage'],augustus_CDS_overlap=rec[0],exonerate_CDS_overlap=rec[1],exon_alignment_support=exon_support,intron_alignment_support=intron_support,strict_unannotated_pass=int(strict and not reasons),decision='PASS' if not reasons else 'EXCLUDE',reason=';'.join(reasons))
                audit.append(item)
                if not reasons:passing.append(dict(item,score=aa['score']))
            if passing:supported[tid]=dict(t,support=sorted(passing,key=lambda z:(-z['protein_identity'],-z['score'],z['reference'])),strict_unannotated=any(x['strict_unannotated_pass'] for x in passing))
        fields=['transcript','alignment','reference','exonerate_identity','protein_identity','query_coverage','prediction_coverage','augustus_CDS_overlap','exonerate_CDS_overlap','exon_alignment_support','intron_alignment_support','strict_unannotated_pass','decision','reason']
        write(d/'protein_support_review',audit,fields)
        return dict(models=supported,alignment_count=len(align))
    def do_liftoff(self,r,d):
        path=d/'liftoff.gff3'
        if r['liftoff_gff']:shutil.copyfile(r['liftoff_gff'],path)
        else:
            need(self.a.reference_fasta,'--reference-fasta required for Liftoff')
            ref=plain(self.a.reference_fasta,d/'reference.fa');gtf=plain(self.a.reference_gtf,d/'reference.gtf')
            localref=d/'liftoff_reference.fa';localgtf=d/'liftoff_reference.gtf'
            shutil.copyfile(ref,localref);shutil.copyfile(gtf,localgtf)
            target=d/'target.fa';shutil.copyfile(self.data(r,'mask')['fasta'],target)
            command(['liftoff','-g',localgtf,'-o',path,'-u',d/'unmapped.txt','-dir',d/'intermediate','-p',str(self.a.threads),'-copies','-sc','0.95','-polish','-exclude_partial',target,localref],d)
            polished=Path(str(path)+'_polished');need(polished.is_file(),'Liftoff -polish did not produce polished annotation');polished.replace(path)
        genes,tx=gff(path,True);fa=Fasta(self.data(r,'mask')['fasta'])
        for g in genes.values():need(g['contig'] in fa.idx and g['end0']<=fa.idx[g['contig']][0],'Liftoff coordinates outside this assembly')
        fa.close();return dict(genes=genes,transcripts=tx)
    def do_compare(self,r,d):
        lif=self.data(r,'liftoff');models=self.data(r,'exonerate')['models'];aug=self.data(r,'augustus');kept={};reviews=[];new={}
        byctg=defaultdict(list)
        for t in models.values():byctg[t['contig']].append(t)
        for gid,g in lif['genes'].items():
            if g['copy']==0 or g['biotype']!='protein_coding':continue
            candidates=[]
            for lt in lif['transcripts'].values():
                if lt['gene_id']!=gid or not lt['cds']:continue
                parts=[[x[0],x[1]] for x in lt['cds']]
                for t in byctg[g['contig']]:
                    if t['strand']!=g['strand'] or not shared(parts,t['cds_parts']):continue
                    ov=ratios(parts,t['cds_parts']);ic=chain(parts,t['cds_parts']);passed=min(ov)>=self.a.liftoff_overlap and ic>=self.a.intron_match
                    reviews.append(dict(liftoff_gene=gid,liftoff_transcript=lt['id'],augustus_transcript=t['id'],liftoff_CDS_overlap=ov[0],augustus_CDS_overlap=ov[1],intron_concordance=ic,decision='PASS' if passed else 'EXCLUDE'))
                    if passed:candidates.append((min(ov),ic,t['support'][0]['protein_identity'],lt['id'],t['id']))
            if candidates:
                best=sorted(candidates,reverse=True)[0];kept[gid]=dict(g,chosen_liftoff_transcript=best[3],augustus_transcript=best[4])
        for tid,t in models.items():
            if not t['strict_unannotated']:continue
            gene=aug['genes'][t['gene_id']]
            hits=[g['id'] for g in lif['genes'].values() if g['contig']==gene['contig'] and max(g['start0'],gene['start0'])<min(g['end0'],gene['end0'])]
            if not hits:new[tid]=dict(t,start0=gene['start0'],end0=gene['end0'])
        write(d/'liftoff_structure_review',reviews,['liftoff_gene','liftoff_transcript','augustus_transcript','liftoff_CDS_overlap','augustus_CDS_overlap','intron_concordance','decision'])
        with (d/'unannotated.proteins.fa').open('w') as f:
            for tid,t in sorted(new.items()):f.write('>'+tid+'\n'+t['aa']+'\n')
        return dict(validated_extra=kept,new_candidates=new)
    def do_blastp(self,r,d):
        candidates=self.data(r,'compare')['new_candidates'];raw=d/'blastp.tsv'
        if r['blastp_tsv']:shutil.copyfile(r['blastp_tsv'],raw)
        elif candidates:
            db=d/'reference_db';command(['makeblastdb','-in',self.protein_fasta,'-dbtype','prot','-out',db],d)
            command(['blastp','-task','blastp-fast','-query',self.folder(r,'compare')/'unannotated.proteins.fa','-db',db,'-num_threads',str(self.a.threads),'-comp_based_stats','0','-seg','yes','-evalue','1e-5','-max_target_seqs','20','-max_hsps','1','-outfmt','6 '+' '.join(BLAST_FIELDS)],d,raw)
        else:raw.touch()
        hits=defaultdict(list)
        with raw.open() as f:
            for n,line in enumerate(f,1):
                if not line.strip() or line.startswith('#'):continue
                a=line.rstrip('\n').split('\t');need(len(a)==len(BLAST_FIELDS),f'{raw}:{n}: use the documented BLAST outfmt')
                b=dict(zip(BLAST_FIELDS,a));q=b['qseqid'];need(q in candidates,f'BLAST query not in current candidates: {q}')
                sid=b['sseqid'];sid=sid[4:] if sid.startswith('lcl|') else sid;need(sid in self.alias,f'BLAST reference absent from FASTA: {sid}')
                pid=self.alias[sid];qlen=int(b['qlen']);slen=int(b['slen']);qs,qe=int(b['qstart']),int(b['qend']);ss,se=int(b['sstart']),int(b['send'])
                need(qlen==len(candidates[q]['aa']) and slen==len(self.ref[pid]['aa']) and 1<=qs<=qe<=qlen and 1<=ss<=se<=slen,'BLAST lengths/coordinates disagree with input FASTAs')
                row=dict(query=q,reference=pid,identity=float(b['pident']),query_coverage=100*(qe-qs+1)/qlen,reference_coverage=100*(se-ss+1)/slen,bitscore=float(b['bitscore']),evalue=float(b['evalue']))
                need(0<=row['identity']<=100 and math.isfinite(row['bitscore']),'Invalid BLAST score')
                hits[q].append(row)
        eligible={};review=[]
        for q,t in candidates.items():
            hs=sorted(hits[q],key=lambda z:(-z['bitscore'],z['evalue'],z['reference']));best=hs[0] if hs else None;reasons=[];gene=''
            if best:
                ref=self.ref[best['reference']];gene=ref['gene'];acc=ref['accession']
                near=[x for x in hs if x['bitscore']>=best['bitscore']*.99]
                symbols={self.ref[x['reference']]['gene'] or 'UNRESOLVED:'+x['reference'] for x in near}
                if not gene:reasons.append('NO_REFERENCE_GENE_SYMBOL')
                if acc in EXCLUDED:reasons.append(EXCLUDED[acc])
                if acc in REVIEW:reasons.append(REVIEW[acc])
                if len(symbols)>1:reasons.append('AMBIGUOUS_NEAR_BEST_GENE_SYMBOLS')
                if len(hs)>=20 and hs[-1]['bitscore']>=best['bitscore']*.99:reasons.append('NEAR_BEST_HIT_LIMIT_REACHED')
                if best['identity']<self.a.blastp_identity:reasons.append('LOW_IDENTITY')
                if best['query_coverage']<95:reasons.append('LOW_QUERY_COVERAGE')
                if best['reference_coverage']<80:reasons.append('LOW_REFERENCE_COVERAGE')
                if best['evalue']>1e-5:reasons.append('EVALUE')
            else:reasons.append('NO_BLAST_HIT')
            if t['lowercase_cds_pct']>=50:reasons.append('REPEAT_DOMINATED_CDS')
            review.append(dict(transcript=q,gene=gene,reference=best['reference'] if best else '',identity=best['identity'] if best else '',query_coverage=best['query_coverage'] if best else '',reference_coverage=best['reference_coverage'] if best else '',lowercase_cds_pct=t['lowercase_cds_pct'],decision='KEEP' if not reasons else 'EXCLUDE',reason=';'.join(reasons)))
            if not reasons:eligible[q]=dict(t,gene=gene,reference=best['reference'])
        write(d/'candidate_review',review,['transcript','gene','reference','identity','query_coverage','reference_coverage','lowercase_cds_pct','decision','reason'])
        return dict(eligible=eligible)
    def do_flagger(self,r,d):
        lif=self.data(r,'liftoff');valid=self.data(r,'compare')['validated_extra'];new=self.data(r,'blastp')['eligible'];bed=Bed(r['flagger_bed'],r['contig_map'],r['sample'],r['flagger_bed'] in self.shared_beds);fa=Fasta(self.data(r,'mask')['fasta']);retained={};new_kept={};review=[]
        for gid,g in lif['genes'].items():
            eligible=g['copy']==0 or g['biotype']!='protein_coding' or gid in valid
            if not eligible:
                decision,reason,cov='EXCLUDE','FAILED_PROTEIN_OR_STRUCTURE_FILTER',Counter()
            else:decision,reason,cov=bed.assess(g['contig'],g['start0'],g['end0'],fa)
            review.append(dict(source='Liftoff',id=gid,gene=g['gene'],contig=g['contig'],start1=g['start0']+1,end1=g['end0'],decision=decision,reason=reason,Hap_bp=cov['Hap'],Issue_bp=cov['Issue'],Covered_bp=cov['Covered']))
            if decision=='KEEP':retained[gid]=g
        for tid,t in new.items():
            decision,reason,cov=bed.assess(t['contig'],t['start0'],t['end0'],fa)
            review.append(dict(source='AUGUSTUS_BLASTp',id=tid,gene=t['gene'],contig=t['contig'],start1=t['start0']+1,end1=t['end0'],decision=decision,reason=reason,Hap_bp=cov['Hap'],Issue_bp=cov['Issue'],Covered_bp=cov['Covered']))
            if decision=='KEEP':new_kept[tid]=t
        fa.close();write(d/'locus_review',review,['source','id','gene','contig','start1','end1','decision','reason','Hap_bp','Issue_bp','Covered_bp'])
        before=[g for gid,g in lif['genes'].items() if g['copy']==0 or g['biotype']!='protein_coding' or gid in valid]
        bcng,l1=cng(before,list(new.values()),self.refcounts,self.reference_genes)
        after,l2=cng(list(retained.values()),list(new_kept.values()),self.refcounts,self.reference_genes)
        write(d/'copy_number_review',l2,['gene','reference_loci','confirmed_loci','CNG','new_prediction_loci'])
        self.export_gff(d/'final_annotation.gff3',lif,retained,valid,new_kept)
        return dict(before_CNG=bcng,CNG=after,retained_liftoff=retained,retained_new=new_kept)
    def export_gff(self,path,lif,kept,valid,new):
        def line(f,seq,kind,s,e,strand,phase,at):
            st=';'.join(quote(str(k),safe='_:.-')+'='+quote(str(v),safe='_:.-|') for k,v in at.items())
            f.write('\t'.join(map(str,[seq,'KPR_pipeline',kind,s+1,e,'.',strand,phase,st]))+'\n')
        with path.open('w') as f:
            f.write('##gff-version 3\n')
            for ix,(gid,g) in enumerate(sorted(kept.items()),1):
                outid='L'+str(ix).zfill(7);ref=canonical(gid,g['copy'])
                line(f,g['contig'],'gene',g['start0'],g['end0'],g['strand'],'.',dict(ID=outid,gene_name=g['gene'],gene_type=g['biotype'],reference_gene_id=ref,source_gene_id=gid,extra_copy_number=g['copy']))
                txs=[t for t in lif['transcripts'].values() if t['gene_id']==gid and (gid not in valid or t['id']==valid[gid]['chosen_liftoff_transcript'])]
                for n,t in enumerate(txs,1):
                    oid=outid+'.t'+str(n);line(f,t['contig'],'mRNA' if t['cds'] else 'transcript',t['start0'],t['end0'],t['strand'],'.',dict(ID=oid,Parent=outid,reference_transcript_id=canonical(t['id'],g['copy'])))
                    for i,(s,e) in enumerate(t['exons'],1):line(f,t['contig'],'exon',s,e,t['strand'],'.',dict(ID=oid+'.exon'+str(i),Parent=oid))
                    for i,(s,e,p) in enumerate(t['cds'],1):line(f,t['contig'],'CDS',s,e,t['strand'],p,dict(ID=oid+'.cds'+str(i),Parent=oid))
            for ix,(tid,t) in enumerate(sorted(new.items()),1):
                oid='N'+str(ix).zfill(7);line(f,t['contig'],'gene',t['start0'],t['end0'],t['strand'],'.',dict(ID=oid,gene_name=t['gene'],gene_type='protein_coding',SwissProt=t['reference']))
                txid=oid+'.t1';line(f,t['contig'],'mRNA',t['start0'],t['end0'],t['strand'],'.',dict(ID=txid,Parent=oid,source_transcript=tid))
                for i,(s,e,p) in enumerate(t['cds'],1):line(f,t['contig'],'CDS',s,e,t['strand'],p,dict(ID=txid+'.cds'+str(i),Parent=txid))
                for kind,key in [('start_codon','starts'),('stop_codon','stops')]:
                    for s,e in t[key]:line(f,t['contig'],kind,s,e,t['strand'],'0',dict(Parent=txid))
    def tables(self):
        for r in self.rows:
            for step in STEPS[:-1]:self.verify(r,step)
        out=self.root/'tables';out.mkdir(exist_ok=True);counts=Counter(r['group'] for r in self.rows);sets={g:defaultdict(set) for g in GROUPS};before={g:defaultdict(set) for g in GROUPS};sample_rows=[];gene_types={}
        details=[]
        for r in self.rows:
            result=self.data(r,'flagger');sample=r['sample'];group=r['group']
            for key,target in [('CNG',sets),('before_CNG',before)]:
                for gene,typ in result[key].items():target[group][gene].add(sample);gene_types[gene]=typ
            newgenes=set(result['CNG'])-set(result['before_CNG'])
            n=lambda d:sum(typ=='protein_coding' and not re.match(r'^ENSG[0-9]',gene) for gene,typ in d.items())
            sample_rows.append(dict(Group=group,ID=sample,CNG_genes_before=n(result['before_CNG']),CNG_genes_after=n(result['CNG']),New_candidate_loci=len(result['retained_new'])))
            details.append(dict(sample=sample,group=group,status='COMPLETE',new_candidate_review='AVAILABLE'))
        write(out/'TableS14',sample_rows,['Group','ID','CNG_genes_before','CNG_genes_after','New_candidate_loci'])
        summary=[]
        for g in GROUPS:
            vals=[r['CNG_genes_after'] for r in sample_rows if r['Group']==g]
            if vals:summary.append(dict(Group=g,n=len(vals),mean=sum(vals)/len(vals),min=min(vals),max=max(vals),CNG_gene_union=sum(gene_types[k]=='protein_coding' and not re.match(r'^ENSG[0-9]',k) for k in sets[g])))
        write(out/'TableS14_group_summary',summary,['Group','n','mean','min','max','CNG_gene_union'])
        hist=defaultdict(Counter)
        for g in GROUPS:
            for gene,ss in sets[g].items():
                for sub,use in [('all',True),('protein',gene_types[gene]=='protein_coding'),('named_protein',gene_types[gene]=='protein_coding' and not re.match(r'^ENSG[0-9]',gene))]:
                    if use:hist[g,sub][len(ss)]+=1
        fields=['CNG_positive_haplotype_assemblies']+[g+'_'+sub for g in GROUPS for sub in ('all','protein','named_protein')]
        freqdist=[]
        for i in range(1,max(counts.values())+1):freqdist.append(dict(CNG_positive_haplotype_assemblies=i,**{g+'_'+sub:hist[g,sub][i] if i<=counts[g] else '' for g in GROUPS for sub in ('all','protein','named_protein')}))
        write(out/'TableS15',freqdist,fields)
        for ext,sep,enc in [('csv',',','utf-8-sig'),('tsv','\t','utf-8')]:
            with (out/('TableS15.two_headers.'+ext)).open('w',encoding=enc,newline='') as f:
                w=csv.writer(f,delimiter=sep);order=['HPRC-nEAS','HPRC-EAS','KPR'];w.writerow(['Group']+[v for g in order for v in [f'{g} (n={counts[g]})','']]);w.writerow(['CNG-positive haplotype assemblies (n)']+['Count_all','Count_protein']*3)
                for i in range(1,max(counts.values())+1):w.writerow([i]+[hist[g,sub][i] if i<=counts[g] else '' for g in order for sub in ('all','protein')])
        names=sorted({gene for g in GROUPS for gene in sets[g] if gene_types[gene]=='protein_coding' and not re.match(r'^ENSG[0-9]',gene)})
        freq=[];members=[];venn=Counter()
        for gene in names:
            nums={g:len(sets[g][gene]) for g in GROUPS};f={g:nums[g]/counts[g] if counts[g] else None for g in GROUPS};d1=f['KPR']-f['HPRC-EAS'] if f['KPR'] is not None and f['HPRC-EAS'] is not None else None;d2=f['KPR']-f['HPRC-nEAS'] if f['KPR'] is not None and f['HPRC-nEAS'] is not None else None
            freq.append(dict(Gene=gene,KPR_count=nums['KPR'],KPR_frequency=f['KPR'],HPRC_EAS_count=nums['HPRC-EAS'],HPRC_EAS_frequency=f['HPRC-EAS'],HPRC_nEAS_count=nums['HPRC-nEAS'],HPRC_nEAS_frequency=f['HPRC-nEAS'],KPR_minus_EAS=d1,KPR_minus_nEAS=d2,Mean_absolute_difference=(abs(d1)+abs(d2))/2 if d1 is not None and d2 is not None else None))
            gs=[g for g in GROUPS if nums[g]];cat=('Only ' if len(gs)==1 else '')+' & '.join(gs);venn[cat]+=1
            members.append(dict(Gene=gene,Category=cat,**{g.replace('-','_')+'_samples':';'.join(sorted(sets[g][gene])) for g in GROUPS}))
        ff=['Gene','KPR_count','KPR_frequency','HPRC_EAS_count','HPRC_EAS_frequency','HPRC_nEAS_count','HPRC_nEAS_frequency','KPR_minus_EAS','KPR_minus_nEAS','Mean_absolute_difference']
        freq.sort(key=lambda x:(-(x['Mean_absolute_difference'] or 0),x['Gene']));write(out/'TableS16.all_genes',freq,ff)
        selected={'PGA4','AMY1C','RASA4','UPK3BL1','USP17L11','USP17L17','GAGE2A','GAGE12E'};write(out/'TableS16.same_selected_genes',[r for r in freq if r['Gene'] in selected],ff)
        mf=['Gene','Category','KPR_samples','HPRC_EAS_samples','HPRC_nEAS_samples'];write(out/'TableS17.membership',members,mf);write(out/'TableS17.KPR_only',[r for r in members if r['Category']=='Only KPR'],mf)
        cats=[]
        for mask in range(1,8):
            gs=[g for i,g in enumerate(GROUPS) if mask&(1<<i)];cat=('Only ' if len(gs)==1 else '')+' & '.join(gs);cats.append(dict(Category=cat,Count=venn[cat]))
        write(out/'TableS17.venn_counts',cats,['Category','Count'])
        columns=[sorted(r['Gene'] for r in members if r['Category']==cat['Category']) for cat in cats]
        from itertools import zip_longest
        with (out/'TableS17.seven_columns.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.writer(f);w.writerow([c['Category'] for c in cats]);w.writerows(zip_longest(*columns,fillvalue=''))
        write(out/'sample_coverage',details,['sample','group','status','new_candidate_review'])
        self.ratios(out)
        jsave(out/'metadata.json',dict(group_denominators=dict(counts),mode=self.a.mode,before_definition='Fresh supported calls before Flagger; NOT the historical baseline snapshot',CNG_definition='Independent confirmed loci exceed GRCh38 reference gene-locus count',noncoding='Liftoff noncoding extra loci are Flagger-filtered but not protein-validated',prediction='Native mode is de novo unless explicit hints are supplied; exon/intron support is recomputed from alignment intervals, not imported as historical P-hint annotations',historical_reproduction='For exact historical S14-S17 use the original completed-output integration package. This portable implementation performs a fresh analysis.'))
        (self.root/'completed.ok').write_text('All requested samples and tables completed.\n');log('[DONE] '+str(self.root))
    def ratios(self,out):
        if not self.a.chm13_annotation:return
        cg,ct=gff(self.a.chm13_annotation);den_ids=defaultdict(set)
        for gid,g in cg.items():
            ref=canonical(gid,g['copy'])
            if ref in self.reference_genes:
                typ=self.reference_genes[ref]['biotype'];den_ids['protein_coding_genes' if typ=='protein_coding' else 'noncoding_genes'].add(ref)
        for tid,t in ct.items():
            ref=canonical(tid,cg[t['gene_id']]['copy'])
            if ref in self.reference_tx:
                typ=self.reference_tx[ref]['biotype'];den_ids['protein_coding_transcripts' if typ=='protein_coding' else 'noncoding_transcripts'].add(ref)
        den=Counter({k:len(v) for k,v in den_ids.items()});need(sum(den.values())>0,'CHM13 annotation has no matching reference GENCODE IDs')
        rows=[]
        for r in self.rows:
            lif=self.data(r,'liftoff');ids=defaultdict(set)
            for gid,g in lif['genes'].items():
                ref=canonical(gid,g['copy'])
                if ref in self.reference_genes:ids['protein_coding_genes' if g['biotype']=='protein_coding' else 'noncoding_genes'].add(ref)
            for tid,t in lif['transcripts'].items():
                g=lif['genes'][t['gene_id']];ref=canonical(tid,g['copy'])
                if ref in self.reference_tx:
                    typ=self.reference_tx[ref]['biotype'];ids['protein_coding_transcripts' if typ=='protein_coding' else 'noncoding_transcripts'].add(ref)
            for category in ('protein_coding_genes','noncoding_genes','protein_coding_transcripts','noncoding_transcripts'):
                rows.append(dict(sample=r['sample'],group=r['group'],category=category,unique_GENCODE_IDs=len(ids[category]),CHM13_count=den[category],annotation_ratio=len(ids[category])/den[category] if den[category] else ''))
        write(out/'TableS13.annotation_ratios',rows,['sample','group','category','unique_GENCODE_IDs','CHM13_count','annotation_ratio'])

class Bed:
    def __init__(self,path,mapfile,sample,shared_bed=False):
        self.shared_bed=shared_bed;self.sample=sample
        match=re.search(r'(?:\.|_hap|_)([12])$',sample);self.haplotype=match.group(1) if match else None
        raw=defaultdict(list);self.mapping={};self.index={}
        if mapfile:
            for r in table(mapfile):
                if r['sample']!=sample:continue
                need(r['assembly_contig'] not in self.mapping,'Duplicate explicit contig mapping')
                self.mapping[r['assembly_contig']]=(r['bed_contig'],int(r.get('offset0','0') or 0))
        op=gzip.open if str(path).endswith('.gz') else open
        with op(path,'rt') as f:
            for n,line in enumerate(f,1):
                if not line.strip() or line.startswith(('#','track','browser')):continue
                a=line.rstrip('\n').split('\t');need(len(a)>=4,f'{path}:{n}: classified BED4 required')
                s,e=int(a[1]),int(a[2]);lab=LABELS.get(a[3].strip().lower());need(0<=s<e and lab,f'{path}:{n}: unknown label/coordinates');raw[a[0]].append((s,e,lab))
        need(raw and any(lab=='Hap' for vals in raw.values() for _,_,lab in vals),'Full classified BED with Hap intervals required')
        for seq,vals in raw.items():
            events=defaultdict(Counter)
            for s,e,lab in vals:events[s][lab]+=1;events[e][lab]-=1
            active=Counter();prev=None;segments=[]
            for pos,delta in sorted(events.items()):
                if prev is not None and pos>prev and active:segments.append((prev,pos,tuple(sorted(active))))
                active.update(delta)
                for k in list(active):
                    if active[k]==0:del active[k]
                prev=pos
            self.index[seq]=(segments,[x[0] for x in segments])
    def assess(self,seq,s,e,fa):
        target,offset=self.mapping.get(seq,(seq,0));parts=target.split('#');source=seq.split('#')
        if len(parts)>=3 and len(source)>=3:need(parts[0:2]==source[0:2],'Explicit mapping attempts to cross sample/haplotype')
        if len(parts)>=3 and self.haplotype is not None:need(parts[1]==self.haplotype,'BED mapping belongs to the other haplotype')
        cov=Counter();need(seq in fa.idx and 0<=s<e<=fa.idx[seq][0],'Gene coordinates outside assembly')
        if self.shared_bed and len(parts)<3:return 'UNASSESSED','UNTAGGED_BED_SHARED_ACROSS_ASSEMBLIES',cov
        entry=self.index.get(target)
        if not entry:return 'UNASSESSED','CONTIG_NOT_MATCHED',cov
        fs,fe=s+offset,e+offset
        if fs<0:return 'UNASSESSED','NEGATIVE_PROJECTED_COORDINATE',cov
        segs,starts=entry;i=max(0,bisect.bisect_right(starts,fs)-1)
        while i<len(segs) and segs[i][0]<fe:
            a,b,labels=segs[i];bp=max(0,min(fe,b)-max(fs,a));cov['Covered']+=bp
            for lab in labels:cov[lab]+=bp
            if any(x in labels for x in ('Dup','Col','Err')):cov['Issue']+=bp
            i+=1
        if cov['Issue']:return 'EXCLUDE','OVERLAP_DUP_COL_ERR',cov
        if cov['Covered']!=e-s:return 'UNASSESSED','INCOMPLETE_BED_COVERAGE',cov
        if cov['Unk'] or cov['Hap']!=e-s:return 'UNASSESSED','NOT_FULLY_HAP',cov
        return 'KEEP','FULLY_HAP',cov

def cng(genes,new,refcounts,reference_genes):
    by=defaultdict(list);types={};newcount=Counter()
    for g in genes:
        if g['gene']:
            by[g['gene']].append(g);types[g['gene']]=g['biotype']
    for t in new:by[t['gene']].append(t);types[t['gene']]='protein_coding';newcount[t['gene']]+=1
    result={};review=[]
    for gene,items in sorted(by.items()):
        keys=defaultdict(list)
        for g in items:keys[g['contig'],g['strand']].append([g['start0'],g['end0']])
        n=0
        for parts in keys.values():
            loci=[]
            for a,b in sorted(parts):
                if loci and a<loci[-1][1]:loci[-1][1]=max(loci[-1][1],b)
                else:loci.append([a,b])
            n+=len(loci)
        refn=len(refcounts.get(gene,set()));positive=refn>0 and n>refn
        if positive:result[gene]=types[gene]
        review.append(dict(gene=gene,reference_loci=refn,confirmed_loci=n,CNG=int(positive),new_prediction_loci=newcount[gene]))
    return result,review

def parser():
    p=argparse.ArgumentParser(prog='2.Ensemble_gene_annotation.sh',description='Seven-step annotation and CNG workflow; full native execution or precomputed-output replay.')
    p.add_argument('step',choices=['all','check']+STEPS,nargs='?',default='check')
    p.add_argument('--samples',type=Path,required=True,help='TSV manifest; paths relative to manifest directory')
    p.add_argument('--reference-gtf',type=Path,required=True)
    p.add_argument('--proteins',type=Path,required=True,help='Human Swiss-Prot FASTA with GN= symbols')
    p.add_argument('--reference-fasta',type=Path)
    p.add_argument('--chm13-annotation',type=Path,help='Optional annotation-ratio denominator')
    p.add_argument('--protein-gene-map',type=Path,help='Optional TSV: protein_id, gene')
    p.add_argument('--satellite-library',type=Path,help='Optional second RepeatMasker library')
    p.add_argument('--out',type=Path,default=Path('gene_annotation_results'))
    p.add_argument('--mode',choices=['tools','outputs'],default='tools')
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--species',default='human')
    p.add_argument('--exonerate-batch-bp',type=int,default=4000000)
    p.add_argument('--protein-identity',type=float,default=95,help='AUGUSTUS/Swiss-Prot protein identity (%)')
    p.add_argument('--blastp-identity',type=float,default=90,help='BLASTp identity (%) for unannotated loci')
    p.add_argument('--exonerate-identity',type=float,default=90)
    p.add_argument('--support-coverage',type=float,default=80)
    p.add_argument('--liftoff-overlap',type=float,default=80)
    p.add_argument('--min-exon-support',type=float,default=80,help='Percentage of coding exons with at least 80 percent alignment coverage')
    p.add_argument('--min-intron-support',type=float,default=80,help='Fraction of AUGUSTUS introns supported by Exonerate')
    p.add_argument('--intron-match',type=float,default=80)
    p.add_argument('--unannotated-overlap',type=float,default=95)
    return p

def main():
    a=parser().parse_args();need(a.threads>0 and a.exonerate_batch_bp>0,'Positive threads and batch size required')
    for key in ('protein_identity','blastp_identity','exonerate_identity','support_coverage','liftoff_overlap','intron_match','unannotated_overlap','min_exon_support','min_intron_support'):need(0<getattr(a,key)<=100,'Percentage out of range: '+key)
    for key in ('reference_gtf','proteins','reference_fasta','chm13_annotation','protein_gene_map','satellite_library'):
        val=getattr(a,key)
        if val:need(val.is_file(),f'Input missing: {val}');setattr(a,key,val.resolve())
    rows=inputs(a)
    if a.mode=='tools':
        needed=set()
        for r in rows:
            if not r['masked_fasta'] and a.step in ('mask','all','check'):needed.add('RepeatMasker')
            if not r['augustus_gff'] and a.step in ('augustus','all','check'):needed.add('augustus')
            if not r['exonerate_gff'] and a.step in ('exonerate','all','check'):needed.add('exonerate')
            if not r['liftoff_gff'] and a.step in ('liftoff','all','check'):needed.add('liftoff');need(a.reference_fasta,'--reference-fasta required for native Liftoff')
            if not r['blastp_tsv'] and a.step in ('blastp','all','check'):needed.update(['blastp','makeblastdb'])
        for exe in sorted(needed):need(shutil.which(exe),f'Executable missing from PATH: {exe}')
    try:import Bio
    except ImportError:raise ValueError('Biopython required: python3 -m pip install biopython')
    if a.step=='check':log('[CHECK OK] '+str(len(rows))+' haplotypes; mode='+a.mode+'; Biopython '+Bio.__version__);return
    p=Pipeline(a)
    for step in STEPS if a.step=='all' else [a.step]:p.run(step)

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,KeyError,csv.Error) as exc:raise SystemExit('[ERROR] '+str(exc))
PY
}


usage() {
  cat <<'EOF'
Usage: bash 2.Ensemble_gene_annotation.sh <step> [options]

  softmask             RepeatMasker soft masking                     (alias: mask)
  predict_orfs         AUGUSTUS coding transcript prediction         (alias: augustus)
  protein_support      Exonerate alignment and candidate validation  (alias: exonerate)
  gene_annotation      Liftoff annotation and extra-copy detection   (alias: liftoff)
  compare_annotations  CDS/intron concordance with Liftoff           (alias: compare)
  search_unannotated   BLASTp review of unannotated candidates       (alias: blastp)
  filter_flagger       assembly-error filtering and final annotation (alias: flagger)
  tables               Table S13 (optional) and Tables S14-S17
  check                input and software checks
  all                  run steps 1-7 and write tables

Configuration variables:
  REF_DIR              reference directory (default: ./reference)
  LIST_DIR             sample-table directory (default: ./lists)
  SAMPLE_TABLE         sample TSV (default: $LIST_DIR/samples.tsv)
  GENCODE_GTF          GRCh38 GENCODE v46 annotation
  SWISSPROT_FA         human Swiss-Prot protein FASTA
  REFERENCE_FA         GRCh38 FASTA for native Liftoff
  CHM13_ANNOTATION     optional GENCODE-based CHM13 annotation for S13
  OUT_DIR              results directory (default: ./results/gene_annotation)
  THREADS              native parallelism setting (default: 8)
  RUN_MODE             tools or outputs (default: tools)
  MIN_PROTEIN_IDENTITY  AUGUSTUS/Swiss-Prot 95% identity
  MIN_BLASTP_IDENTITY   BLASTp 90% identity

Command-line overrides:
  --samples FILE             --reference-gtf FILE      --proteins FILE
  --reference-fasta FILE     --chm13-annotation FILE   --out DIR
  --mode tools|outputs       --threads N               --species NAME
  --satellite-library FILE   --protein-gene-map FILE
  --protein-identity PCT     --blastp-identity PCT
  --exonerate-identity PCT
  --support-coverage PCT     --liftoff-overlap PCT
  --min-exon-support PCT     --min-intron-support PCT
  --intron-match PCT         --unannotated-overlap PCT
  --exonerate-batch-bp N

Examples:
  export REF_DIR=/path/to/reference
  export SAMPLE_TABLE=/path/to/samples.tsv
  export OUT_DIR=/path/to/results
  export REFERENCE_FA=/path/to/GRCh38.fa
  bash 2.Ensemble_gene_annotation.sh check
  bash 2.Ensemble_gene_annotation.sh all

See src/Ensemble_gene_annotation.md for sample manifests and workflow details.
EOF
}

main() {
  local step
  case "${1:-}" in
    ''|-h|--help) usage; return 0 ;;
    --*) step=check ;;
    *) step=$1; shift ;;
  esac

  parse_options "$@"
  case "${step}" in
    softmask|mask|predict_orfs|augustus|protein_support|exonerate|\
    gene_annotation|liftoff|compare_annotations|compare|search_unannotated|blastp|\
    filter_flagger|flagger|tables|check|all) ;;
    *) usage; die "unknown step: ${step}" ;;
  esac
  configure_runtime

  case "${step}" in
    softmask|mask)                 softmask ;;
    predict_orfs|augustus)         predict_orfs ;;
    protein_support|exonerate)     protein_support ;;
    gene_annotation|liftoff)       gene_annotation ;;
    compare_annotations|compare)  compare_annotations ;;
    search_unannotated|blastp)     search_unannotated ;;
    filter_flagger|flagger)        filter_flagger ;;
    tables)                       make_tables ;;
    check)                        check_inputs ;;
    all)                          run_all ;;
  esac
}

main "$@"
