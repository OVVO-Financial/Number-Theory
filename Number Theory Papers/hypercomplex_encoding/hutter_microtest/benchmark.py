#!/usr/bin/env python3
"""Causal enwik8 prediction microtest: Hamilton quaternion vs matched controls."""
import argparse
from collections import Counter,defaultdict
from hashlib import sha256
import json,math,time,zlib
from pathlib import Path
import brotli

MIX=(0,0.05,0.1,0.2,0.4,0.8,1.0)
FEATURES=('hamilton','xor_pairs','add_pairs','hash8')

def hamilton(a,b):
    w,x,y,z=a
    r,s,t,u=b
    return (w*r-x*s-y*t-z*u,w*s+x*r+y*u-z*t,
            w*t-x*u+y*r+z*s,w*u+x*t-y*s+z*r)

def feature_ids(history):
    h=history[-8:]
    a=[int(x)-128 for x in h[:4]]
    b=[int(x)-128 for x in h[4:]]
    q=hamilton(a,b)
    hk=(((q[0]>>10)&7)<<6)|(((q[1]>>10)&7)<<3)|((q[2]>>10)&7)
    xor=(((h[0]^h[4])>>5)<<6)|(((h[1]^h[5])>>5)<<3)|((h[2]^h[6])>>5)
    add=(((h[0]+h[4])>>6)<<6)|(((h[1]+h[5])>>6)<<3)|((h[2]+h[6])>>6)
    hv=2166136261
    for bval in h: hv=((hv^int(bval))*16777619)&0xffffffff
    return hk,xor,add,hv&511

def run(data,warmup,validation,heldout,output=None):
    assert len(data)>=warmup+validation+heldout
    limit=warmup+validation+heldout
    data=data[:limit]
    ngrams=[defaultdict(Counter) for _ in range(4)]
    tables=[{} for _ in FEATURES]
    totals=[{} for _ in FEATURES]
    chars=Counter()
    n=0
    losses={s:{f:[0.0]*len(MIX) for f in FEATURES} for s in ('validation','heldout')}
    baselines={'validation':0.0,'heldout':0.0}
    elapsed_start=time.perf_counter()
    for i in range(limit):
        cur=int(data[i]); hist=data[max(0,i-8):i]
        region=('validation' if warmup<=i<warmup+validation else
                'heldout' if i>=warmup+validation else None)
        if region:
            p=(chars[cur]+0.5)/(n+128)
            for order,tau in ((1,12.0),(2,8.0),(3,6.0),(4,5.0)):
                if i>=order:
                    bucket=ngrams[order-1].get(bytes(data[i-order:i]))
                    if bucket:
                        p=(bucket[cur]+tau*p)/(bucket.total()+tau)
            baselines[region]-=math.log2(p)
            for name,key,table,total in zip(FEATURES,feature_ids(hist),tables,totals):
                bucket=table.get(key)
                count=bucket[cur] if bucket else 0
                qty=total.get(key,0)
                pfeature=(count+48.0*p)/(qty+48.0)
                for j,mix in enumerate(MIX):
                    losses[region][name][j]-=math.log2((1-mix)*p+mix*pfeature)
        chars[cur]+=1
        n+=1
        for order in range(1,5):
            if i>=order:
                ngrams[order-1][bytes(data[i-order:i])][cur]+=1
        if i>=8:
            for name,key,table,total in zip(FEATURES,feature_ids(hist),tables,totals):
                if key not in table: table[key]=Counter()
                table[key][cur]+=1
                total[key]=total.get(key,0)+1
    out={'sample_bytes':limit,'sample_sha256':sha256(data).hexdigest(),
         'warmup_bytes':warmup,'validation_bytes':validation,'heldout_bytes':heldout,
         'baseline_ideal_bytes':baselines['heldout']/8,
         'baseline_bpb':baselines['heldout']/heldout,
         'direct_brotli11_bytes':len(brotli.compress(data[-heldout:],quality=11)),
         'direct_zlib9_bytes':len(zlib.compress(data[-heldout:],9)),
         'time_seconds':time.perf_counter()-elapsed_start,
         'methods':[]}
    for name in FEATURES:
        chosen=min(range(len(MIX)),key=lambda j:losses['validation'][name][j])
        bits=losses['heldout'][name][chosen]
        out['methods'].append({
           'name':name,'validation_frozen_weight':MIX[chosen],
           'validation_gain_ideal_bytes':(baselines['validation']-losses['validation'][name][chosen])/8,
           'heldout_ideal_bytes':bits/8,
           'heldout_bpb':bits/heldout,
           'heldout_gain_vs_baseline_bytes':(baselines['heldout']-bits)/8})
    if output:
        path=Path(output)
        path.mkdir(parents=True,exist_ok=True)
        (path/'result.json').write_text(json.dumps(out,indent=2))
        lines=['# Hypercomplex context ablation on official enwik8','',
             'Measured ideal adaptive coding lengths; not standalone archives.',
             f'Corpus prefix: {limit:,} bytes, SHA-256 {out["sample_sha256"]}.',
             f'Split: {warmup:,} warmup / {validation:,} validation / {heldout:,} heldout.',
             '','| Predictor | Frozen weight | Held-out ideal bytes | Gain vs ngram |',
             '|---|---:|---:|---:|',
             f'| Adaptive 4-gram | - | {out["baseline_ideal_bytes"]:,.1f} | - |']
        for m in out['methods']:
            lines.append(f'| {m["name"]} | {m["validation_frozen_weight"]} | {m["heldout_ideal_bytes"]:,.1f} | {m["heldout_gain_vs_baseline_bytes"]:+,.1f} |')
        lines += ['',
             f'Standalone Brotli-11 on same held-out bytes: {out["direct_brotli11_bytes"]:,} bytes.',
             f'Standalone zlib-9: {out["direct_zlib9_bytes"]:,} bytes.',
             '',
             'All models are causal; weights frozen using validation. '
             'Each comparator has the same 512-bin key budget. '
             'Model implementation/entropy-coder bytes and warmup-prefix cost excluded. '
             'Not a Hutter Prize score.']
        (path/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True)
    parser.add_argument('--warmup',type=int,default=200000)
    parser.add_argument('--validation',type=int,default=60000)
    parser.add_argument('--heldout',type=int,default=120000)
    parser.add_argument('--output',default='results')
    args=parser.parse_args()
    with Path(args.input).open('rb') as handle:
        data=handle.read(args.warmup+args.validation+args.heldout)
    run(data,args.warmup,args.validation,args.heldout,args.output)
    print(Path(args.output,'REPORT.md').read_text(),flush=True)

if __name__=='__main__': main()
