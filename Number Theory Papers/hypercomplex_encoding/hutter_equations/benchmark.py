#!/usr/bin/env python3
"""Discover direct arithmetic byte constraints, test on chronological enwik8 holdouts.

The equations predict byte values, not hashed context IDs. Each relation produces
an exactly invertible (observed-predicted) mod-256 correction. Report two separate
measurements: ideal prequential conditional code length and real standalone Brotli
of the residual streams (with seed, model metadata and CRC). Hyperparameters and
model choices are selected on validation only. No claims of Hutter competitiveness.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha1, sha256
import json
import math
from pathlib import Path
import time
import numpy as np
import brotli
import zlib
from equations import NAMES,FAMILIES,SEED_LENGTH,prev_class,predictions,encode,decode

MIX=(0., .025,.05,.10,.20,.40,.60,.80,1.)
WGRID=(0.,.025,.05,.10,.20,.40)
PRI=0.25
SPLIT=(200_000,60_000,120_000)


def baseline_probabilities(data: bytes) -> np.ndarray:
    """Causal smoothed byte n-gram of orders 1..4 with globally adaptive prior."""
    tables=[defaultdict(Counter) for _ in range(4)]
    char_counts=np.zeros(256,dtype=np.int64)
    p=np.zeros(len(data),dtype=np.float64)
    for i,x in enumerate(data):
        pb=(char_counts[x]+0.5)/(i+128.0)
        for k,tau in ((1,12.),(2,8.),(3,6.),(4,5.)):
            if i>=k:
                ctx=data[i-k:i]
                counts=tables[k-1].get(ctx)
                if counts:
                    pb=(counts[x]+tau*pb)/(counts.total()+tau)
        p[i]=pb
        char_counts[x]+=1
        for k in range(1,5):
            if i>=k:
                tables[k-1][data[i-k:i]][x]+=1
    return p


def prior_classes(data: bytes) -> np.ndarray:
    ar=np.frombuffer(data,dtype=np.uint8)
    preceding=np.empty(len(ar),dtype=np.uint8)
    preceding[0]=0
    preceding[1:]=ar[:-1]
    return np.fromiter((prev_class(int(v)) for v in preceding),count=len(data),dtype=np.uint8)


def causal_residual_probabilities(data: np.ndarray, pred: np.ndarray, classes:np.ndarray) -> np.ndarray:
    """Conditional probability of actual residual, updated after each observed byte.

    For every equation f and class c: q(b | f,c) = smoothed count of residue
    ((b - f(history)) mod 256). Mapping b -> residue is bijective for fixed f.
    """
    residual=((data.astype(np.int16)-pred.astype(np.int16))&255).astype(np.uint8)
    counts=np.zeros((8,256),dtype=np.int32)
    totals=np.zeros(8,dtype=np.int32)
    result=np.empty(len(data),dtype=np.float64)
    for i in range(len(data)):
        cl=int(classes[i]);r=int(residual[i])
        result[i]=(counts[cl,r]+PRI)/(totals[cl]+256*PRI)
        counts[cl,r]+=1
        totals[cl]+=1
    return result


def combo_fit(p0, q1, q2, va, ve, te, stop):
    valmin=float('inf');weights=None
    for a in WGRID:
        for b in WGRID:
            if a+b>0.8+1e-9:continue
            vv=((1-a-b)*p0[va:ve]+a*q1[va:ve]+b*q2[va:ve])
            loss=float(-np.log2(vv).sum())
            if loss<valmin:valmin=loss;weights=(a,b)
    a,b=weights
    test=float(-np.log2((1-a-b)*p0[te:stop]+a*q1[te:stop]+b*q2[te:stop]).sum())
    return {'weights':list(weights),'validation_bits':valmin,'heldout_bits':test}


def run(data:bytes, warmup=200_000, validation=60_000, heldout=120_000, source_offset=0, output=None, source_label="provided byte input (not independently authenticated)"):
    if len(data)<warmup+validation+heldout:
        raise ValueError('input does not cover declared split')
    t0=time.perf_counter()
    stop=warmup+validation+heldout
    data=data[:stop]
    b=np.frombuffer(data,dtype=np.uint8)
    classes=prior_classes(data)
    p0=baseline_probabilities(data)
    pred=predictions(data)
    va,ve,te=warmup,warmup+validation,warmup+validation
    base_valid=float(-np.log2(p0[va:ve]).sum())
    base_test=float(-np.log2(p0[te:]).sum())
    model_rows=[]; qdict={}
    for j,(name,family) in enumerate(zip(NAMES,FAMILIES)):
        q=causal_residual_probabilities(b,pred[name],classes)
        qdict[name]=q
        scores=np.array([-np.log2((1-w)*p0[va:ve]+w*q[va:ve]).sum() for w in MIX])
        idx=int(np.argmin(scores));w=MIX[idx]
        test_loss=float(-np.log2((1-w)*p0[te:]+w*q[te:]).sum())
        exact=float(np.mean(b[te:]==pred[name][te:]))
        row={'name':name,'family':family,'mix_frozen':w,
             'validation_bits':float(scores[idx]),'heldout_bits':test_loss,
             'heldout_gain_bytes':(base_test-test_loss)/8,
             'test_exact_equation_pct':round(100*exact,4),
             'validation_gain_bytes':(base_valid-scores[idx])/8}
        model_rows.append(row)
    select={}
    for fam in sorted(set(FAMILIES)):
        select[fam]=min((r for r in model_rows if r['family']==fam),key=lambda r:r['validation_bits'])
    ordinary=[r for r in model_rows if r['family'] not in ('complex','quaternion')]
    topordinary=min(ordinary,key=lambda r:r['validation_bits'])
    best_nonquat=min((r for r in model_rows if r['family']!='quaternion'),key=lambda r:r['validation_bits'])
    best_quat=select['quaternion']
    extra={}
    for key,a,bm in (
        ('ordinary_plus_quaternion',topordinary['name'],best_quat['name']),
        ('ordinary_plus_complex',topordinary['name'],select['complex']['name']),
        ('ordinary_plus_affine',topordinary['name'],select['affine']['name']),
        ('ordinary_plus_xor',topordinary['name'],select['xor']['name']),
    ):
        combo=combo_fit(p0,qdict[a],qdict[bm],va,ve,te,stop)
        combo['equations']=[a,bm]
        combo['heldout_gain_bytes']=(base_test-combo['heldout_bits'])/8
        extra[key]=combo
    # Discover different equations for each byte class *only using validation*.
    class_map=[];class_support=[]
    for cl in range(8):
        where=np.flatnonzero(classes[va:ve]==cl)+va
        class_support.append(int(len(where)))
        if len(where)<200:
            winner=topordinary['name']
        else:
            winner=min(NAMES,key=lambda name:float(-np.log2(qdict[name][where]).sum()))
        class_map.append(NAMES.index(winner))
    suffix=data[te:]
    raw_baselines={'brotli11':len(brotli.compress(suffix,quality=11)),
                   'brotli9':len(brotli.compress(suffix,quality=9)),
                   'zlib9':len(zlib.compress(suffix,9))}
    framed={}
    for fam,row in select.items():
        mid=NAMES.index(row['name'])
        frame=encode(suffix,mid,quality=11)
        if decode(frame)!=suffix:raise AssertionError(f'nonexact {row["name"]}')
        framed[fam]={'equation':row['name'],'bytes':len(frame),
                     'vs_brotli11_bytes':len(frame)-raw_baselines['brotli11']}
    classframe=encode(suffix, 0, bytes(class_map),quality=11)
    if decode(classframe)!=suffix: raise AssertionError('class model did not roundtrip')
    framed['class_selected']={'equations':[NAMES[x] for x in class_map],
             'validation_class_support':class_support,'bytes':len(classframe),
             'vs_brotli11_bytes':len(classframe)-raw_baselines['brotli11']}
    result={'source':source_label,
            'offset':source_offset, 'window_sha256':sha256(data).hexdigest(),
            'window_bytes':len(data),'warmup_bytes':warmup,'validation_bytes':validation,
            'heldout_bytes':heldout,'heldout_sha256':sha256(suffix).hexdigest(),
            'baseline_ideal_bytes':base_test/8,'validation_baseline_ideal_bytes':base_valid/8,
            'baseline_actual_bytes':raw_baselines,
            'model_rows':model_rows,'family_winners':select,
            'best_ordinary_equation':topordinary['name'],
            'best_nonquaternion_equation':best_nonquat['name'],
            'joint_ablations':extra,'actual_framed_bytes':framed,
            'elapsed_sec':round(time.perf_counter()-t0,2),
            'caveats':'Prequential ideal code length excludes arithmetic coder and full decoder source; framed Brotli residuals include 14-byte model header, 16 seed bytes and 8-byte class map if used. Only the equation class is tuned on validation; online model may update on decoded test bytes. This tests equations on natural Wiki text, not an enwik9 prize score.'}
    if output is not None:
        write_result(result,Path(output))
    return result


def write_result(out,path):
    path.mkdir(parents=True,exist_ok=True)
    (path/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    with (path/'equations.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(out['model_rows'][0]))
        writer.writeheader();writer.writerows(out['model_rows'])
    a=out['actual_framed_bytes'];lines=[
        '# Direct algebraic-equation discovery on real Wikipedia bytes',
        '',f'- Source: {out["source"]}; offset: {out["offset"]:,}; window SHA-256: `{out["window_sha256"]}`',
        f'- Chronological split: {out["warmup_bytes"]:,} warmup / {out["validation_bytes"]:,} validation / {out["heldout_bytes"]:,} heldout.',
        f'- Held-out SHA-256: `{out["heldout_sha256"]}`',
        f'- Causal order-4 byte-context baseline ideal length: **{out["baseline_ideal_bytes"]:,.1f} bytes**.',
        f'- Direct standalone Brotli-11: **{out["baseline_actual_bytes"]["brotli11"]:,} bytes**.',
        '', '## Exact algebraic constraints on the next byte',
        '', 'Each candidate equation directly predicts b_t from previous bytes, and transmits the exact modular residual (b_t - predicted) mod 256. A known context category conditions the residual probability. Validation selects a fixed blend with the same byte-n-gram baseline. No equation is assumed to hold universally.',
        '', '| Family | Selected equation | Equation exact on heldout | Ideal bytes | Gain vs ngram | Standalone framed Brotli bytes |',
        '|---|---|---:|---:|---:|---:|',]
    for fam,row in sorted(out['family_winners'].items()):
        entry=a[fam];lines.append(f'| {fam} | `{row["name"]}` | {row["test_exact_equation_pct"]:.2f}% | {row["heldout_bits"]/8:,.1f} | {row["heldout_gain_bytes"]:+,.1f} | {entry["bytes"]:,} |')
    lines.extend(['','### Multiple-equation ablations', '', '| Combination | Validation-frozen weights | Test ideal bytes | Gain vs ngram |', '|---|---|---:|---:|'])
    for name, row in out['joint_ablations'].items():
        lines.append(f'| {name} (`{row["equations"][0]}` + `{row["equations"][1]}`) | {row["weights"]} | {row["heldout_bits"]/8:,.1f} | {row["heldout_gain_bytes"]:+,.1f} |')
    classes=a['class_selected'];lines.extend(['',f'**Validation-discovered class-conditional equation map:** {classes["bytes"]:,} bytes in a standalone stream (vs direct Brotli-11 {out["baseline_actual_bytes"]["brotli11"]:,} bytes).',
                 '', 'Class-selected equations in order lowercase, uppercase, digit, whitespace, XML markup, control, non-ASCII, other:',
                 '',f'`{classes["equations"]}`',
                 '', '## Interpretation and scope',
                 '', '- Formulas are explicit modular arithmetic, complex products and true Hamilton quaternion products. They predict *the next byte itself*, not a hash of preceding bytes.',
                 '- Residuals are exactly invertible and every standalone stream is round-trip verified, including its checksum.',
                 '- Ideal bits are a diagnostic model comparison, **not** compressed file bytes. Framed residual streams are actual compressed files; direct Brotli is compared on the same heldout bytes.',
                 '- Training/validation and holdout are chronological. No heldout bytes select formulas, feature weights or class map.',
                 '- Fixed formula-pool search and coding overhead must be charged before claiming competitiveness. Other character-based compressors and dictionary methods may capture the same dependencies more efficiently.',
                 f'- Wall-clock experiment: {out["elapsed_sec"]:.1f} seconds.', ''])
    (path/'REPORT.md').write_text('\n'.join(lines))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--offsets',nargs='+',type=int,default=[0,2_000_000])
    ap.add_argument('--warmup',type=int,default=200_000)
    ap.add_argument('--validation',type=int,default=60_000)
    ap.add_argument('--heldout',type=int,default=120_000)
    ap.add_argument('--output',default='results')
    ap.add_argument('--source-label',default='provided byte input (not independently authenticated)')
    args=ap.parse_args()
    src=Path(args.input)
    outputs=[]
    for offset in args.offsets:
        with src.open('rb') as f:
            f.seek(offset)
            data=f.read(args.warmup+args.validation+args.heldout)
        result=run(data,args.warmup,args.validation,args.heldout,offset,Path(args.output)/f'offset_{offset}',args.source_label)
        outputs.append(result)
        print(f'OFFSET {offset}: direct Brotli-11={result["baseline_actual_bytes"]["brotli11"]} ' +
              f'class-map={result["actual_framed_bytes"]["class_selected"]["bytes"]} ' +
              f'best ordinary={result["best_ordinary_equation"]} ' +
              f'quaternion={result["family_winners"]["quaternion"]["name"]}',flush=True)
        print((Path(args.output)/f'offset_{offset}'/'REPORT.md').read_text(),flush=True)
    summary={ 'runs':[{'offset':r['offset'],'sha256':r['window_sha256'],
                         'standalone_brotli11':r['baseline_actual_bytes']['brotli11'],
                         'standalone_class_rules':r['actual_framed_bytes']['class_selected']['bytes'],
                         'best_ordinary':r['best_ordinary_equation'],
                         'quaternion':r['family_winners']['quaternion']['name'],
                         'ordinary_ideal_gain':next(row for row in r['model_rows'] if row['name']==r['best_ordinary_equation'])['heldout_gain_bytes'],
                         'quaternion_ideal_gain':r['family_winners']['quaternion']['heldout_gain_bytes']} for r in outputs],
               'corpus_full_sha1_required':'57b8363b814821dc9d47aa4d41f58733519076b2'}
    folder=Path(args.output);folder.mkdir(parents=True,exist_ok=True)
    (folder/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('RUN SUMMARY',json.dumps(summary),flush=True)

if __name__=='__main__':main()
