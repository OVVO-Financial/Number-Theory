#!/usr/bin/env python3
"""Causal, lossless predictive encoding vs matched nonlinear controls.

The model is trained ONLY on a transmitted prefix (at most 25%, capped at 256 frames). All models
predict each block from strictly previous blocks. Full model and frame bytes
are counted, and the decoder reconstructs sequentially from actual residues.

Benchmarks are exploratory diagnostics, not evidence of SOTA competitiveness.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import lzma
import os
import platform
import statistics
import time
from pathlib import Path
import numpy as np
import brotli
from sota_benchmark.codec import Zstd
from predictive_lab import predictor as pr

PROJECT=Path(__file__).resolve().parents[1]
PAPER=PROJECT.parent/'Real_Pair_Representation_of_Complex_Numbers_and_Quaternions.pdf'
FILES=[
    ('algebra_source',PROJECT/'hypercomplex_highdim_test.py','heldout_natural'),
    ('cascade_source',PROJECT/'hypercomplex_cascade_test.py','heldout_natural'),
    ('sota_codec_source',PROJECT/'sota_benchmark/codec.py','heldout_natural'),
    ('sota_benchmark_source',PROJECT/'sota_benchmark/benchmark.py','heldout_natural'),
    ('paper_text',PROJECT/'realpair_text.txt','heldout_natural_derived'),
    ('paper_pdf',PAPER,'heldout_natural'),
]


def generate_cd(frames=1024,D=16,seed=20261118)->bytes:
    rng=np.random.default_rng(seed)
    x=np.zeros((frames,D),dtype=np.uint8)
    x[0]=rng.integers(90,167,size=D).astype(np.uint8)
    x[1]=rng.integers(90,167,size=D).astype(np.uint8)
    for t in range(2,frames):
        p1=x[t-1].astype(np.int64)-128
        p2=x[t-2].astype(np.int64)-128
        product=pr.pair_feature(p1,p2,1)
        values=128+0.50*p1+0.04*p2+3.0*product+rng.normal(0,14,size=D)
        x[t]=np.clip(np.rint(values),0,255).astype(np.uint8)
    return x.tobytes()


def generate_ar(frames=1024,D=16,seed=20261119)->bytes:
    rng=np.random.default_rng(seed)
    x=np.zeros((frames,D),dtype=np.uint8)
    x[0]=rng.integers(80,176,size=D).astype(np.uint8)
    x[1]=rng.integers(80,176,size=D).astype(np.uint8)
    for t in range(2,frames):
        a=x[t-1].astype(np.int64)-128
        b=x[t-2].astype(np.int64)-128
        x[t]=np.clip(np.rint(128+0.75*a-0.18*b+rng.normal(0,8,size=D)),0,255).astype(np.uint8)
    return x.tobytes()


def corpus(smoke=False,natural_only=False)->list[dict]:
    sources=[]
    for name,path,label in FILES:
        if not path.exists():
            raise FileNotFoundError(f'missing independently tracked input: {path}')
        raw=path.read_bytes()
        if smoke:raw=raw[:min(len(raw),4096)]
        sources.append({'id':name,'bytes':len(raw),'label':label,'path':str(path.relative_to(PROJECT)) if path.is_relative_to(PROJECT) else '../'+path.name,
                        'sha256':hashlib.sha256(raw).hexdigest(),'raw':raw,'truncated':smoke})
    if not natural_only:
        synthetics={'cd_coupled_control':generate_cd(frames=200 if smoke else 1024),
                    'linear_ar_control':generate_ar(frames=200 if smoke else 1024),
                    'random_control':np.random.default_rng(20261120).integers(0,256,size=4096 if smoke else 16384,dtype=np.uint8).tobytes()}
        for name,raw in synthetics.items():
            sources.append({'id':name,'bytes':len(raw),'label':'synthetic_diagnostic',
                            'path':'deterministic generator in predictive_lab/benchmark.py',
                            'sha256':hashlib.sha256(raw).hexdigest(),'raw':raw,'truncated':smoke})
    return sources


def entropy(data:bytes)->float:
    if not data:return 0.
    c=np.bincount(np.frombuffer(data,dtype=np.uint8),minlength=256)
    p=c[c>0]/len(data)
    return float(-(p*np.log2(p)).sum())


def measured(fn,repeat=2):
    out=fn()
    elapsed=[]
    for _ in range(repeat):
        t=time.perf_counter_ns()
        res=fn()
        assert res==out
        elapsed.append((time.perf_counter_ns()-t)/1e6)
    return out,round(statistics.median(elapsed),3)


def benchmark_one(s:dict,models, D:int, backend:str,repeats=2)->list[dict]:
    data=s['raw']
    zstd=Zstd()
    direct={
       'direct_zstd3':lambda:zstd.compress(data,3),
       'direct_zstd9':lambda:zstd.compress(data,9),
       'direct_brotli11':lambda:brotli.compress(data,quality=11),
       'direct_lzma6':lambda:lzma.compress(data,preset=6),
    }
    rows=[]
    for name,enc in direct.items():
        coded,e=measured(enc,repeats)
        decoder={
          'direct_zstd3':lambda:zstd.decompress(coded,len(data)),
          'direct_zstd9':lambda:zstd.decompress(coded,len(data)),
          'direct_brotli11':lambda:brotli.decompress(coded),
          'direct_lzma6':lambda:lzma.decompress(coded),
        }[name]
        out,d=measured(decoder,repeats)
        assert out==data
        rows.append({'id':s['id'],'label':s['label'],'sha256':s['sha256'],'input_bytes':len(data),
                     'D':D,'method':name,'encoded_bytes':len(coded),'encode_ms':e,'decode_ms':d,
                     'encode_MBps':round((len(data)/1e6)/(e/1e3),3) if e else None,
                     'decode_MBps':round((len(data)/1e6)/(d/1e3),3) if d else None,
                     'suffix_zstd9_bytes':None,'suffix_entropy':None,'residual_zstd9_bytes':None,
                     'residual_entropy':None,'model_coeff_bytes':0,'prefix_bytes':0,'model_selected':None,'exact':True})
    for model in models:
        if model=='auto':
            fun=lambda:pr.encode_adaptive(data,D,backend)
        else:fun=lambda model=model:pr.encode(data,model,D,backend)
        (coded,meta),e=measured(fun,repeats)
        out,d=measured(lambda:pr.decode(coded),repeats)
        assert out==data
        mode=meta['model'];prefix_frames=meta['prefix_frames']
        cf=pr.fit_coefficients(data,D,prefix_frames,mode)
        _,residual=pr.residual_stream(data,mode,D,prefix_frames,cf)
        suffix=data[prefix_frames*D:]
        row={'id':s['id'],'label':s['label'],'sha256':s['sha256'],'input_bytes':len(data),
             'D':D,'method':'predict_'+model+'_'+backend,'encoded_bytes':len(coded),
             'encode_ms':e,'decode_ms':d,
             'encode_MBps':round((len(data)/1e6)/(e/1e3),3) if e else None,
             'decode_MBps':round((len(data)/1e6)/(d/1e3),3) if d else None,
             'prefix_bytes':meta['prefix_bytes'],'suffix_entropy':round(entropy(suffix),5),
             'residual_entropy':round(entropy(residual),5),
             'suffix_zstd9_bytes':len(zstd.compress(suffix,9)),
             'residual_zstd9_bytes':len(zstd.compress(residual,9)),
             'model_coeff_bytes':meta['coefficient_bytes'],'model_selected':mode,
             'exact':True}
        rows.append(row)
    return rows


def summarize(rows:list[dict],metadata:dict)->str:
    lines=['# Causal hypercomplex multiplication prediction: benchmark results','',
           f'**Run**: {metadata["date"]} | Python {metadata["python"]} | NumPy {metadata["numpy"]}',
           f'**Protocol**: {metadata["protocol"]}',
           '', 'All codecs round-trip exactly. All trained predictors use ONLY a raw initial prefix ',
           'as training data (at most 25% of full blocks, capped at 256); the coefficients, compressed prefix, framing, and CRC are transmitted.',
           'Residual predictions at time t depend only on t-1 and t-2. No oracle per-file selection is counted as a deployed codec.',
           '', '| File | Direct Brotli11 (B) | Direct zstd9 (B) | Auto causal predictor (B) | CD nonlinear predictor (B) | CD vs Brotli (%) |',
           '|---|---:|---:|---:|---:|---:|']
    by={}
    for r in rows:by.setdefault((r['id'],r['D']),{})[r['method']]=r
    for (name,D),m in by.items():
        a=m['direct_brotli11']['encoded_bytes'];z=m['direct_zstd9']['encoded_bytes']
        b=next((v['encoded_bytes'] for k,v in m.items() if k.startswith('predict_auto_')),None)
        cd=next((v['encoded_bytes'] for k,v in m.items() if k.startswith('predict_ar2_cd_')),None)
        if cd is None:continue
        lines.append(f'| {name}, D={D} | {a:,} | {z:,} | {b:,} | {cd:,} | {(cd/a-1)*100:+.1f}% |')
    lines.extend(['','## Matched nonlinear ablations (D fixed, one coefficient per feature)',
                  '', '| Dataset | Linear AR2 | CD product | Unsigned XOR-convolution | Pointwise product | Circular convolution |',
                  '|---|---:|---:|---:|---:|---:|'])
    for (name,D),m in by.items():
        modes=['ar2','ar2_cd','ar2_xorconv','ar2_pointwise','ar2_circular']
        vals=[]
        for mod in modes:
            v=next((r['encoded_bytes'] for k,r in m.items() if k.startswith('predict_'+mod+'_')),None)
            vals.append(f'{v:,}' if v is not None else '—')
        lines.append('| '+name+', D='+str(D)+' | '+' | '.join(vals)+' |')
    lines.extend(['','## Important interpretation','','The synthetic CD-coupled control is generated using an actual Cayley-Dickson bilinear product between prior blocks. It is **not** held-out natural-world evidence. Its role is to test whether the model can recover a known nonlinear generative dependency. A model can show a lower held-out residual entropy but still lose on full frame size and throughput.','',
                  'The expected general-purpose baseline is full-file Brotli11/Zstd9/LZMA6; these are not matched to the predictive model\'s internal 25% prefix cost. This is intentionally conservative and tests deployable compression rather than only prediction error.','',
                  'The spectral construction in PR #5 was a signed Hadamard-type linear transformation. The CD bilinear feature here is different: it multiplies *two previous data vectors* using the actual Cayley-Dickson multiplication rules. Unsigned XOR-convolution, pointwise multiplication, and cyclic convolution use matched-capacity extra features.',''])
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=Path('/tmp/hc-predictive-result'))
    p.add_argument('--smoke',action='store_true')
    p.add_argument('--dim',type=int,default=16)
    p.add_argument('--backend',choices=pr.COMPRESSORS,default='zstd9')
    p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--natural-only',action='store_true')
    p.add_argument('--models',nargs='*',default=['auto',*pr.MODELS])
    args=p.parse_args()
    sources=corpus(args.smoke,args.natural_only)
    rows=[]
    # JIT compilation is excluded from timings, after the same initial warmup for all variants.
    test=np.arange(args.dim,dtype=np.int64)-128
    for x in (1,2,3,4,5):pr.pair_feature(test,test,x)
    for source in sources:
        print(f'[{source["label"]}] {source["id"]}: {source["bytes"]:,} bytes, {len(args.models)} models',flush=True)
        rows.extend(benchmark_one(source,args.models,args.dim,args.backend,args.repeats))
    metadata={'date':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
              'protocol':'causal-prefix-min-25pct-256frames-quantized-v1','python':platform.python_version(),
              'numpy':np.__version__,'platform':platform.platform(),
              'run':{'dim':args.dim,'backend':args.backend,'smoke':args.smoke,'repeats':args.repeats,'models':args.models},
              'source_manifest':[{k:v for k,v in a.items() if k!='raw'} for a in sources]}
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/'results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (args.out/'results.json').write_text(json.dumps({'meta':metadata,'measurements':rows},indent=2)+'\n')
    report=summarize(rows,metadata)
    (args.out/'SUMMARY.md').write_text(report)
    print(report,flush=True)

if __name__=='__main__':main()
