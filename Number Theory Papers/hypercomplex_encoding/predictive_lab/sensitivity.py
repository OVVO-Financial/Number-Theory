#!/usr/bin/env python3
"""Predeclared dimension sweep; exploratory, no per-file oracle configuration."""
from __future__ import annotations
import csv
from pathlib import Path
from predictive_lab import predictor as p
from predictive_lab.benchmark import corpus


def run(output:Path,smoke=False):
    natural=corpus(smoke=smoke,natural_only=True)
    dims=(4,8,16,32,64)
    models=('ar2','ar2_cd','ar2_xorconv','ar2_pointwise')
    rows=[]
    for x in natural:
        data=x['raw']
        for D in dims:
            for mode in models:
                frame,info=p.encode(data,mode,D,'zstd9')
                assert p.decode(frame)==data
                rows.append({'dataset':x['id'],'label':x['label'],'original_bytes':len(data),'D':D,
                            'method':mode,'output_bytes':len(frame),'lossless':True})
    output.mkdir(parents=True,exist_ok=True)
    with (output/'sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    lines=['# Dimension sweep (exploratory, not a selected deployable codec)','',
           'Fixed candidate set D=4,8,16,32,64. Every input is encoded with every setting, and every file round-tripped.','',
           '| Natural source | Smallest AR2 size (D) | Smallest CD size (D) | CD benefit at its best D |',
           '|---|---:|---:|---:|']
    for x in natural:
        a=[r for r in rows if r['dataset']==x['id'] and r['method']=='ar2']
        b=[r for r in rows if r['dataset']==x['id'] and r['method']=='ar2_cd']
        ra=min(a,key=lambda z:z['output_bytes']); rb=min(b,key=lambda z:z['output_bytes'])
        lines.append(f'| {x["id"]} | {ra["output_bytes"]:,} (D={ra["D"]}) | {rb["output_bytes"]:,} (D={rb["D"]}) | {rb["output_bytes"]-ra["output_bytes"]:+,} B |')
    lines+=['', 'These minima use the evaluation data to choose D and are **oracle diagnostics**, not valid deployable results. They are provided only to check whether larger hypercomplex dimensions plausibly change the conclusion. Use a single frozen D in the primary benchmark.','']
    (output/'SENSITIVITY.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
    return rows

if __name__=='__main__':
    import argparse
    q=argparse.ArgumentParser();q.add_argument('--out',type=Path,default=Path('/tmp/hc-predictive-sensitivity'));q.add_argument('--smoke',action='store_true')
    a=q.parse_args();run(a.out,a.smoke)
