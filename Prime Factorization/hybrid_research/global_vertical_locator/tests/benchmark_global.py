#!/usr/bin/env python3
"""Matched balanced-domain factorization benchmark; no factors supplied to binaries.
Degree is chosen from the public N-only count, and is capped to bound memory.
A successful run must verify a nontrivial factor.
"""
from pathlib import Path
import json,subprocess,statistics,math
import sympy
ROOT=Path(__file__).resolve().parents[1]
BIN=ROOT/'global_vertical'
WHEEL=ROOT/'y_fermat_wheel'
CASES=[('manuscript-40',978508015703)]
for b in (40,44,48,50):
    p=int(sympy.nextprime(2**(b//2-1)))
    q=int(sympy.prevprime(1.9*p))
    assert 1<p<q<2*p and p*q%5 != 0
    CASES.append((f'balanced-{b}',p*q))

def run(cmd,timeout=45):
    proc=subprocess.run(list(map(str,cmd)),capture_output=True,text=True,check=True,timeout=timeout)
    return json.loads(proc.stdout.strip().splitlines()[0])
results=[]
for label,n in CASES:
    y=run([WHEEL,n,420,0,(3*n)//4])
    total=int(y['qunit10'])
    b=[]
    for method in ('wheel','wheel840','batch','batch840'):
        mtot=int(y['qunit2m']) if method.endswith('840') else total
        degree=min(2048,max(1,math.isqrt(mtot)))
        trials=[]
        for j in range(3):
            rec=run([BIN,method,n,2048,100000000,degree],timeout=90)
            assert rec['status']=='factor_found',(label,method,j,rec)
            f=int(rec['factor']);assert 1<f<n and n%f==0
            trials.append(rec['seconds'])
        b.append({'method':method,'degree':degree,'median_seconds':statistics.median(trials),'times_seconds':trials})
    results.append({'name':label,'n':str(n),'actual_bits':n.bit_length(),'balanced_assumption':True,'degrees':{x['method']:x['degree'] for x in b},
                    'Y_wheel':{'m':420,'Yclasses':y['Yresidue_count'],'paired_classes':y['paired_residue_count'],
                                'Y_density':y['Yresidue_count']/420,'total_integer_Y':y['Ycount'],
                                'Y_wheel_positions':y['Ywheel'],'q_all':y['qall'],'q_unit10':y['qunit10'],
                                'q_unit840':y['qunit2m']},'methods':b,'speedup_batch_over_wheel':b[0]['median_seconds']/b[2]['median_seconds'],'speedup_batch840_over_wheel840':b[1]['median_seconds']/b[3]['median_seconds']})
    print(label, 'bits',n.bit_length(),[(x['method'],round(1000*x['median_seconds'],3)) for x in b],'speedup batch10',round(results[-1]['speedup_batch_over_wheel'],2),'batch840',round(results[-1]['speedup_batch840_over_wheel840'],2),flush=True)
(ROOT/'results').mkdir(exist_ok=True)
(ROOT/'results'/'global_benchmark.json').write_text(json.dumps({'protocol':{'repeats':3,'same_balanced_q_interval':True,'same_factor_candidates':True,'baseline':'wheel-10 divisibility with GMP','batch':'shared normalized product tree GMP polynomial','comparison_scope':'not a GNFS comparison'},'cases':results},indent=2))
