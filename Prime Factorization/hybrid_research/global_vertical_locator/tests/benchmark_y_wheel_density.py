#!/usr/bin/env python3
"""Exact full-Y-domain counts for direct Fermat sieve: one input, many moduli."""
import json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BIN=ROOT/'y_fermat_wheel'
N=978508015703
results=[]
for m in (20,60,420,840,2520,9240):
    out=subprocess.run([str(BIN),str(N),str(m),'0',str(3*N//4)],text=True,capture_output=True,check=True)
    x=json.loads(out.stdout)
    assert x['status']=='ok'
    yclasses=x['Yresidue_count']
    assert int(x['Ywheel'])<=int(x['Ycount'])
    results.append({'modulus':m,'classes_Y':yclasses,'density_Y':yclasses/m,'Y_wheel_count':x['Ywheel'],
                    'total_integer_Y':x['Ycount'],'paired_classes':x['paired_residue_count'],
                    'candidate_q_unit_2m':x['qunit2m'],'candidate_q_all':x['qall']})
    print('m',m,'Y classes',yclasses,'fraction',f'{yclasses/m:.6f}','q units',x['qunit2m'])
(ROOT/'results'/'y_wheel_density.json').write_text(json.dumps({'N':str(N),'Y_domain':['0',str(3*N//4)],
      'scope':'conditional balanced factor range only; modular exact factor-point sieve, not a residual sieve',
      'results':results},indent=2))
