#!/usr/bin/env python3
"""Repeated end-to-end timings for the standalone corridor and fair wheel-10 trial controls."""
import json, statistics, subprocess
from pathlib import Path
exe=Path(__file__).resolve().parent/'corridor_gmp'
inputs={
 'close_40b':1000036000099,
 'balanced_40b':978508015703,
 'balanced_46b':54000165000121,
 'balanced_54b':12000001540000049,
}
all_rows=[]
for label,N in inputs.items():
  pair={}
  for method in ['corridor','trial']:
    samples=[]
    for i in range(7):
      o=json.loads(subprocess.check_output([str(exe),str(N),'100000000',method],text=True))
      f=int(o['factor'])
      assert f>1 and f<N and N%f==0,(label,method,o)
      samples.append(o)
    median=statistics.median(s['seconds'] for s in samples)
    pair[method]=median
    all_rows.append({'case':label,'N':str(N),'method':method,'seconds_median':median,'samples':samples,'nodes_median':statistics.median(x['nodes'] for x in samples)})
  print(label,'corridor',round(pair['corridor']*1000,3),'ms','trial',round(pair['trial']*1000,3),'ms','trial/corridor',round(pair['trial']/pair['corridor'],3),'x',flush=True)
(Path(__file__).resolve().parent/'benchmarks.json').write_text(json.dumps(all_rows,indent=2))
