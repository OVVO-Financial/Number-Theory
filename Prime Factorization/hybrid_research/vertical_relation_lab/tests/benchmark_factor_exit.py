#!/usr/bin/env python3
"""Stop at first independently verified nontrivial factor from a parity dependency."""
import json,subprocess,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'results'/'relation_benchmark.json').read_text())
EXE=ROOT/'relation_bench';K=32768;B=1021
modes=('QS','CTM','CTMWalk','Grouping','GroupingBand','NearPlus','NearMinus')
results=[]
for case in data['cases']:
 N=int(case['N']);name=case['case']
 for mode in modes:
  repetitions=[]
  for _ in range(3):
   d=json.loads(subprocess.check_output([str(EXE),str(N),str(K),str(B),mode,'--stop-on-factor'],text=True,timeout=120))
   if int(d['factor_recovered']):assert 1<int(d['factor_recovered'])<N and N%int(d['factor_recovered'])==0
   repetitions.append(d)
  d=repetitions[0]
  for dd in repetitions[1:]:
   assert dd['factor_recovered']==d['factor_recovered'] and dd['candidates_evaluated']==d['candidates_evaluated']
  r={'case':name,'bits':N.bit_length(),'method':mode,'factor_found':int(d['factor_recovered'])>0,'candidates_evaluated':d['candidates_evaluated'],'valid_relations_before_stop':d['valid_relations'],'dependencies_before_stop':d['dependencies'],'median_kernel_time_ms':round(statistics.median(dd['time_us'] for dd in repetitions)/1000,4)}
  print(r,flush=True);results.append(r)
(ROOT/'results'/'factor_exit_benchmark.json').write_text(json.dumps({'K_max':K,'B':B,'repeats':3,'criteria':'nontrivial, externally verified factor through direct hit or a valid double-square dependency','results':results},indent=2)+'\n')
