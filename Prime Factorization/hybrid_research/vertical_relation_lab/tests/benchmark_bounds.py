#!/usr/bin/env python3
import json,subprocess,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];EXE=ROOT/'relation_bench';K=32768
cases=json.loads((ROOT/'results'/'relation_benchmark.json').read_text())['cases']
chosen=[c for c in cases if c['case'] in ('40b_manuscript','48b','56b')]
records=[]
for c in chosen:
 for B in (257,1021,4093):
  for method in ('QS','CTM','CTMWalk','GroupingBand','SquareY','NearPlus'):
   samples=[]
   for _ in range(3):
    x=json.loads(subprocess.check_output([str(EXE),c['N'],str(K),str(B),method],text=True))
    if int(x['factor_recovered']):assert int(c['N'])%int(x['factor_recovered'])==0
    samples.append(x)
   d=samples[0]
   r={'case':c['case'],'bits':c['N_bits'],'B':B,'K':K,'method':method,'evaluated':d['candidates_evaluated'],'residual_smooth':d['residual_smooth'],'valid':d['valid_relations'],'rank':d['rank'],'dependencies':d['dependencies'],'factor':d['factor_recovered'],'median_ms':round(statistics.median(y['time_us'] for y in samples)/1000,3)}
   records.append(r)
   print(r,flush=True)
(ROOT/'results'/'smoothness_bound_sensitivity.json').write_text(json.dumps({'methods':'Same smoothness predicate and matrix for all, varying B on 3 balanced inputs','results':records},indent=2)+'\n')
