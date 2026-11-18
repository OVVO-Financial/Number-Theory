#!/usr/bin/env python3
import json,subprocess,time,statistics,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXE=ROOT/'vertical_lattice'
def run(*args):
 t=time.perf_counter()
 d=json.loads(subprocess.check_output([str(EXE),*map(str,args)],text=True))
 return int(d['envelope_lattice_count']), (time.perf_counter()-t)*1000
cases={
 '40b':978508015703,
 '64b':2200000013*3400000033,
 '128b':14791069684487230471*14791069684487230513,
}
rows=[]
for label,N in cases.items():
 m=420; span=100000;total=3*N//4
 for step in (0,1,2):
  counts=[];wall=[]
  for i in range(80):
   A=i*(total-span)//79
   n,ms=run('count',N,m,A,A+span,step)
   counts.append(n);wall.append(ms)
  row={'case':label,'N_bits':N.bit_length(),'Y_windows':80,'Y_window_width':span,'wheel_modulus':m,'newton_steps':step,'empty_tiles':sum(c==0 for c in counts),'candidate_envelope_points':sum(counts),'median_process_ms_per_tile':round(statistics.median(wall),3),'max_count_per_tile':max(counts)}
  print(row,flush=True)
  rows.append(row)
(ROOT/'results'/'newton_window_benchmark.json').write_text(json.dumps({'experiment':'N-only deterministic 80 evenly spaced windows across Y in [0,3N/4], each width100000; no factor information used','note':'Process-startup included; counts are enclosure candidates, not genuine factors or full-domain complexity','results':rows},indent=2)+'\n')
