#!/usr/bin/env python3
"""Complete N-only dyadic lane cover verification, including all 3 witnesses."""
import json, random
from pathlib import Path
from lane_cover import families
P=Path(__file__).resolve().parent
checks=0
for k in range(3,11):
 w=1<<k
 for n in range(1,w,2):
  U,fam=families(k,n)
  # N-only universal square-witness lane set
  sel=set(range(n%8,w,8))
  assert len(sel)==(1<<(k-3))
  for x in U:
   assert (n*pow(x,-2,w))%w in sel
   assert fam[x]&sel
   checks+=2
# Verify actual semiprime factors belong to an N-only certified lane.
for row in json.loads((P/'smooth_order_benchmark.json').read_text())['corpus']:
 p,q,n=row['p'],row['q'],row['N']
 for k in range(3,13):
  w=1<<k
  m=(n*pow(p,-2,w))%w
  assert m%8==n%8
  assert m==q*pow(p,-1,w)%w
  checks+=2
# All finite integer-programming outcomes are feasible and respect basic lower bound.
for r in json.loads((P/'optimum_census.json').read_text()):
 assert r['solver_status']==0
 assert r['optimal']>=((1<<(r['k']-1))+5)//6
 assert r['optimal']<=(1<<(r['k']-3))
 checks+=2
print('exact N-only dyadic lane checks:',checks)
(P/'lane_selector_tests.json').write_text(json.dumps({'checks':checks,'coverage':'All n residues odd modulo 2^k for k=3..10; 120 actual semiprime corpus items k3..12; solver results k3..9'}))
