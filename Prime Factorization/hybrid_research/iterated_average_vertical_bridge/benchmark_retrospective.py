#!/usr/bin/env python3
"""Arithmetic operation-rank census, NOT a measured runtime benchmark."""
import bisect
import json
import math
import random
import statistics
from pathlib import Path
def primes_below(n):
    sieve=bytearray(b"\x01")*(n+1)
    sieve[0:2]=b"\x00\x00"
    for p in range(2,math.isqrt(n)+1):
        if sieve[p]: sieve[p*p:n+1:p]=b"\x00"*((n-p*p)//p+1)
    return [i for i in range(3,n+1,2) if sieve[i]]

def ceil_isqrt(n):
    r=math.isqrt(n)
    return r + (r*r<n)

def med(rows,key):return statistics.median(z[key] for z in rows)
primes=primes_below(1400000)
lo=bisect.bisect_left(primes,430000)
hi=bisect.bisect_left(primes,670000)
left=primes[lo:hi]
rng=random.Random(20261117)
rows=[]
for _ in range(600):
 p=left[rng.randrange(len(left))]
 low=bisect.bisect_right(primes,p)
 high=bisect.bisect_left(primes,2*p)
 q=primes[rng.randrange(low,high)]
 n=p*q;r=ceil_isqrt(n);a=(p+q)//2
 d=[];fractions=[]
 for k in range(1,15):
  w=1<<k;f,s=divmod((w-1)*n+r,w)
  d.append(min(f%p,p-f%p,f%q,q-f%q))
  fractions.append((s,w))
 chosen=max(2, max(range(1,7),key=lambda k: fractions[k-1][0]/fractions[k-1][1]))
 rows.append({'N':str(n),'p':p,'q':q,'best_k':d.index(min(d))+1,'j1':d[0], 'j6':d[5], 'j8':d[7], 'best6':min(d[:6]),'best14':min(d), 'chosen_k':chosen, 'chosen_j':d[chosen-1], 'fermat':a-r,'small_factor_scan_odd':(r-p)//2})
summary={'scope':'600 seeded distinct balanced odd semiprimes p ~430k..670k and p<q<2p',
 'method':'Exact retrospective nearest multiple of known p or q from each average; oracle measurement of ranks; not a factoring runtime benchmark; best_k cannot be selected without a locator',
 'medians':{k:med(rows,k) for k in ('j1','j6','j8','best6','best14','chosen_j','fermat','small_factor_scan_odd')},
 'proportions':{'source_fraction_rule_better_than_k1':sum(z['chosen_j']<z['j1'] for z in rows)/len(rows),
                'source_fraction_rule_better_than_fermat_positions':sum(z['chosen_j']<z['fermat'] for z in rows)/len(rows),
                'source_rule_chose_best_of_six':sum(z['chosen_j']==z['best6'] for z in rows)/len(rows)},
 'parallel_work':['For K levels tested to identical offset budget J, work is order 2*K*J candidate GCDs without batching','GCD arithmetic cost differs from Fermat square test cost'],
 'seed':20261117,'cases':len(rows)}
out=Path(__file__).with_name('results.json')
out.write_text(json.dumps({'summary':summary,'cases':rows},indent=2)+'\n')
print(json.dumps(summary,indent=2))
