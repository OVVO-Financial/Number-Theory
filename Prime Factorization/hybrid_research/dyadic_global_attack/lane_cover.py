#!/usr/bin/env python3
"""N-only dyadic lane covering experiments. No factors provided to policies.
For units x=p mod 2^k and n=N mod 2^k, q=n/x, and the
three witness multipliers {p,q,q/p}={x,n/x,n/x^2}.
"""
from math import ceil
import json, argparse

def families(k, n):
    M=1<<k
    U=list(range(1,M,2))
    inv={x:pow(x,-1,M) for x in U}
    return U, {x:{x,n*inv[x]%M,n*inv[x]*inv[x]%M} for x in U}

def optimal_greedy(k,n):
    M=1<<k; U, fam=families(k,n)
    inv={x:set() for x in U}
    for x, opts in fam.items():
        for m in opts: inv[m].add(x)
    remaining=set(U)
    chosen=[]
    while remaining:
        m=max(U,key=lambda t:len(inv[t]&remaining))
        gain=inv[m]&remaining
        if not gain: raise ValueError('no cover')
        chosen.append(m)
        remaining.difference_update(gain)
    return len(chosen),chosen

def report():
    rows=[]
    for k in range(3,13):
        Ms=[1,3,5,7] if k<=10 else [1,3]
        for n in Ms:
            M=1<<k
            U,fam=families(k,n)
            s={x for x in U if x%8==n%8}
            assert all(fam[x]&s for x in U)
            count, sel=optimal_greedy(k,n)
            possible=set().union(*(fam[x] for x in U))
            coverage=[sum(m in fam[x] for x in U) for m in U]
            rootbound=6
            assert max(coverage)<=rootbound
            assert count>=ceil(len(U)/6)
            rows.append({'k':k,'N_mod_8':n%8,'unit_lanes':len(U),'square_lanes':len(s),'universal_lower_bound':ceil(len(U)/6),'greedy_cover_lanes':count, 'max_p_residues_per_lane':max(coverage),'sample_selected_lanes':sel[:24]})
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='lane_census.json');a=p.parse_args()
    rows=report()
    with open(a.output,'w') as f:json.dump({'interpretation':'greedy is not optimal; local unit residue coverage only, not factoring','rows':rows},f,indent=2)
    for r in rows:
        if r['N_mod_8']==1:print(r['k'],r['unit_lanes'],r['square_lanes'],r['universal_lower_bound'],r['greedy_cover_lanes'])
