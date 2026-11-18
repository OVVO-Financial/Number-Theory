#!/usr/bin/env python3
"""Independent executable-level tests for batch localizers and inverse Y map."""
from pathlib import Path
import json,math,random,subprocess

BIN=Path(__file__).resolve().parents[1]/'global_vertical'
RNG=random.Random(613819)

def run(mode,n,*args):
    cp=subprocess.run([str(BIN),mode,str(n),*map(str,args)],capture_output=True,text=True,check=True,timeout=15)
    return json.loads(cp.stdout)

def prime(x):
    return x>=2 and all(x%d for d in range(2,math.isqrt(x)+1))

def check_factor(rec,n):
    factor=int(rec['factor'])
    if rec['status']=='factor_found':
        assert 1<factor<n and n%factor==0,rec
        p,q=sorted((factor,n//factor))
        assert int(rec['p'])==p and int(rec['q'])==q,rec
        assert int(rec['Y'])==(q*q-p*p)//2,rec
        assert int(rec['H'])==(q*q+p*p)//2,rec
        assert int(rec['Y'])%4==0
    else:
        assert factor==0,rec
        assert rec['status'] in ('budget_exhausted','no_factor_in_balanced_interval','small_prime'),rec
    if rec['status']=='no_factor_in_balanced_interval':
        assert int(rec['covered_qhi'])==int(rec['candidate_qhi']),rec
    elif rec['status']=='budget_exhausted':
        assert int(rec['covered_qhi']) < int(rec['candidate_qhi']),rec


def main():
    checks=0
    small_primes=[p for p in range(3,120) if prime(p)]
    semiprimes=[]
    for i,p in enumerate(small_primes):
        for q in small_primes[i:]:
            if q<2*p:semiprimes.append((p*q,p,q))
    RNG.shuffle(semiprimes)
    for n,p,q in semiprimes[:125]:
        for method in ('wheel','wheel840','batch','batch840'):
            r=run(method,n,128,10000,8)
            check_factor(r,n)
            assert r['status']=='factor_found',(n,p,q,method,r)
            checks+=5
    for n in (3,5,7,11,13,17,19,23,97,101,103,109):
        for method in ('wheel','wheel840','batch','batch840'):
            r=run(method,n,128,10000,8)
            check_factor(r,n)
            assert r['status']!='factor_found'
            checks+=2
    for n in (978508015703,522309809783,8356402666259,51505):
        for method in ('wheel','wheel840','batch','batch840'):
            r=run(method,n,2048,100000000,256)
            check_factor(r,n)
            assert r['status']=='factor_found'
            checks+=3
    for n in (978508015703,8356402666259):
        for method in ('wheel','wheel840','batch','batch840'):
            r=run(method,n,2048,1,8)
            check_factor(r,n)
            assert r['status']=='budget_exhausted'
            checks+=3
    for N in (143,187,8051,978508015703,522309809783):
        root=math.isqrt(N)
        for i in range(50):
            q=root+1+RNG.randrange(0,40000)
            record=run('map',N,q)
            yf=int(record['Yfloor']);yc=int(record['Yceil']);qb=int(record['qback_floor'])
            assert yf<=yc<=yf+1
            assert q-1<=qb<=q
            def cmp(q,y):return q**4-2*y*q*q-N*N
            assert cmp(q,yf)>=0 and (yc==yf or cmp(q,yc)<=0)
            checks+=4
    print(f'PASS {checks} full-coverage, factor recovery, budget, CRT, map and control checks')

if __name__=='__main__':main()
