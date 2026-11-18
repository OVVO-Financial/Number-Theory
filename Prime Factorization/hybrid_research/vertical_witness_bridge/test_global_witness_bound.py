#!/usr/bin/env python3
"""Exact integer-only check of N-only Y<=O(N^(4/3)) witness theorem."""
from math import gcd
import json

def icbrt(n):
    assert n>=0
    lo,hi=0,1 << ((n.bit_length()+2)//3)
    while lo+1<hi:
        m=(lo+hi)//2
        if m*m*m<=n:lo=m
        else:hi=m
    return lo

def ccbrt(n):
    k=icbrt(n)
    return k if k*k*k==n else k+1

def bound(N):
    c=ccbrt(N)
    return N*(c*c-1)//(2*c)

def primes(n):
    out=[]
    for k in range(3,n+1,2):
        if all(k%p for p in out if p*p<=k):out.append(k)
    return out

ps=primes(151)
checks=0; output=[]
for i,p in enumerate(ps):
    for q in ps[i+1:]:
        N=p*q
        square=(q*q-p*p)//2
        extra=q*(p*p-1)//2
        y=min(square,extra)
        u=max(q,p*p)
        assert y<=bound(N),(p,q,N,y,bound(N))
        assert u>=ccbrt(N*N),(p,q,u,ccbrt(N*N))
        assert u<N and gcd(u,N) in (p,q)
        assert ((extra<square)==(q>p*p))
        checks+=1
for p,q in [(3,103),(11,17),(83,97),(101,7907),(752867,1299709),(1009,1000000007),(11,127),(101,10211),(1009,1018091)]:
    N=p*q
    val=min((q*q-p*p)//2,q*(p*p-1)//2)
    b=bound(N)
    assert val<=b
    output.append({'p':p,'q':q,'N':N,'nearest_witness_Y':val,'N_only_integer_Y_bound':b,
                   'ratio_nearest_to_bound':round(val/b,9), 'threshold_branch':'nonsquare_q' if q>p*p else 'original_square'})
print(json.dumps({'checks':checks,'sample_cases':output},indent=2))
