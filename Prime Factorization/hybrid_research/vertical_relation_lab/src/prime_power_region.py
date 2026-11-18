#!/usr/bin/env python3
"""Exact regionwise prime-power valuation mass for F(R)=R**2-N.

A proof-level, N-only arithmetic certificate. Root lifting and progression counts
avoid testing each R in the region. The *optional* exact smooth-count upper bound
materializes a potentially huge integer P; keep intervals moderate for that mode.
"""
import argparse
import json
from math import gcd, isqrt


def primes_upto(B):
    if B<2: return []
    mask=bytearray(b'\x01')*(B+1);mask[0:2]=b'\x00\x00'
    for p in range(2,isqrt(B)+1):
        if mask[p]: mask[p*p:B+1:p]=b'\x00'*(((B-p*p)//p)+1)
    return [p for p in range(2,B+1) if mask[p]]


def roots_mod_prime(N,p):
    if p==2: return [N&1] if N&1 else [0]
    return [r for r in range(p) if (r*r-N)%p==0]


def lift_roots(N,p,mod,roots):
    out=[]
    for r in roots:
        if p==2:
            for t in (0,1):
                cand=r+t*mod
                if (cand*cand-N)%(mod*2)==0:out.append(cand)
        else:
            # p does not divide N: derivative 2r is invertible mod p.
            t=(-((r*r-N)//mod)*pow(2*r,-1,p))%p
            out.append(r+t*mod)
    return sorted(set(out))


def ap_count(L,U,mod,r):
    return (U-r)//mod-(L-1-r)//mod


def region_masses(N,L,U,B):
    if N<3 or N%2==0 or L>U or L*L<=N:
        raise ValueError('N odd >1, L<=U, and L^2>N required')
    max_f=U*U-N
    masses={}
    for p in primes_upto(B):
        if gcd(N,p)>1:
            raise ValueError(f'factor precheck: {p} divides N')
        roots=roots_mod_prime(N,p)
        mod=p
        val=0
        while mod<=max_f and roots:
            val+=sum(ap_count(L,U,mod,r) for r in roots)
            if mod>max_f//p:break
            roots=lift_roots(N,p,mod,roots)
            mod*=p
        masses[p]=val
    return masses


def no_smooth_certificate(N,L,U,masses):
    """True implies zero B-smooth F(R) anywhere in [L,U], exactly."""
    min_f=L*L-N
    acc=1
    # Saturate at min_f; never construct giant products for a zero certificate.
    for p,v in masses.items():
        if not v:continue
        b=p
        while v:
            if v&1:acc=min(min_f,acc*b)
            if acc>=min_f:return False
            v//=2
            if v:b=min(min_f,b*b)
    return acc<min_f


def upper_bound_smooth_count(N,L,U,masses,max_span=20000):
    if U-L+1>max_span:
        raise ValueError('exact product bound disabled for large windows; use no_smooth_certificate')
    # P = product over R of (B-smooth part of F(R)), computed without visiting R.
    P=1
    for p,v in masses.items():
        if v:P*=pow(p,v)
    min_f=L*L-N
    if min_f<=1:return U-L+1
    bound=0
    while P>=min_f and bound <U-L+1:
        P//=min_f
        bound+=1
    return bound


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('N',type=int);p.add_argument('L',type=int);p.add_argument('U',type=int)
    p.add_argument('B',type=int);p.add_argument('--bound',action='store_true')
    a=p.parse_args()
    m=region_masses(a.N,a.L,a.U,a.B)
    res={'N':str(a.N),'L':str(a.L),'U':str(a.U),'B':a.B,
         'prime_power_valuation_mass':{str(k):str(v) for k,v in m.items()},
         'zero_smooth_certified':no_smooth_certificate(a.N,a.L,a.U,m)}
    if a.bound:res['max_smooth_relations']=upper_bound_smooth_count(a.N,a.L,a.U,m)
    print(json.dumps(res,sort_keys=True))

if __name__=='__main__':main()
