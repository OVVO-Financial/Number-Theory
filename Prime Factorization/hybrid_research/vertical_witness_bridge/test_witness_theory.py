#!/usr/bin/env python3
"""Independent exact oracle for N^2 Pythagorean witnesses on the squared complex line.

Factor pairs are supplied to TEST the classification, never to claim factor localization.
The two modular sieves deliberately differ: the original square-factor wheel
cannot be used when admitting extra (non-square) factor-revealing witnesses.
"""
from math import gcd, isqrt
import json

def primes_upto(M):
    out=[]
    for n in range(3,M+1,2):
        if all(n%p for p in out if p*p<=n):out.append(n)
    return out

def wheel(N,m,kind):
    z=2*m
    assert gcd(N,z)==1
    allowed=set()
    for r in range(1,z,2):
        if gcd(r,z)>1:continue
        if kind=='square':
            s=N*pow(r,-1,z)%z
            # r=p, s=q, so (H,Y)=((s^2+r^2)/2,(s^2-r^2)/2).
            h=(s*s+r*r)//2; y=(s*s-r*r)//2
        elif kind=='pythagorean':
            s=N*N*pow(r,-1,z)%z
            # r=H-Y, s=H+Y
            h=(s+r)//2; y=(s-r)//2
        else:raise ValueError(kind)
        allowed.add((y%m,h%m))
    return allowed

def witness_values(p,q):
    assert 3<=p<q and p%2 and q%2
    N=p*q; ulist=[1,p,q,p*p,N]
    assert len(set(ulist))==5
    rows=[]
    for u in sorted(ulist):
        assert N*N%u==0
        v=N*N//u
        H=(u+v)//2;Y=(v-u)//2
        assert H>=Y>=0 and H*H-Y*Y==N*N
        g=gcd(N,H-Y)
        assert g==gcd(N,u)
        assert (1<g<N)==(u in (p,q,p*p))
        rows.append({'u':u,'H':H,'Y':Y,'gcd_N_u':g,'gcd_N_Y':gcd(N,Y),
                     'square_ends':isqrt(u)**2==u and isqrt(v)**2==v})
    assert sum(1 for x in rows if 1<x['gcd_N_u']<N)==3
    return rows

def main():
    ps=primes_upto(151)
    cases=0; wheel_checks=0; extra_square_wheel_excluded=0; gcdYnontrivial=0
    for i,p in enumerate(ps):
        for q in ps[i+1:]:
            N=p*q;rows=witness_values(p,q)
            cases+=1
            for m in (20,60,420):
                if gcd(N,2*m)!=1:continue
                good=wheel(N,m,'pythagorean')
                strict=wheel(N,m,'square')
                for x in rows:
                    assert (x['Y']%m,x['H']%m) in good,(N,m,x)
                    wheel_checks+=1
                    if x['u'] in (p,q) and (x['Y']%m,x['H']%m) not in strict:
                        extra_square_wheel_excluded+=1
            gcdYnontrivial+=sum(1 for x in rows if x['u'] in (p,q) and 1<x['gcd_N_Y']<N)
    examples={}
    for p,q in [(3,103),(11,17),(83,97),(101,7907),(752867,1299709),(1009,1000000007)]:
        casestr=str(p*q)
        rows=witness_values(p,q)
        examples[casestr]={'p':p,'q':q,'original_Y':(q*q-p*p)//2,
                           'lowest_nontrivial_Y':min(x['Y'] for x in rows if 1<x['gcd_N_u']<p*q),
                           'three_nontrivial_witnesses':[x for x in rows if 1<x['gcd_N_u']<p*q]}
    res={'distinct_semiprime_cases':cases,'pythagorean_wheel_point_checks':wheel_checks,
         'extra_pythagorean_witnesses_rejected_by_original_square_wheel':extra_square_wheel_excluded,
         'nontrivial_gcd_Y_from_u_p_or_q':gcdYnontrivial,'examples':examples,
         'wheel_density_example_N_187':{str(m):{'square':len(wheel(187,m,'square')),
           'pythagorean':len(wheel(187,m,'pythagorean')),
           'pythagorean_Y_classes':len(set(y for y,h in wheel(187,m,'pythagorean')))} for m in (20,60,420)},
         'statement':'oracle verifies extra witnesses, not an N-only factor locator'}
    print(json.dumps(res,indent=2))
if __name__=='__main__':main()
