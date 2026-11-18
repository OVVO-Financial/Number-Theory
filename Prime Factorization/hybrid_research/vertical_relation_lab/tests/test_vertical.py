#!/usr/bin/env python3
"""Independent integer-oracle tests for the vertical wheel and interval count."""
import json, math, random, subprocess
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXE=ROOT/'vertical_lattice'
def run(*args):
    return json.loads(subprocess.check_output([str(EXE),*map(str,args)],text=True))
def ceildiv(a,b):return -((-a)//b)
def cnt_brute(n,m,A,B,pairs):
    s0=math.isqrt(n*n+A*A);s0+=s0*s0<n*n+A*A
    s1=math.isqrt(n*n+B*B);s1+=s1*s1<n*n+B*B
    total=0
    for y in range(A,B+1):
        if A==B:
            lower=Fraction(n*n+A*y,s0);upper=Fraction(s0,1)
        else:
            lower=Fraction(n*n+A*y,s0)
            upper=Fraction(s0*(B-A)+(s1-s0)*(y-A),B-A)
        lh=ceildiv(lower.numerator,lower.denominator)
        uh=upper.numerator//upper.denominator
        for h in range(lh,uh+1):
            if (y%m,h%m) in pairs:total+=1
    return total
rng=random.Random(1424)
test_count=0
for n in [187,8051,978508015703,1009*10007,1000003*1000033]:
    for m in [4,12,20,60,140,420]:
        if math.gcd(n,2*m)!=1:continue
        pairs=set(map(tuple,run('sieve',n,m)['pairs']))
        for _ in range(17):
            A=rng.randrange(0,4000)
            B=A+rng.randrange(0,18)
            actual=cnt_brute(n,m,A,B,pairs)
            got=int(run('count',n,m,A,B)['envelope_lattice_count'])
            assert actual==got,(n,m,A,B,actual,got)
            test_count+=1
        for _ in range(20):
            p=rng.randrange(3,4000,2);q=rng.randrange(3,4000,2)
            if math.gcd(p*q,2*m)!=1 or (p*q)%5==0:continue
            y=(q*q-p*p)//2;h=(q*q+p*p)//2
            assert(y%m,h%m) in set(map(tuple,run('sieve',p*q,m)['pairs']))
            test_count+=1
for n,p,q in [(187,11,17),(8051,83,97),(309,3,103),(978508015703,752867,1299709)]:
    d=run('verify',n,p,q)
    assert d['valid'] and int(d['E'])==0 and int(d['squared_residual'])==0
    if n<10000 and math.gcd(n,840)==1:
        d=run('search',n,420,200000)
        assert d['status']=='factor_found' and int(d['p'])*int(d['q'])==n,d
    test_count+=1
for n in [15,25,35,65,105,5005]:
    for action, args in [('sieve', (420,)), ('search', (420, 1000))]:
        d=run(action,n,*args)
        assert d['status']=='precheck_factor' and int(d['factor'])==5 and int(d['cofactor'])*5==n, (n,d)
        test_count+=1
print('vertical lattice tests passed:',test_count)
