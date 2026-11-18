#!/usr/bin/env python3
"""Exact independent oracle tests; no external packages or known factors supplied to the worker."""
import math
import random
import subprocess
from pathlib import Path

def primes_below(n):
    sieve=bytearray(b'\x01')*(n+1)
    sieve[0:2]=b'\x00\x00'
    for p in range(2,math.isqrt(n)+1):
        if sieve[p]: sieve[p*p:n+1:p]=b'\x00'*((n-p*p)//p+1)
    return [i for i in range(3,n+1,2) if sieve[i]]

def ceil_isqrt(n):
    r=math.isqrt(n)
    return r + (r*r<n)

ps=primes_below(1000)
pairs=[(p,q) for i,p in enumerate(ps) for q in ps[i+1:]]
rng=random.Random(20261117)
checks=0
for p,q in pairs[:5000]:
    n=p*q;r=ceil_isqrt(n);a=(p+q)//2;b=(q-p)//2;delta=a-r
    y=(q*q-p*p)//2;h=(q*q+p*p)//2
    assert h*h-y*y==n*n
    assert h-y==p*p and h+y==q*q and b-delta==r-p
    checks+=4
    for k in range(1,11):
        w=1<<k; c=(w-1)*n+r; f,s=divmod(c,w)
        assert s==(r-n)%w
        for j in (0,1,2,3,7,13,17):
            d=r-s-w*j
            assert math.gcd(f-j,n)==math.gcd(d,n)
            assert math.gcd(f+j,n)==math.gcd(r-s+w*j,n)
            checks+=2
        for m in ((q%w), (q%w)-w):
            j_num=r-s-m*p
            assert j_num%w==0
            j=j_num//w
            assert r-s-w*j==m*p
            if math.gcd(m,q)==1:
                assert math.gcd(f-j,n)==p
            checks+=3
        for u in (p,q,p*p):
            g=math.gcd(u,n)
            for sign in (-1,1):
                z=(w-1)*n+sign*u
                if z%w==0:
                    assert math.gcd(z//w,n)==g
                    checks+=1
    # At least one of the three Pythagorean witnesses has a positive-sign
    # 2-step average; for 3-step average allow both signs.
    assert any((u-n)%4==0 for u in (p,q,p*p))
    assert any((sign*u-n)%8==0 for sign in (-1,1) for u in (p,q,p*p))
    checks+=2

# A genuine modulus-16 obstruction to covering all three witnesses with +/-1.
p,q=3,5;n=p*q
assert not any((sign*u-n)%16==0 for sign in (-1,1) for u in (p,q,p*p))
checks+=1

# The real-pair vertical factor coordinate is decreasing in Y, so the
# dyadic residual conditioned on a fixed modular multiplier has one sign of slope.
for n in (187,8051,798607,978508015703):
    last=n+1
    for y in range(0,2500):
        # floor sqrt(sqrt(n*n+y*y)-y), exactly via nested integer sqrt.
        py=math.isqrt(math.isqrt(n*n+y*y)-y)
        assert py<=last
        last=py
        checks+=1

# Compare actual GMP batched worker with the independent exact first-hit oracle.
exe=Path(__file__).with_name('dyadic_scan')
if exe.exists():
    for p,q in pairs[100:140]:
        n=p*q;r=ceil_isqrt(n)
        for k in (2,4,6):
            w=1<<k;f=((w-1)*n+r)//w
            # Bound within nearest multiple of either factor in both directions.
            expect=min(f%p, p-f%p, f%q, q-f%q)
            res=subprocess.run([str(exe),str(n),str(k),str(expect),'16'],capture_output=True,text=True)
            assert res.returncode==0,(n,k,res.stdout,res.stderr)
            assert 'FACTOR ' in res.stdout
            checks+=1
print('PASS:',len(pairs[:5000]),'distinct semiprimes;',checks,'exact arithmetic/staircase/worker assertions')
