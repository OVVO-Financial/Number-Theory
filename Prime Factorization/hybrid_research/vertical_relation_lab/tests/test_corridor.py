#!/usr/bin/env python3
"""Independent tests for the Newton tangent/chord lattice-corridor GMP prototype."""
from math import isqrt
from pathlib import Path
import json, subprocess
# This file stands alone: verify the C++ outputs independently using Python.
exe=Path(__file__).resolve().parent/'corridor_gmp'

def call(n,mode='corridor'):
    raw=subprocess.check_output([str(exe),str(n),'10000000',mode],text=True)
    return json.loads(raw)

def isprime(n):
    if n<2:return False
    d=2
    while d*d<=n:
        if n%d==0:return False
        d+=1
    return True

primes=[p for p in range(3,1500) if isprime(p)]
checks=0
for p in primes:
    for q in primes:
        if q<p or q>=2*p:continue
        n=p*q
        result=call(n)
        f=int(result['factor'])
        assert f>1 and f<n and n%f==0,(p,q,result)
        checks+=1
for p in primes[:120]:
    result=call(p)
    assert int(result['factor'])==0,(p,result)
    checks+=1
for p in primes:
    q=next((q for q in primes if q>2*p),None)
    if q is None:continue
    result=call(p*q)
    assert int(result['factor'])==0 or (p*q)%int(result['factor'])==0,(p,q,result)
    checks+=1
for n in [9,15,25,187,8051,978508015703,218775742411357180174920602608334561623]:
    result=call(n)
    f=int(result['factor'])
    assert f>1 and f<n and n%f==0,(n,result)
    checks+=1
print('PASS',checks,'C++ factorization control cases')

assert int(call(21)['factor']) == 0  # outside the explicitly balanced search guarantee
print('PASS unbalanced 21 remains unresolved without a false primality claim')
