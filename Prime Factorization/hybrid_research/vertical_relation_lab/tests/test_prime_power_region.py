#!/usr/bin/env python3
"""Independent brute oracle for exact region mass + zero/upper certificates."""
import sys
from pathlib import Path
from math import gcd, isqrt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from prime_power_region import primes_upto,region_masses,no_smooth_certificate,upper_bound_smooth_count


def brute_masses(N,L,U,B):
    out={p:0 for p in primes_upto(B)}
    for r in range(L,U+1):
        z=r*r-N
        for p in out:
            while z%p==0:
                out[p]+=1
                z//=p
    return out


def brute_smooth(N,L,U,B):
    primes=primes_upto(B)
    count=0
    for r in range(L,U+1):
        z=r*r-N
        for p in primes:
            while z%p==0:z//=p
        if z==1:count+=1
    return count


def main():
    checked=0;zeros=0
    for N in (187,8051,798607,978508015703):
        s=isqrt(N)+1
        for B in (3,7,19,41):
            if any(N%p==0 for p in primes_upto(B)):continue
            for offset in (0,3,9,20,31,100):
                for span in (0,1,2,5,16,47):
                    L=s+offset;U=L+span
                    m=region_masses(N,L,U,B)
                    assert m==brute_masses(N,L,U,B),(N,L,U,B,m,brute_masses(N,L,U,B))
                    actual=brute_smooth(N,L,U,B)
                    if no_smooth_certificate(N,L,U,m):
                        zeros+=1
                        assert actual==0,(N,L,U,B,actual)
                    ub=upper_bound_smooth_count(N,L,U,m)
                    assert actual<=ub<=span+1,(N,L,U,B,actual,ub)
                    checked+=1
    print('prime_power_region_cases=',checked,'certified_empty=',zeros)

if __name__=='__main__':main()
