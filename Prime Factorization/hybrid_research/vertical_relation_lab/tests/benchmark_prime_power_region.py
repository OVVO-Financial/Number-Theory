#!/usr/bin/env python3
"""Finite performance/control test: exact mass regions versus brute smoothness scan.

Small-window microbenchmark, not a theorem about relation distribution or a GNFS comparison.
"""
import json,sys,time
from pathlib import Path
from math import isqrt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from prime_power_region import region_masses,no_smooth_certificate,primes_upto


def brute_empty(N,L,U,B,ps):
    for x in range(L,U+1):
        z=x*x-N
        for p in ps:
            while z%p==0:z//=p
        if z==1:return False
    return True


def main():
    N=978508015703
    R=isqrt(N)+1
    cases=[]
    for B in (7,101):
        if any(N%p==0 for p in primes_upto(B)):continue
        ps=primes_upto(B)
        for width in (1,4,16,64):
            windows=[(R+1000+83*j, R+1000+83*j+width-1) for j in range(120)]
            start=time.perf_counter()
            certificates=[no_smooth_certificate(N,L,U,region_masses(N,L,U,B)) for L,U in windows]
            elapsed_cert=time.perf_counter()-start
            start=time.perf_counter()
            actual=[brute_empty(N,L,U,B,ps) for L,U in windows]
            elapsed_brute=time.perf_counter()-start
            assert all(not c or e for c,e in zip(certificates,actual))
            item={'B':B,'window_width':width,'windows':len(windows),'certified_empty':sum(certificates),'actually_empty':sum(actual),'cert_time_ms':round(elapsed_cert*1000,3),'brute_time_ms':round(elapsed_brute*1000,3)}
            cases.append(item)
            print(item,flush=True)
    (ROOT/'results'/'prime_power_region_benchmark.json').write_text(json.dumps({'input':str(N),'cases':cases},indent=2)+'\n')

if __name__=='__main__':main()
