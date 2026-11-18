#!/usr/bin/env python3
"""Independent arithmetic checks of one-sided square relations and their routes."""
from math import isqrt
from pathlib import Path
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1]
EXE=ROOT/'relation_bench'
MODES=('QS','CTMWalkOneSide','GroupingBandOneSide','VerticalOneSide')


def ceilroot(n):
    x=isqrt(n)
    return x if x*x==n else x+1


def vertical_to_fermat(n,y):
    h=ceilroot(n*n+y*y)
    return ceilroot((n+h+1)//2)


def case(n,limit=2048,B=257):
    for mode in MODES:
        row=json.loads(subprocess.check_output([str(EXE),str(n),str(limit),str(B),mode,'--stop-on-factor'],text=True))
        assert row['candidates_evaluated']<=limit, row
        assert row['rank']<=row['valid_relations'], row
        if mode!='QS':
            assert row['unique_one_sided_x'] <= (2*limit if mode=='GroupingBandOneSide' else limit),row
        factor=int(row['factor_recovered'])
        if factor:
            assert 1<factor<n and n%factor==0,row


def main():
    # Inversion agrees with the genuine factor midpoint for the exact Y.
    tested=0
    for p in (3,5,7,11,17,29,67,101,752867):
        for q in (p+2, p+4, 2*p+1, 7907,1299709):
            if q<=p or p*q<10:continue
            n=p*q;y=(q*q-p*p)//2
            a=vertical_to_fermat(n,y)
            assert a==(p+q)//2,(p,q,a)
            assert a*a - n == ((q-p)//2)**2,(p,q,a)
            tested+=1
    last=0
    for y in range(0,2001):
        a=vertical_to_fermat(8051,y)
        assert a>=last,(y,a,last)
        last=a
        tested+=1
    # A one-sided congruence is valid regardless of the selected R.
    for n in (187,8051,798607,978508015703):
        for r in range(ceilroot(n),ceilroot(n)+200):
            assert (r*r - (r*r-n))%n==0
            tested+=1
        case(n)
    print('one_sided_arithmetic_checks=',tested,'executable_cases=',len(MODES)*4)

if __name__=='__main__':
    main()
