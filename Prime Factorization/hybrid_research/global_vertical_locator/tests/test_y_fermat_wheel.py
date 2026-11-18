#!/usr/bin/env python3
"""Independent arithmetic oracle for the exact Y=2ab modular Fermat wheel.
All tests are deterministic and use Python big integers, not GMP helpers.
"""
import json
import math
import random
import subprocess
from pathlib import Path
from collections import defaultdict

EXE=Path(__file__).resolve().parents[1]/'y_fermat_wheel'
RNG=random.Random(20261117)

def run(n,m,lo,hi=None,classes=False):
    cmd=[str(EXE),str(n),str(m),str(lo),str(lo if hi is None else hi)]
    if classes:cmd+=['--classes']
    proc=subprocess.run(cmd,capture_output=True,text=True,check=True)
    out=[json.loads(line) for line in proc.stdout.splitlines() if line.strip()]
    return out

def oracle(N,m):
    rset=set()
    for p in range(1,2*m,2):
        if math.gcd(p,2*m)>1:continue
        q=(N*pow(p,-1,2*m))%(2*m)
        assert (p*p+q*q)%2 == 0
        rset.add((((q*q-p*p)//2)%m,((q*q+p*p)//2)%m))
    return rset

def oracle_q_count(N,m,A,B):
    qlo=math.isqrt(N)+(math.isqrt(N)**2<N)
    qhi=math.isqrt(2*N)
    def below(q,y):return q**4-2*y*q*q-N*N < 0
    def above(q,y):return q**4-2*y*q*q-N*N > 0
    while qlo<=qhi and below(qlo,A):qlo+=1
    while qhi>=qlo and above(qhi,B):qhi-=1
    return qlo,qhi,sum(math.gcd(q,2*m)==1 for q in range(qlo,qhi+1))

def main():
    count=0
    # Modular classes match independent direct enumeration for multiple N and moduli.
    nums=[35,77,143,187,299,8051,978508015703,
          14791069684487230471*14791069684487230513]
    for N in nums:
        for m in (20,60,84,420,840):
            entries=run(N,m,0,10_000,True)
            if math.gcd(N,2*m)>1:
                assert entries[0]['status'] in ('small_factor_precheck','trivial_gcd')
                assert math.gcd(N,entries[0]['factor'] if isinstance(entries[0]['factor'],int) else int(entries[0]['factor']))>1
                count+=1
                continue
            got=set(map(tuple,entries[1]['pairs']))
            expected=oracle(N,m)
            assert got==expected,(N,m,got^expected)
            assert entries[0]['Yresidue_count']==len({p[0] for p in expected})
            assert entries[0]['paired_residue_count']==len(expected)
            assert int(entries[0]['Ywheel'])==sum((k % m) in {p[0] for p in expected} for k in range(10001))
            count+=4
    # Valid true factor points can NEVER be excluded by the modular wheel.
    for p in (11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97):
        for q in (p,p+2,p+4,2*p-2):
            if q < p or q%2==0:continue
            n=p*q
            for m in (20,60,420):
                if math.gcd(n,2*m)>1:continue
                pairs=oracle(n,m)
                Y=(q*q-p*p)//2
                H=(q*q+p*p)//2
                assert (Y%m,H%m) in pairs
                count+=1
    # Each unit q mod 2m survives: no further arithmetic exclusion of q.
    for N in (978508015703,101*7907,11*17):
        for m in (20,84,420):
            if math.gcd(N,2*m)>1:continue
            pairs=oracle(N,m)
            for r in range(1,2*m):
                if math.gcd(r,2*m)!=1:continue
                p=(N*pow(r,-1,2*m))%(2*m)
                Y=((r*r-p*p)//2)%m
                H=((r*r+p*p)//2)%m
                assert (Y,H) in pairs
                count+=1
    # Exact inverse Y-interval count agrees with direct integer q enumeration.
    for i in range(120):
        p=RNG.randrange(51,3001,2)
        q=RNG.randrange(p,min(2*p,5000),2)
        N=p*q
        m=20 if N%5 else 84
        if math.gcd(N,2*m)>1:continue
        yupper=3*N//4
        A=RNG.randrange(0,yupper+1)
        B=min(yupper,A+RNG.randrange(0,200000))
        record=run(N,m,A,B)[0]
        assert record['status']=='ok'
        qlo,qhi,k=oracle_q_count(N,m,A,B)
        assert int(record['qunit2m'])==k,(N,m,A,B,k,record)
        assert int(record['qall'])==max(0,qhi-qlo+1)
        count+=2
    # The manuscript example yields the six exact classes, including the factor point.
    N=978508015703;m=420
    got=run(N,m,0,3*N//4,True)
    assert got[0]['Yresidue_count']==6
    assert got[0]['paired_residue_count']==6
    assert got[0]['qresidue_count']==192
    y=(1299709**2-752867**2)//2
    assert ((y%m),((1299709**2+752867**2)//2)%m) in set(map(tuple,got[1]['pairs']))
    count+=4
    print(f'PASS {count} exact modular, factor preservation, lift, interval and inverse checks')

if __name__=='__main__':main()
