#!/usr/bin/env python3
from __future__ import annotations
import csv, math, argparse
import numpy as np
from pathlib import Path


def mobius_sieve(n:int)->np.ndarray:
    mu=np.zeros(n+1,dtype=np.int8); primes=[]; comp=np.zeros(n+1,dtype=bool);mu[1]=1
    for i in range(2,n+1):
        if not comp[i]: primes.append(i);mu[i]=-1
        for p in primes:
            v=i*p
            if v>n: break
            comp[v]=True
            if i%p==0: mu[v]=0;break
            mu[v]=-mu[i]
    return mu

def mangoldt_sieve(n:int)->np.ndarray:
    lam=np.zeros(n+1,float); prime=np.ones(n+1,bool);prime[:2]=False
    for p in range(2,n+1):
        if prime[p]:
            lp=math.log(p);q=p
            while q<=n:
                lam[q]=lp
                if q>n//p:break
                q*=p
            if p*p<=n: prime[p*p:n+1:p]=False
    return lam

def weights(M:int,sigma:float,delta:float,nv:int):
    m=np.arange(M,2*M+1,dtype=float);H=2*m+1;ell=np.log(m/M)
    v=np.linspace(-4*delta,4*delta,nv);g=np.exp(-.5*(v/delta)**2);g/=np.trapezoid(g,v)
    center=math.sqrt(M*(2*M+1));x=np.log(m/center)
    raw=np.exp(-.5*((x[:,None]-v[None,:])/sigma)**2)
    rawv=raw*((x[:,None]-v[None,:])/sigma**2)
    z=(m-1.5*M)/(.55*M);chi0=np.exp(-.5*z*z);chi1=chi0*z
    C=np.column_stack([chi0,chi1]);moment=np.vstack([H,H*ell]);A=moment@C
    coeff=np.linalg.solve(A,moment@raw);coeffv=np.linalg.solve(A,moment@rawv)
    W=raw-C@coeff;Wv=rawv-C@coeffv
    return m.astype(np.int64),v,g,W,Wv

def hnorm(signal,deriv,v,g,delta):
    return float(np.trapezoid((signal*signal+delta*delta*deriv*deriv)*g,v))

def run(M:int,sigma=.22,delta=.10,nv=161):
    U=2*M; X=(U+1)**2
    mi,v,g,W,Wv=weights(M,sigma,delta,nv)
    mu=mobius_sieve(X).astype(float); Mmu=np.cumsum(mu)
    lam=mangoldt_sieve(X); psi=np.cumsum(lam)
    lo=mi*mi; hi=(mi+1)*(mi+1)
    block=psi[hi]-psi[lo]
    # small-a Type I
    T1=np.zeros(len(mi))
    for j,(l,h) in enumerate(zip(lo,hi)):
        s=0.0
        for a in range(1,U+1):
            if mu[a]==0: continue
            q0=l//a+1; q1=h//a
            if q1>=q0:
                s += mu[a]*(math.lgamma(q1+1)-math.lgamma(q0))
        T1[j]=s
    # exact contracted large-a tail, resolved by complementary quotient q
    nq=U+1
    Z=np.zeros((nq+1,len(mi)))
    for q in range(2,nq+1):
        lq=lo//q
        hq=hi//q
        # a > U and l < aq <= h
        lower=np.maximum(lq,U)
        counts=np.where(hq>lower, Mmu[hq]-Mmu[lower], 0.0)
        Z[q]=math.log(q)*counts
    tail=Z.sum(axis=0)
    identity=float(np.max(np.abs(block-(T1+tail))))
    # wave energies
    sig=(tail@W)/M; dsig=(tail@Wv)/M
    totalE=hnorm(sig,dsig,v,g,delta)
    diagE=0.0
    qener=[]
    for q in range(2,nq+1):
        sq=(Z[q]@W)/M; dsq=(Z[q]@Wv)/M
        e=hnorm(sq,dsq,v,g,delta);diagE+=e;qener.append((q,e))
    cross=totalE-diagE
    # low/high q buckets
    buckets=[(2,max(2,U//8)),(max(3,U//8+1),U//4),(U//4+1,U//2),(U//2+1,U),(U+1,U+1)]
    bucket_rows=[]
    for a,b in buckets:
        if a>b: continue
        vec=Z[a:b+1].sum(axis=0)
        e=hnorm((vec@W)/M,(vec@Wv)/M,v,g,delta)
        d=sum(e0 for q,e0 in qener if a<=q<=b)
        bucket_rows.append((a,b,e,d,e-d))
    return {
        'M':M,'U':U,'X':X,'identity_error':identity,
        'tail_energy':totalE,'q_diagonal_energy':diagE,
        'tail_over_qdiag':totalE/diagE if diagE else float('nan'),
        'cross_q_energy':cross,
        'max_q_energy':max(e for _,e in qener),
        'max_q':max(qener,key=lambda z:z[1])[0],
        'bucket_rows':bucket_rows,
    }

def main():
    p=argparse.ArgumentParser();p.add_argument('--Ms',default='30,50,80,120,200');p.add_argument('--out',default='/mnt/data/lambda_mobius_contraction_results.csv');a=p.parse_args()
    rows=[]; buckets=[]
    for M in map(int,a.Ms.split(',')):
        r=run(M);b=r.pop('bucket_rows');rows.append(r)
        for x in b:buckets.append({'M':M,'q_lo':x[0],'q_hi':x[1],'bucket_energy':x[2],'bucket_diag':x[3],'bucket_cross':x[4]})
        print(M,r)
    with open(a.out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    with open('/mnt/data/lambda_mobius_contraction_buckets.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=buckets[0].keys());w.writeheader();w.writerows(buckets)
if __name__=='__main__':main()
