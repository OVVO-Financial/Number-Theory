#!/usr/bin/env python3
from __future__ import annotations
import math, csv, argparse
import numpy as np


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
    return mu.astype(float)


def weights(M:int,sigma:float=.22,delta:float=.10,nv:int=161):
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


def hnorm(blockvec,W,Wv,M,v,g,delta):
    s=(blockvec@W)/M; ds=(blockvec@Wv)/M
    return float(np.trapezoid((s*s+delta*delta*ds*ds)*g,v)),s,ds


def bernoulli_vals(x):
    t=x-np.floor(x)
    b1=t-0.5
    b2=t*t-t+1/6
    b3=t*t*t-1.5*t*t+0.5*t
    b4=t**4-2*t**3+t*t-1/30
    return b1,b2,b3,b4


def phi_exact(x):
    # sum_{1<=q<=x} log q = log Gamma(floor(x)+1)
    n=np.floor(x).astype(np.int64)
    # vector lgamma
    return np.array([math.lgamma(int(k)+1) for k in n],dtype=float)


def components(x,Hfourier=8):
    # Euler-Bernoulli components for Phi(x)=log Gamma(floor x+1), valid x>0.
    b1,b2,b3,b4=bernoulli_vals(x)
    smooth=x*np.log(x)-x+0.5*math.log(2*math.pi)
    c1=-b1*np.log(x)
    c2=b2/(2*x)
    c3=b3/(6*x*x)
    c4=b4/(12*x**3)
    exact=phi_exact(x)
    # Fourier split of right-continuous B1: B1^* + atom, with B1^*=0 at integers.
    low=np.zeros_like(x)
    for h in range(1,Hfourier+1):
        low -= np.sin(2*math.pi*h*x)/(math.pi*h)
    isint=np.isclose(x,np.rint(x),atol=1e-12,rtol=0)
    atom=-0.5*isint.astype(float)
    tail=b1-low-atom
    c1_low=-low*np.log(x)
    c1_atom=-atom*np.log(x)
    c1_tail=-tail*np.log(x)
    return {
        'exact':exact,'smooth':smooth,'B1':c1,'B2':c2,'B3':c3,'B4':c4,
        'B1_low':c1_low,'B1_atom':c1_atom,'B1_tail':c1_tail
    }


def run(M:int,Hfourier:int=8,sigma=.22,delta=.10,nv=161):
    U=2*M; X=(U+1)**2; amax=X//2
    mi,v,g,W,Wv=weights(M,sigma,delta,nv)
    mu=mobius_sieve(amax)
    names=['exact','smooth','B1','B2','B3','B4','B1_low','B1_atom','B1_tail']
    blocks={k:np.zeros(len(mi),float) for k in names}
    # Also residuals after successive Bernoulli orders.
    blocks.update({k:np.zeros(len(mi),float) for k in ['rem1','rem2','rem3','rem4']})
    # quotient-size buckets via Q=floor(X/a)
    bucket_edges=[(2,4),(5,8),(9,16),(17,32),(33,64),(65,U+1)]
    bucket_exact={b:np.zeros(len(mi),float) for b in bucket_edges}
    for a in range(U+1,amax+1):
        muv=mu[a]
        if muv==0: continue
        x=mi.astype(float)**2/a
        y=(mi.astype(float)+1)**2/a
        cx=components(x,Hfourier); cy=components(y,Hfourier)
        dif={k:cy[k]-cx[k] for k in names}
        for k in names: blocks[k]+=muv*dif[k]
        # exact residuals after smooth+B1+... orders
        blocks['rem1'] += muv*(dif['exact']-dif['smooth']-dif['B1'])
        blocks['rem2'] += muv*(dif['exact']-dif['smooth']-dif['B1']-dif['B2'])
        blocks['rem3'] += muv*(dif['exact']-dif['smooth']-dif['B1']-dif['B2']-dif['B3'])
        blocks['rem4'] += muv*(dif['exact']-dif['smooth']-dif['B1']-dif['B2']-dif['B3']-dif['B4'])
        Q=X//a
        for b in bucket_edges:
            if b[0]<=Q<=b[1]:
                bucket_exact[b]+=muv*dif['exact'];break
    energies={}
    signals={}
    for k,bv in blocks.items():
        energies[k],signals[k],_=hnorm(bv,W,Wv,M,v,g,delta)
    # Recombined approximations and their errors
    combos={
        'smooth+B1':blocks['smooth']+blocks['B1'],
        'smooth+B1+B2':blocks['smooth']+blocks['B1']+blocks['B2'],
        'smooth+B1+B2+B3':blocks['smooth']+blocks['B1']+blocks['B2']+blocks['B3'],
        'smooth+B1+B2+B3+B4':blocks['smooth']+blocks['B1']+blocks['B2']+blocks['B3']+blocks['B4'],
        'B1_low+atom':blocks['B1_low']+blocks['B1_atom'],
    }
    combo_rows=[]
    for k,bv in combos.items():
        e,_,_=hnorm(bv,W,Wv,M,v,g,delta)
        err,_,_=hnorm(blocks['exact']-bv,W,Wv,M,v,g,delta)
        combo_rows.append((k,e,err))
    bucket_rows=[]
    for b,bv in bucket_exact.items():
        e,_,_=hnorm(bv,W,Wv,M,v,g,delta)
        bucket_rows.append((b[0],b[1],e))
    # Gram cross table for main pieces: energy of sum vs sum energies
    main=['smooth','B1','B2','B3','rem3']
    diag=sum(energies[k] for k in main)
    return {'M':M,'U':U,'X':X,'amax':amax,'Hfourier':Hfourier,
            **{f'E_{k}':v for k,v in energies.items()},
            'main_diag':diag,'exact_over_main_diag':energies['exact']/diag if diag else np.nan,
            'combo_rows':combo_rows,'bucket_rows':bucket_rows}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--Ms',default='30,50,80,120,200');ap.add_argument('--H',type=int,default=8);ap.add_argument('--out',default='/mnt/data/bernoulli_extraction_results.csv');args=ap.parse_args()
    rows=[]; combos=[]; buckets=[]
    for M in map(int,args.Ms.split(',')):
        r=run(M,args.H); cr=r.pop('combo_rows');br=r.pop('bucket_rows');rows.append(r)
        for name,e,err in cr: combos.append({'M':M,'component':name,'energy':e,'error_energy':err})
        for lo,hi,e in br: buckets.append({'M':M,'Q_lo':lo,'Q_hi':hi,'energy':e})
        print('M',M,'exact',r['E_exact'],'smooth',r['E_smooth'],'B1',r['E_B1'],'B2',r['E_B2'],'rem2',r['E_rem2'])
    with open(args.out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    with open('/mnt/data/bernoulli_extraction_combos.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=combos[0].keys());w.writeheader();w.writerows(combos)
    with open('/mnt/data/bernoulli_extraction_buckets.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=buckets[0].keys());w.writeheader();w.writerows(buckets)

if __name__=='__main__': main()
