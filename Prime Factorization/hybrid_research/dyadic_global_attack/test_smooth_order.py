#!/usr/bin/env python3
"""Independent exact verification and controlled experiment; factors used only in oracle."""
import sympy as s,subprocess,random,time,statistics,json,collections,math,os
from pathlib import Path
BASE=Path(__file__).resolve().parent
BIN=BASE/'smooth_order'
rnd=random.Random(20261118)

def run(n,B=256,mode='hybrid',budget=200000):
    x=subprocess.run([str(BIN),str(n),str(B),mode,str(budget)],text=True,capture_output=True,check=True)
    result=json.loads(x.stdout)
    factor=int(result['factor'])
    assert factor==0 or (1<factor<n and n%factor==0), (n,result)
    return result

def smooth_prime(bits,plus=False):
    allowed=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]
    for _ in range(15000):
        k=rnd.randint(5,min(bits-3,13))
        m=1<<k
        for j in range(20):
            m*=rnd.choice(allowed)
            if m.bit_length()>=bits:break
        if m.bit_length() not in (bits,bits+1):continue
        p=m-1 if plus else m+1
        if s.isprime(p):return int(p)
    raise ValueError('unable to construct smooth prime')

def corpus():
    rows=[]
    for bits in [18,22,26,30]:
        for kind in ['random','pminus_smooth','pplus_smooth']:
            for idx in range(10):
                if kind=='random':p=int(s.randprime(1<<(bits-1),1<<bits)) # randprime not reproducibly seeded
                else:p=smooth_prime(bits,plus=(kind=='pplus_smooth'))
                q=int(s.nextprime(p+rnd.randint(max(10,p//40),max(100,p//2))))
                if q>2*p:continue
                n=p*q
                rows.append(dict(kind=kind,bits=n.bit_length(),p=p,q=q,N=n, v2_qminus=int(s.factorint(q-1).get(2,0)), v2_pminus=int(s.factorint(p-1).get(2,0))))
    return rows

def bench():
    data=corpus()
    summary=collections.defaultdict(list)
    for r in data:
        n=r['N']
        for mode in ['hybrid','pminus','lucas','rho']:
            result=run(n,256,mode,400000)
            r[mode]=dict(factor=int(result['factor']),us=result['us'],gcd_calls=result['gcd_calls'],mod_steps=result['mod_steps'],method=result['method'])
            summary[(r['kind'],mode)].append(result)
    report=[]
    for (kind,mode),results in sorted(summary.items()):
        successes=[r for r in results if int(r['factor'])>1]
        report.append(dict(kind=kind,mode=mode,n=len(results),hits=len(successes),median_us=statistics.median(int(r['us']) for r in results),median_success_us=statistics.median(int(r['us']) for r in successes) if successes else None))
    result={'B':256,'rho_budget':400000,'corpus':data,'summary':report,'notes':'All measured times inside C++ process; corpus p,q supplied only to oracle, not worker; standard Pollard p-1/Lucas p+1/Brent rho methods, not a new algorithmic bound.'}
    (BASE/'smooth_order_benchmark.json').write_text(json.dumps(result,indent=2))
    print(*report,sep='\n')

if __name__=='__main__':bench()
