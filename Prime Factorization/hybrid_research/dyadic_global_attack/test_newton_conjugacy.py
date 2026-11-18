#!/usr/bin/env python3
import json,random,math,subprocess,statistics,collections
from pathlib import Path
P=Path(__file__).resolve().parent
rnd=random.Random(20261118)

def check_conjugacy():
    cases=0
    for row in json.loads((P/'smooth_order_benchmark.json').read_text())['corpus']:
        n=row['N']
        for _ in range(80):
            a=rnd.randrange(2,5000)
            b=rnd.randrange(1,a)
            H=a*a+b*b;Y=2*a*b
            if math.gcd(H,n)!=1 or math.gcd(H+Y,n)!=1 or math.gcd(Y,n)!=1:continue
            s=Y*pow(H,-1,n)%n
            assert math.gcd(s,n)==1 and math.gcd(s+1,n)==1
            t=(s-1)*pow(s+1,-1,n)%n
            d=-(a-b)**2*pow(a+b,-2,n)%n
            assert t==d
            s_next=(s+pow(s,-1,n))*pow(2,-1,n)%n
            assert (s_next-1)*pow(s_next+1,-1,n)%n==(t*t)%n
            # p- and q- Legendre structure of t is fixed by the Gaussian -1.
            for prime in [row['p'],row['q']]:
                if t%prime==0:continue
                symbol=pow(t,(prime-1)//2,prime)
                assert symbol==(-1 if prime%4==3 else 1)%prime
            cases+=5
    return cases

def run(n,k=48,ns=4):
    r=json.loads(subprocess.check_output([str(P/'newton2_primary'),str(n),str(k),str(ns)],text=True))
    d=int(r['factor']);assert d==0 or (1<d<n and n%d==0)
    return r

def bench():
    data=json.loads((P/'smooth_order_benchmark.json').read_text())['corpus']
    for r in data:r['newton2']=run(r['N'])
    report=[]
    for kind in ['random','pminus_smooth','pplus_smooth']:
        rows=[x['newton2'] for x in data if x['kind']==kind]
        succ=[x for x in rows if int(x['factor'])>1]
        report.append(dict(kind=kind,n=len(rows),hits=len(succ),median_us=statistics.median(x['us'] for x in rows),median_success_us=statistics.median(x['us'] for x in succ) if succ else None))
    out=dict(conjugacy_checks=check_conjugacy(),summary=report,corpus=data,special={})
    for p,q in [(257,263),(65537,65543),(65537,65599)]:
        n=p*q
        out['special'][str(n)]=run(n)
    (P/'newton2_benchmark.json').write_text(json.dumps(out,indent=2))
    print(out['conjugacy_checks'],report,out['special'])
if __name__=='__main__':bench()
