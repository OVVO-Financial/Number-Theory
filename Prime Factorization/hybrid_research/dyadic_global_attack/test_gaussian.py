from pathlib import Path
import subprocess, json,random,statistics,collections
from sympy import primerange,legendre_symbol
B=Path(__file__).resolve().parent
r=random.Random(1118)

def run(n,bound=256,seeds=4):
    v=subprocess.run([str(B/'gaussian_norm_one'),str(n),str(bound),str(seeds)],capture_output=True,text=True,check=True)
    o=json.loads(v.stdout)
    p=int(o['factor'])
    assert p==0 or (1<p<n and n%p==0)
    return o

def verify_gaussian():
    checks=0
    for p in list(primerange(13,2500)):
        for _ in range(4):
            a,b=r.randint(1,200),r.randint(1,200)
            if (a*a+b*b)%p==0: continue
            c=pow(a*a+b*b,-1,p)
            x=((a*a-b*b)*c)%p
            y=((2*a*b)*c)%p
            assert (x*x+y*y)%p==1
            order=p-(1 if p%4==1 else -1)
            A,Z=1,0
            U,V=x,y
            e=order
            while e:
                if e&1:A,Z=(A*U-Z*V)%p,(A*V+Z*U)%p
                U,V=(U*U-V*V)%p,(2*U*V)%p
                e//=2
            assert (A,Z)==(1,0), (p,a,b,order,A,Z)
            checks+=1
    return checks

def bench():
    corpus=json.loads((B/'smooth_order_benchmark.json').read_text())['corpus']
    for row in corpus:
        for bound in [256,1024]:
            row['gauss'+str(bound)]=run(row['N'],bound,4)
    counts=[]
    for kind in ['random','pminus_smooth','pplus_smooth']:
        for bound in [256,1024]:
            rows=[row['gauss'+str(bound)] for row in corpus if row['kind']==kind]
            hits=sum(int(x['factor'])>1 for x in rows)
            counts.append(dict(kind=kind,B=bound,n=len(rows),hits=hits,median_us=statistics.median(x['us'] for x in rows),median_success_us=statistics.median(x['us'] for x in rows if int(x['factor'])>1) if hits else None))
    (B/'gaussian_benchmark.json').write_text(json.dumps({'checks':verify_gaussian(),'summary':counts,'corpus':corpus,'warning':'These are standard Gaussian norm-one / p+1 smooth-order methods; not a factor-location breakthrough for generic balanced RSA.'},indent=2))
    print(*counts,sep='\n')

if __name__=='__main__':
    print('Gaussian prime identity checks:',verify_gaussian());bench()
