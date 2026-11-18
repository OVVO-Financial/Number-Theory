import math,random,json,subprocess,statistics
from pathlib import Path
P=Path(__file__).resolve().parent
rng=random.Random(7625)

def test():
    corpus=json.loads((P/'smooth_order_benchmark.json').read_text())['corpus']
    checks=0
    for row in corpus:
        p,q=row['p'],row['q'];n=p*q
        a=(p+q)//2;b=(q-p)//2
        h=a*a+b*b;y=2*a*b
        assert a*a-b*b==n and h*h==n*n+y*y
        assert math.gcd(h,n)==1
        inv=pow(h,-1,n)
        zr=n*inv%n;zi=y*inv%n
        assert zr==0 and zi*zi%n==1 and math.gcd(zi-1,n)==p and math.gcd(zi+1,n)==q
        # Projective norm-one coordinate from W=N+iY collapses to -1 mod N.
        norm=n*n+y*y
        invnorm=pow(norm,-1,n)
        Wr=(n*n-y*y)*invnorm%n
        Wi=2*n*y*invnorm%n
        assert Wr==n-1 and Wi==0
        checks+=7
        for _ in range(25):
            aa=a+rng.randint(-100,100)
            bb=b+rng.randint(-100,100)
            if aa<1 or bb<0:continue
            E=aa*aa-bb*bb-n;H=aa*aa+bb*bb;Y=2*aa*bb
            if math.gcd(H,n)!=1:continue
            invH=pow(H,-1,n)
            Ar=(n+E)*invH%n;Ai=Y*invH%n
            lhs_real=(Ar*Ar-Ai*Ai+1)%n
            lhs_im=(2*Ar*Ai)%n
            invH2=invH*invH%n
            rhs_real=(2*E*E*invH2)%n
            rhs_im=(2*E*Y*invH2)%n
            assert lhs_real==rhs_real and lhs_im==rhs_im
            # Since n is a distinct-prime semiprime, squaring does not alter
            # the prime support of the residual: the normalized 2ab test
            # contains precisely the same GCD information as the CTM error.
            Sp=Y*invH%n
            assert math.gcd(Sp*Sp-1,n)==math.gcd(E,n)
            checks+=3
    return checks

def run(n,upper,B=256,ns=4):
    o=subprocess.check_output([str(P/'gaussian_norm_one'),str(n),str(B),str(ns),upper],text=True)
    r=json.loads(o)
    factor=int(r['factor']);assert factor==0 or (factor>1 and factor<n and n%factor==0)
    return r

def bench():
    corpus=json.loads((P/'smooth_order_benchmark.json').read_text())['corpus']
    report=[]
    for row in corpus:
        for geometry in ['fixed','near']:
            for B in [256,1024]:
                row[f'gaussian_{geometry}_{B}']=run(row['N'],geometry,B,4)
    for kind in ['random','pminus_smooth','pplus_smooth']:
        for B in [256,1024]:
            for geometry in ['fixed','near']:
                entries=[row[f'gaussian_{geometry}_{B}'] for row in corpus if row['kind']==kind]
                hits=sum(int(r['factor'])>0 for r in entries)
                report.append(dict(kind=kind,B=B,geometry=geometry,hits=hits,total=len(entries),median_us=statistics.median(r['us'] for r in entries),median_hit_us=statistics.median(r['us'] for r in entries if int(r['factor'])>0) if hits else None))
    (P/'geometry_vs_fixed.json').write_text(json.dumps(dict(exact_identity_checks=test(),summary=report,corpus=corpus),indent=2))
    print(*report,sep='\n')

if __name__=='__main__':
 print('exact checks',test());bench()
