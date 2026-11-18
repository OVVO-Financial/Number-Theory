from pathlib import Path
import json,math,random,collections
import sympy as s
P=Path(__file__).resolve().parent
rng=random.Random(220917)
checks=0
extra_hits=[]
controls=[]
for p in list(s.primerange(13,1500)):
    q=int(s.nextprime(p+7));n=p*q
    if p%4==1:
        witness=None
        for b in range(1,math.isqrt(p)+1):
            a=math.isqrt(p-b*b) if b*b<=p else 0
            if a*a+b*b==p:witness=(a,b);break
        assert witness is not None,p
        a,b=witness
        H=a*a+b*b
        E=a*a-b*b-n
        assert math.gcd(H,n)==p
        checks+=3
        if math.gcd(E,n)==1 and len(extra_hits)<20:extra_hits.append(dict(p=p,q=q,N=n,a=a,b=b,H=H,E=E))
    if p%4==3:
        for _ in range(50):
            a=rng.randint(1,p-1);b=rng.randint(1,p-1)
            assert (a*a+b*b)%p!=0
            checks+=1
    if p%4==1:
        roots=[x for x in range(1,p) if (x*x+1)%p==0]
        assert len(roots)==2
        checks+=1
for n in [8051,978508015703,309,798607]:
    p,q=[int(y) for y in s.factorint(n)];assert p*q==n
    r=math.isqrt(n);r+=(r*r<n)
    countE=0;countH=0;both=0
    # Sample a fixed-size, N-independent curve-near boundary band.
    for da in range(1,10001):
        a=r+da
        b=math.isqrt(a*a-n)
        E=a*a-b*b-n
        H=a*a+b*b
        e=math.gcd(E,n)
        h=math.gcd(H,n)
        okE=1<e<n;okH=1<h<n
        countE+=okE;countH+=okH;both+=okE and okH
        checks+=2
    controls.append(dict(N=n,countE=countE,countH=countH,overlap=both,window=10000))
(P/'norm_channel_checks.json').write_text(json.dumps(dict(exact_checks=checks,extra_factor_witnesses=extra_hits,near_curve_controls=controls),indent=2))
print('checks',checks,'extra witness',extra_hits[:2],'near curve',controls)
