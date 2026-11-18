import sympy as s, random, math, json, statistics
from pathlib import Path
P=Path(__file__).resolve().parent
rng=random.Random(74282)
rows=[]
for bits in [13,16,19]:
 for _ in range(100):
  p=int(s.nextprime(rng.randrange(1<<(bits-1),1<<bits)))
  q=int(s.nextprime(rng.randrange(p+10,min(2*p,p+p//2))))
  if q>=2*p:continue
  n=p*q;r=math.isqrt(n)+(math.isqrt(n)**2<n)
  limit=12000
  fE=fH=None;extra=0
  for j in range(limit):
   a=r+j;b=math.isqrt(a*a-n)
   E=a*a-b*b-n;H=a*a+b*b
   ge=math.gcd(E,n);gh=math.gcd(H,n)
   if fE is None and (1<ge<n or (E==0 and 1<a-b<n and n%(a-b)==0)):
    fE=(j,(a-b) if E==0 else ge)
   if fH is None and 1<gh<n:fH=(j,gh)
   if fE and fH:break
  dual=min([x[0] for x in (fE,fH) if x is not None],default=None)
  rows.append(dict(bits=bits,N=n,p=p,q=q,E=fE,H=fH,dual_idx=dual,H_wins=bool(fH and (not fE or fH[0]<fE[0]))))
summary=[]
for bits in [13,16,19]:
 cases=[r for r in rows if r['bits']==bits]
 E=[r for r in cases if r['E']]
 H=[r for r in cases if r['H']]
 wins=[r for r in cases if r['H_wins']]
 summary.append(dict(bits=bits,cases=len(cases),E_hits=len(E),H_hits=len(H),H_before_E=len(wins),median_E_index=statistics.median(r['E'][0] for r in E) if E else None,median_dual_index=statistics.median(r['dual_idx'] for r in cases if r['dual_idx'] is not None) if any(r['dual_idx'] is not None for r in cases) else None))
print(*summary,sep='\n')
(P/'dual_residual_census.json').write_text(json.dumps(dict(summary=summary,raw=rows,notes='N-only scan; factors supplied only to oracle to classify 1<gcd<N; 12,000 near-boundary states max per input; not optimized GMP benchmark.'),indent=2))
