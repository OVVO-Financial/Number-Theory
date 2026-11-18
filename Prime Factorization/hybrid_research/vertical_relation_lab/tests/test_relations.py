#!/usr/bin/env python3
"""Independent Python oracle for smooth relation counts, modulo identities and GF(2) rank."""
import json,math,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXE=ROOT/'relation_bench'
def primes(B):
 return [p for p in range(2,B+1) if all(p%d for d in range(2,math.isqrt(p)+1))]
def exps(x,ps):
 if x==0:return None
 e=[];v=abs(x)
 for p in ps:
  c=0
  while v%p==0:v//=p;c+=1
  e.append(c)
 return e if v==1 else None
def oracle(N,K,B,mode):
 ps=primes(B);r0=math.isqrt(N);r0+=r0*r0<N
 nearplus=(mode=='NearPlus');nearminus=(mode=='NearMinus')
 low=r0-1;lo=low+(-low)%4;hi=math.isqrt(2*N);hi+=(hi*hi<2*N);hi+=1
 seen=set();raw=0;valid=0;direct=0;duplicates=0;evaluated=0;rows=[]
 lanes=[]
 if mode=='CTMWalk':
  for ar in range(10):
   for br in range(10):
    if (ar*ar-br*br-N)%20:continue
    a=r0+(ar-r0)%10;b=br
    if a<=b+1:continue
    lanes.append([a,b,a*a-b*b-N])
 for j in range(K):
  if mode in ('QS','CTM'):
   a=r0+j;f=a*a-N
   if mode=='QS':
    if math.isqrt(f)**2==f:
     bb=math.isqrt(f);d=a-bb
     if 1<d<N and N%d==0:direct+=1
     evaluated+=1;continue
    left=a*a;E=f
   else:
    b=math.isqrt(f);left=a*a-b*b;E=left-N
  elif mode=='CTMWalk':
   w=lanes[j%len(lanes)];a,b,E=w
   if a-b<=1:continue
   left=a*a-b*b
   if E<0:w[0]+=10;w[2]+=20*a+100
   else:w[1]+=10;w[2]-=20*b+100
  elif mode=='Grouping':
   g=2*j+3;left=g*(N//g);E=left-N
  elif mode=='GroupingBand':
   g=lo+4*j
   if g>hi:break
   left=g*(N//g);E=left-N
  elif mode=='SquareY':
   y=(3*N*(j+1)//(4*(K+1)));y-=y%4
   if y==0:continue
   h=math.isqrt(N*N+y*y);left=h*h-y*y;E=left-N*N
  else:
   E=j+1 if nearplus else -(j+1);left=N+E
  evaluated+=1
  if E==0:direct+=1;continue
  if mode!='QS':
   if E in seen:duplicates+=1;continue
   seen.add(E)
  ef=exps(E,ps)
  if ef is None:continue
  raw+=1
  lf=None if mode=='QS' else exps(left,ps)
  if mode!='QS' and (left<=0 or lf is None):continue
  assert (left-E)%N==0
  valid+=1
  if mode=='QS':mask=sum((e%2)<<i for i,e in enumerate(ef))
  else:
   mask=sum((e%2)<<i for i,e in enumerate(lf))
   mask+=sum((e%2)<<(len(ps)+i) for i,e in enumerate(ef))
  mask|=(E<0)<<(len(ps) if mode=='QS' else 2*len(ps))
  rows.append(mask)
 basis={}
 for x in rows:
  while x:
   j=x.bit_length()-1
   if j not in basis:basis[j]=x;break
   x^=basis[j]
 return dict(candidates_evaluated=evaluated,duplicate_residuals=duplicates,residual_smooth=raw,valid_relations=valid,rank=len(basis),dependencies=valid-len(basis),direct_exact_hits=direct)
checks=0
for N in (187,8051,798607,1009*10007,978508015703):
 for B in (23,53):
  for mode in ('QS','CTM','CTMWalk','Grouping','GroupingBand','SquareY','NearPlus','NearMinus'):
   K=200
   actual=json.loads(subprocess.check_output([str(EXE),str(N),str(K),str(B),mode],text=True))
   expected=oracle(N,K,B,mode)
   for key in expected:
    assert actual[key]==expected[key],(N,B,mode,key,expected[key],actual[key])
   recovered=int(actual['factor_recovered'])
   if recovered:assert 1<recovered<N and N%recovered==0
   checks+=1
print('Independent smooth relation + GF2 oracle cases passed:',checks)
