#!/usr/bin/env python3
"""Newton rational tangent/chord exact independent oracle."""
import json, math,random,subprocess
from fractions import Fraction
from pathlib import Path
EXE=Path(__file__).resolve().parents[1]/'vertical_lattice'
def run(*args):return json.loads(subprocess.check_output([str(EXE),*map(str,args)],text=True))
def upper(n,y,t):
 x=n*n+y*y;u=math.isqrt(x);u+=u*u<x;v=1
 for _ in range(t):u,v=u*u+x*v*v,2*u*v
 return Fraction(u,v)
def ceil(fr):return -((-fr.numerator)//fr.denominator)
def brute(n,m,A,B,steps,pairs):
 U0=upper(n,A,steps);U1=upper(n,B,steps);total=0
 for y in range(A,B+1):
  L=Fraction(n*n+A*y,1)/U0
  U=U0 if A==B else (U0*(B-y)+U1*(y-A))/(B-A)
  for h in range(ceil(L),int(U)+1):
   if (y%m,h%m) in pairs:total+=1
 return total
rng=random.Random(3021);checks=0
for n in [187,8051,978508015703,1009*10007,1000003*1000033]:
 for m in [12,20,60,140,420]:
  if math.gcd(n,2*m)!=1:continue
  pairs=set(map(tuple,run('sieve',n,m)['pairs']))
  for i in range(22):
   A=rng.randrange(0,100000)
   B=A+rng.randrange(0,50)
   for step in (0,1,2):
    expected=brute(n,m,A,B,step,pairs)
    result=int(run('count',n,m,A,B,step)['envelope_lattice_count'])
    assert result==expected,(n,m,A,B,step,result,expected)
    checks+=1
 for p,q in [(83,97),(752867,1299709)]:
  if p*q!=n:continue
  y=(q*q-p*p)//2
  for step in (0,1,2):
   assert int(run('count',n,m,y,y,step)['envelope_lattice_count'])>=1
   checks+=1
print('Newton enclosure independent tests passed:',checks)
