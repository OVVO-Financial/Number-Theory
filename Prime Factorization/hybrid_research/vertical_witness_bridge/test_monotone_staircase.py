#!/usr/bin/env python3
"""Exact no-floating-point verification of monotone u_Y and inverse Y thresholds."""
from math import isqrt,gcd
import json,random
rng=random.Random(14683)
vals=[3,15,21,187,309,8051,978508015703, 1009000007063]
checks=0

def u_level(N,Y): return isqrt(N*N+Y*Y)-Y

def t(N,k):return (N*N-k*k)//(2*k)

for N in vals:
 assert N%2==1
 for Y in range(0,500):
  a=u_level(N,Y);b=u_level(N,Y+1)
  assert b in (a,a-1)
  checks+=1
 for j in range(14000):
  k=rng.randint(1,N)
  threshold=t(N,k)
  for Y in [max(0,threshold-1),threshold,threshold+1, max(0,threshold+int(rng.randrange(-5,6)))]:
   assert (u_level(N,Y)>=k)==(Y<=threshold),(N,k,Y,threshold,u_level(N,Y))
   checks+=1
 for _ in range(800):
  k=rng.randint(1,N-1)
  low=t(N,k+1)+1;hi=t(N,k)
  if low<=hi:
   assert u_level(N,low)==k and u_level(N,hi)==k
   checks+=2
examples={}
for N,ks in [(187,[11,17,121]),(309,[3,9,103]),(978508015703,[1299709,566808719689])]:
 examples[str(N)]={}
 for k in ks:
  lo=t(N,k+1)+1;hi=t(N,k)
  examples[str(N)][str(k)]={'Y_first':lo,'Y_last':hi,'Y_multiplicity':hi-lo+1}
print(json.dumps({'checked_identities':checks,'examples':examples},indent=2))
