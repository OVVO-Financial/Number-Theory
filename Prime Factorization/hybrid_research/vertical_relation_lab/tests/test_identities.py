#!/usr/bin/env python3
import math,random
from sympy import divisor_count,divisors
rng=random.Random(8711)
for _ in range(25000):
 n=rng.randrange(3,10**20,2)
 a=rng.randrange(1,10**11)
 b=rng.randrange(0,a)
 e=a*a-b*b-n
 y=2*a*b
 h=a*a+b*b
 w=h*h-n*n-y*y
 assert w==e*(2*n+e)
 assert math.gcd(w,n)==math.gcd(e*e,n)
for n in [9,15,21,45,143,187,309,8051,978508015703]:
 points=[(u,n*n//u) for u in divisors(n*n) if u<=n]
 assert len(points)==(int(divisor_count(n*n))+1)//2
 true=sum(math.isqrt(u)**2==u and math.isqrt(v)**2==v for u,v in points)
 assert true==(int(divisor_count(n))+1)//2
print('25,000 exact residual/gcd identities + 9 divisor-lattice counts passed')
