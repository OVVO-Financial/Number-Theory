#!/usr/bin/env python3
import json, subprocess, math

def run(n,mode,*args):
    p=subprocess.run(['./ng_gmp',str(n),mode,*map(str,args)],text=True,capture_output=True,check=True)
    return json.loads(p.stdout)

for n in [15,21,25,187,798607,8051,978508015703,143,173*179,101*103,101*113,101*151]:
    s=run(n,'seeds',500000,2000000)
    assert s['status']=='factor_found',(n,s)
    d=int(s['factor']);assert 1<d<n and n%d==0
for n,seed in [(187,1),(798607,59),(798607,67),(978508015703,324927)]:
    s=run(n,'seed',seed,1000)
    assert s['status']=='factor_found',(n,seed,s)
for n in [187,8051,978508015703,101*103,101*113,101*151,1009*1013]:
    s=run(n,'groups',500000,500000)
    assert s['status']=='factor_found',(n,s)
for n in [187,8051,978508015703]:
    s=run(n,'seed',1,1)
    if s['status']=='factor_found':
        assert n % int(s['factor']) == 0
print('naive grouping GMP: all tests passed')
