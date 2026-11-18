import subprocess,json
from pathlib import Path
P=Path(__file__).resolve().parent
checks=0
for n,expected in [(8051,{83,97}),(978508015703,{752867,1299709}),(4295491591,{65537,65543})]:
    for mode in ['hybrid','rho']:
        o=json.loads(subprocess.check_output([str(P/'hybrid_native'),str(n),mode,'256','400000'],text=True))
        d=int(o['factor']);assert d in expected and n%d==0
        checks+=1
for n,k in [(67591,48),(4295491591,48),(4299161663,48)]:
    o=json.loads(subprocess.check_output([str(P/'newton2_primary'),str(n),str(k),'4'],text=True))
    assert 1<int(o['factor'])<n and n%int(o['factor'])==0
    checks+=1
for mode in ['E','H','both']:
    o=json.loads(subprocess.check_output([str(P/'dual_residual'),'243866858267','12000','256',mode],text=True))
    d = int(o['factor'])
    assert 1 < d < 243866858267 and 243866858267 % d == 0
    assert (mode=='E' and d==410323) or (mode!='E' and d==594329)
    checks+=1
for n in [101,1000003]:
    o=json.loads(subprocess.check_output([str(P/'hybrid_native'),str(n),'hybrid','32','2000'],text=True))
    assert int(o['factor'])==0 and o['method']=='unresolved'
    checks+=1
print('Native executable regression cases passed:',checks)
