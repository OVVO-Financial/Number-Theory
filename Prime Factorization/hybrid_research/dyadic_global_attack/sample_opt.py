from lane_cover import families
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
import numpy as np, json, random
rnd=random.Random(3042)
rows=[]
for k in range(5,10):
    mod=1<<k
    for residue in [1,3,5,7]:
        pool=list(range(residue,mod,8));rnd.shuffle(pool)
        for n in pool[:min(6,len(pool))]:
            U,fam=families(k,n);index={m:i for i,m in enumerate(U)}
            a=lil_matrix((len(U),len(U)))
            for i,x in enumerate(U):
                for m in fam[x]:a[i,index[m]]=1
            r=milp(np.ones(len(U)),integrality=np.ones(len(U)),bounds=Bounds(0,1),constraints=LinearConstraint(a.tocsr(),np.ones(len(U)),np.inf),options={'time_limit':10.0})
            rows.append((k,n,residue,int(round(r.fun)) if r.fun is not None else None,r.status))
            print(k,n,residue,rows[-1][-2],flush=True)
with open('/mnt/data/dyadic_global_attack/sample_optimum.json','w') as f:json.dump(rows,f,indent=2)
