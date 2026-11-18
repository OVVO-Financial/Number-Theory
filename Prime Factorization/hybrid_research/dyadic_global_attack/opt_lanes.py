from lane_cover import families
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
import numpy as np, json

rows=[]
for k in range(3,10):
    mod=1<<k
    for n in [1,3,5,7]:
        U,fam=families(k,n); P={m:i for i,m in enumerate(U)}
        a=lil_matrix((len(U),len(U)))
        for i,x in enumerate(U):
            for m in fam[x]: a[i,P[m]]=1
        res=milp(np.ones(len(U)),integrality=np.ones(len(U)),bounds=Bounds(0,1),constraints=LinearConstraint(a.tocsr(),np.ones(len(U)),np.inf),options={'time_limit':15.0,'mip_rel_gap':0.0})
        obj=int(round(res.fun)) if res.fun is not None else None
        print(k,n,obj,res.message,flush=True)
        rows.append(dict(k=k,n_mod8=n,optimal=obj,lower_bound=res.mip_gap, solver_status=res.status))
with open('/mnt/data/dyadic_global_attack/optimum_census.json','w') as f:json.dump(rows,f,indent=2)
