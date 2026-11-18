from lane_cover import families
from scipy.optimize import milp,LinearConstraint,Bounds
from scipy.sparse import lil_matrix
import numpy as np
for k,n in [(6,63),(7,127),(8,255),(9,511),(8,223),(8,7)]:
 U,fam=families(k,n);ind={m:i for i,m in enumerate(U)}
 a=lil_matrix((len(U),len(U)))
 for i,x in enumerate(U):
  for m in fam[x]:a[i,ind[m]]=1
 res=milp(np.ones(len(U)),integrality=np.ones(len(U)),bounds=Bounds(0,1),constraints=LinearConstraint(a.tocsr(),np.ones(len(U)),np.inf),options={'time_limit':15.0})
 sel=[m for i,m in enumerate(U) if res.x[i]>.5] if res.x is not None else []
 print('k',k,'n',n,'opt',res.fun,'selected',sel,'mod8-counts',{i:sum(m%8==i for m in sel) for i in [1,3,5,7]},flush=True)
