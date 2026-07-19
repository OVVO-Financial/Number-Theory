import sys, math, numpy as np
sys.path.append('/mnt/data')
from square_root_vaughan_attack import mobius_sieve, mangoldt_sieve, weights
M=120;sigma=.22;delta=.10;nv=161;U=2*M;Q=U*U
mi,H,ell,v,g,W,Wv,err=weights(M,sigma,delta,nv)
mm=mi*mi;mm1=(mi+1)*(mi+1)
dv=np.diff(v);tw=np.empty(nv);tw[0]=dv[0]/2;tw[-1]=dv[-1]/2;tw[1:-1]=(dv[:-1]+dv[1:U+1])/2
wg=np.sqrt(tw*g)
Aq=np.zeros((Q+1,nv));Aqv=np.zeros((Q+1,nv))
for q in range(1,Q+1):
 d=((mm%q)-(mm1%q))/q;Aq[q]=d@W;Aqv[q]=d@Wv
lam=mangoldt_sieve(U)
K=np.zeros((U,U))
for b in range(1,U+1):
 if lam[b]==0:continue
 idx=b*np.arange(1,U+1)
 F=np.concatenate([Aq[idx]*wg,delta*Aqv[idx]*wg],axis=1)
 K+=lam[b]*(F@F.T)
ev,E=np.linalg.eigh(K);o=np.argsort(ev)[::-1];ev=ev[o];E=E[:,o]
a=np.arange(1,U+1,dtype=float);la=np.log(a/U);ts=np.linspace(0,80,1601)
mu=mobius_sieve(U+1).astype(float)[1:U+1]
print('trace',np.trace(K),'muKmu',mu@K@mu)
for j in range(12):
 u=E[:,j];best=(0,None,None)
 for t in ts:
  for typ,x in [('c',np.cos(t*la)),('s',np.sin(t*la))]:
   x=x-x.mean();n=np.linalg.norm(x)
   if n:
    c=abs(u@x/n)
    if c>best[0]:best=(c,t,typ)
 print(j,ev[j],best,'mu coeff2',(u@mu)**2)
