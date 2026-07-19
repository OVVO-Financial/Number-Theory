import sys,math,numpy as np
sys.path.append('/mnt/data')
from square_root_vaughan_attack import mobius_sieve,mangoldt_sieve,weights
M=120;U=2*M;sigma=.22;delta=.10;nv=161
mi,H,ell,v,g,W,Wv,err=weights(M,sigma,delta,nv)
# map m to row index
m0=M;m1=U
# weighted function features for each block m
dv=np.diff(v);tw=np.empty(nv);tw[0]=dv[0]/2;tw[-1]=dv[-1]/2;tw[1:-1]=(dv[:-1]+dv[1:])/2
wg=np.sqrt(tw*g)
Fblock=np.concatenate([W*wg,delta*Wv*wg],axis=1)
rs=np.arange(U+1,2*U+1)
nr=len(rs)
lam=mangoldt_sieve(U)
K=np.zeros((nr,nr))
for b in range(1,U+1):
 if lam[b]==0:continue
 idx=[];valid=[]
 for j,r in enumerate(rs):
  n=b*r
  m=int(math.isqrt(n-1)) # interval (m^2,(m+1)^2] for n; if n square, belongs previous m=n^.5-1
  if M<=m<=U:
   idx.append(m-M);valid.append(j)
 if not valid:continue
 Fb=Fblock[np.array(idx)]
 K[np.ix_(valid,valid)] += lam[b]*(Fb@Fb.T)
mu=mobius_sieve(2*U+1).astype(float)[U+1:2*U+1]
ev,E=np.linalg.eigh(K);o=np.argsort(ev)[::-1];ev=ev[o];E=E[:,o]
print('trace',np.trace(K),'muKmu',mu@K@mu,'ratio',mu@K@mu/np.trace(K),'munorm',mu@mu)
la=np.log(rs/U);ts=np.linspace(0,100,2001)
for j in range(12):
 u=E[:,j];best=(0,None,None)
 for t in ts:
  for typ,x in [('c',np.cos(t*la)),('s',np.sin(t*la))]:
   x=x-x.mean();n=np.linalg.norm(x)
   if n:
    c=abs(u@x/n)
    if c>best[0]:best=(c,t,typ)
 print(j,ev[j],best,'mucoeff2',(u@mu)**2)
