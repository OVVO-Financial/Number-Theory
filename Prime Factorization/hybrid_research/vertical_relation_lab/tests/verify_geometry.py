from math import isqrt
# Independent Python reference for the exact rational corridor count.


def fs(n,m,a,b):
    # Sum_j floor((a*j+b)/m) for j=0,...,n-1, m>0.
    aa,a=divmod(a,m)
    bb,b=divmod(b,m)
    ans=aa*n*(n-1)//2 + bb*n
    while 1:
        if a >= m:
            ans += (n-1)*n*(a//m)//2
            a %= m
        if b >= m:
            ans += n*(b//m)
            b %= m
        y_max=a*n+b
        if y_max < m: break
        n=y_max//m
        b=y_max % m
        m,a=a,m
    return ans

def count_corridor(N, Q, H, p_mod=None, q_mod=None, mod=10):
    # q in [Q,Q+H] with q == q_mod mod mod; and p == p_mod.
    if p_mod is None:
        n=H+1
        return n + fs(n, Q*(Q+H), -N, N*(Q+H)) + fs(n,Q*Q,N,-N*Q)
    q0=Q+(q_mod-Q)%mod
    if q0>Q+H: return 0
    n=(Q+H-q0)//mod+1
    j0=q0-Q
    # upper: N*(Q+H-j)/(Q*(Q+H))
    # lower: N*(Q-j)/Q^2
    # count p == p_mod mod mod in [lower,upper]
    # floor((u-p_mod)/mod)+floor((p_mod-lower)/mod)+1
    den_up=mod*Q*(Q+H)
    den_low=mod*Q*Q
    A1=-N*mod
    B1=N*(Q+H-j0)-p_mod*Q*(Q+H)
    A2=N*mod
    B2=p_mod*Q*Q-N*(Q-j0)
    return n+fs(n,den_up,A1,B1)+fs(n,den_low,A2,B2)

def direct(N,Q,H,p_mod=None,q_mod=None,mod=10):
    c=0; good=[]
    for q in range(Q,Q+H+1):
        if q_mod is not None and q%mod != q_mod: continue
        lower_num=N*(Q-(q-Q)); lower_den=Q*Q
        upper_num=N*(Q+H-(q-Q)); upper_den=Q*(Q+H)
        lo=(lower_num+lower_den-1)//lower_den
        hi=upper_num//upper_den
        for p in range(max(1,lo),hi+1):
            if p_mod is not None and p%mod != p_mod: continue
            c+=1
            if p*q==N: good.append((p,q))
    return c,good


if __name__=='__main__':
    import json, random, subprocess
    from pathlib import Path
    rng=random.Random(20261117)
    exe=str(Path(__file__).resolve().parent/'corridor_gmp')
    count_direct=0
    for k in range(2500):
        N=rng.randrange(27,1000001,2)
        if N%5==0: continue
        Q=rng.randrange(isqrt(N)+1,10*isqrt(N)+1)
        H=rng.randrange(0,min(Q,90))
        C=count_corridor(N,Q,H)
        D,_=direct(N,Q,H)
        assert C==D,(N,Q,H,C,D)
        for rq in (1,3,7,9):
            rp=(N%10)*pow(rq,-1,10)%10
            lane=count_corridor(N,Q,H,rp,rq)
            d,_=direct(N,Q,H,rp,rq)
            assert lane==d,(N,Q,H,rp,rq,lane,d)
        if k<800:
            raw=subprocess.check_output([exe,str(N),'100000','count',str(Q),str(Q+H)],text=True)
            got=int(json.loads(raw)['count'])
            expected=sum(count_corridor(N,Q,H,(N%10)*pow(rq,-1,10)%10,rq) for rq in (1,3,7,9))
            assert got==expected,(N,Q,H,got,expected)
        count_direct+=1
    print('PASS:',count_direct,'independent unfiltered corridor checks +',4*count_direct,'paired-lane checks, plus C++ cross-checks')
