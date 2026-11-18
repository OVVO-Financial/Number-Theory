#!/usr/bin/env python3
"""N-only research dispatcher. Bounded, deterministic reproducible, NOT a primality certificate.
Stages: Newton(2-primary), Gaussian norm-one (p±1), standard Pollard p−1/Lucas p+1,
then bounded Brent-rho. Known polynomial product-tree fallback is available in v4 separately.
"""
import argparse,json,subprocess,time,pathlib
HERE=pathlib.Path(__file__).resolve().parent

def execute(n,B=256,newton_k=48,gauss_B=512,rho_budget=400000,quick_only=False):
    stages=[('newton_power2',['newton2_primary',str(n),str(newton_k),'4']),
            ('gaussian_norm_one',['gaussian_norm_one',str(n),str(gauss_B),'4','fixed']),
            ('pminus_stage1',['smooth_order',str(n),str(B),'pminus',str(rho_budget)]),
            ('lucas_stage1',['smooth_order',str(n),str(B),'lucas',str(rho_budget)])]
    if not quick_only: stages.append(('brent_rho',['smooth_order',str(n),str(B),'rho',str(rho_budget)]))
    start=time.perf_counter_ns()
    attempts=[]
    for stage,args in stages:
        t=time.perf_counter_ns()
        result=json.loads(subprocess.check_output([str(HERE/args[0])]+args[1:],text=True))
        ns=time.perf_counter_ns()-t
        factor=int(result['factor'])
        assert factor==0 or (1<factor<n and n%factor==0)
        attempts.append(dict(stage=stage,wall_us=round(ns/1000,1),kernel_us=result['us'],factor=str(factor)))
        if factor>1:break
    return {'N':str(n),'factor':attempts[-1]['factor'] if attempts else '0', 'elapsed_ms':round((time.perf_counter_ns()-start)/1e6,3), 'attempts':attempts, 'status':'factored' if attempts and int(attempts[-1]['factor'])>1 else 'unresolved'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('N',type=int)
    parser.add_argument('--B',type=int,default=256);parser.add_argument('--gauss-B',type=int,default=512)
    parser.add_argument('--rho-budget',type=int,default=400000);parser.add_argument('--newton-k',type=int,default=48)
    parser.add_argument('--quick-only',action='store_true')
    a=parser.parse_args();print(json.dumps(execute(a.N,B=a.B,newton_k=a.newton_k,gauss_B=a.gauss_B,rho_budget=a.rho_budget,quick_only=a.quick_only),indent=2))
