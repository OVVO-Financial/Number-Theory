#!/usr/bin/env python3
"""Balanced semiprime, N-only relation-generator comparisons; no secret factors passed to kernels."""
import hashlib,json,os,platform,statistics,subprocess,time
from pathlib import Path
from sympy import nextprime,isprime
ROOT=Path(__file__).resolve().parents[1]
EXE=ROOT/'relation_bench'
K=32768
B=1021
methods=['QS','CTM','CTMWalk','Grouping','GroupingBand','SquareY','NearPlus','NearMinus']
def p(x):return int(nextprime(x))
cases={
 '28b':(p(11500),p(19000)),
 '32b':(p(42000),p(79000)),
 '36b':(p(190000),p(300000)),
 '40b_manuscript':(752867,1299709),
 '44b':(p(2400000),p(3900000)),
 '48b':(p(7900000),p(13200000)),
 '56b':(p(105000000),p(205000000)),
 '64b':(p(2200000000),p(3400000000)),
}
assert all(isprime(a) and isprime(b) and a<b<2*a for a,b in cases.values())
results=[]
for label,(pf,qf) in cases.items():
    N=pf*qf
    for mode in methods:
        outputs=[]
        wall=[]
        for rep in range(3):
            start=time.perf_counter()
            data=json.loads(subprocess.check_output([str(EXE),str(N),str(K),str(B),mode],text=True,timeout=120))
            wall.append((time.perf_counter()-start)*1e6)
            recovered=int(data['factor_recovered'])
            if recovered:assert 1<recovered<N and N%recovered==0
            outputs.append(data)
        nominal=outputs[0]
        for other in outputs[1:]:
            for key in ('residual_smooth','valid_relations','rank','dependencies','direct_exact_hits','factor_recovered','unique_residuals','duplicate_residuals'):
                assert other[key]==nominal[key],(label,mode,key,other[key],nominal[key])
        nominal['case']=label
        nominal['bits']=N.bit_length()
        nominal['kernel_time_ms_median']=round(statistics.median(d['time_us'] for d in outputs)/1000,4)
        nominal['process_wall_ms_median']=round(statistics.median(wall)/1000,4)
        nominal['independent_rank_per_kernel_second']=round(nominal['rank']*1e6/statistics.median(d['time_us'] for d in outputs),2)
        nominal['valid_relations_per_kernel_second']=round(nominal['valid_relations']*1e6/statistics.median(d['time_us'] for d in outputs),2)
        results.append(nominal)
        print(label,mode,'smooth',nominal['residual_smooth'],'valid',nominal['valid_relations'],'rank',nominal['rank'],'deps',nominal['dependencies'],'factor',nominal['factor_recovered'],'ms',nominal['kernel_time_ms_median'],flush=True)
obj={
 'methodology':{'input_generation':'sympy.nextprime deterministic seeds; factors withheld from program invocations','composites':len(cases),'K_per_method_case':K,'smoothness_bound':B,'timing_repetitions':3,'timing_scope':'C++ kernel: candidate construction, full smoothness trial division, rank/dependencies and up to 64 dependency extraction attempts; prime list and Engine allocation outside kernel timing','factorization_criteria':'only independently checked nontrivial gcd of valid square congruence is a factor; raw smooth residual is NOT a valid relation','comparator':'QS x^2-N smooth values; NearPlus/NearMinus compare direct residual enumeration; NOT a full CADO-NFS/GNFS comparison','geometry_assumption':'all selected semiprimes satisfy p<q<2p; only GroupingBand uses this to restrict group range','note':'SquareY uses uniformly spaced Y in [0,3N/4]; CTM and QS start at ceil(sqrt N); these policies are not matched factor-geography priors'},
 'environment':{'system':platform.platform(),'python':platform.python_version(),'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'executable_sha256':hashlib.sha256(EXE.read_bytes()).hexdigest(),'processor':platform.processor()},
 'cases':[{ 'case':label,'N':str(a*b),'N_bits':(a*b).bit_length(),'p_hidden_from_kernels':str(a),'q_hidden_from_kernels':str(b)} for label,(a,b) in cases.items()],
 'results':results
}
out=ROOT/'results'/'relation_benchmark.json';out.write_text(json.dumps(obj,indent=2)+'\n')
print('Saved',out)
