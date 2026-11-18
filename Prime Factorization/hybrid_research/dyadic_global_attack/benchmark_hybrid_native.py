import json,statistics,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent
rows=json.loads((P/'smooth_order_benchmark.json').read_text())['corpus']
for row in rows:
    n=row['N']
    for mode in ['hybrid','rho','gauss','newton','pminus','lucas']:
        res=json.loads(subprocess.check_output([str(P/'hybrid_native'),str(n),mode,'256','400000'],text=True))
        p=int(res['factor']);assert p==0 or(1<p<n and n%p==0)
        row['native_'+mode]=res
results=[]
for kind in ['random','pminus_smooth','pplus_smooth']:
    group=[x for x in rows if x['kind']==kind]
    for mode in ['hybrid','rho','gauss','newton','pminus','lucas']:
        v=[x['native_'+mode] for x in group]
        results.append(dict(kind=kind,mode=mode,n=len(v),factored=sum(int(x['factor'])>0 for x in v),median_us=statistics.median(x['us'] for x in v),median_success_us=statistics.median(x['us'] for x in v if int(x['factor'])>0) if any(int(x['factor'])>0 for x in v) else None))
print(*results,sep='\n')
(P/'hybrid_native_benchmark.json').write_text(json.dumps({'source_corpus':'smooth_order_benchmark.json','B':256,'rho_budget':400000,'summary':results,'cases':rows,'caveats':'same in-process C++ executable per mode; compiled with GCC -O3 GMP, no GPU, no GNFS comparison. Corpus includes planted p±1 smoothness.'},indent=2))
