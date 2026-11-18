#!/usr/bin/env python3
"""Same K and B across conventional QS and three geometry-selected QS streams.

Times include relation factor-base trial-smoothness, matrix/gcd and candidate
mapping; excludes process startup. The methods do NOT have equal distinct-x
coverage. We explicitly report both raw candidates and distinct x.
"""
from pathlib import Path
import json
import subprocess
from statistics import median

ROOT=Path(__file__).resolve().parents[1]
EXE=ROOT/'relation_bench'
CASES=[187,8051,798607,978508015703,10967535067]
MODES=['QS','CTMWalkOneSide','GroupingBandOneSide','VerticalOneSide']


def run(n,k,b,mode):
    result=subprocess.run([str(EXE),str(n),str(k),str(b),mode,'--stop-on-factor'],capture_output=True,text=True,check=True)
    return json.loads(result.stdout)


def main():
    rows=[]
    for n in CASES:
        for mode in MODES:
            runs=[run(n,32768,1021,mode) for _ in range(3)]
            r=dict(runs[1]);r['time_us_median']=median(x['time_us'] for x in runs)
            r['repetitions']=len(runs);r['all_factor_results_agree']=len(set(x['factor_recovered'] for x in runs))==1
            f=int(r['factor_recovered']);r['factor_verified']=f!=0 and 1<f<n and n%f==0
            rows.append(r)
            print(n,mode,'distinct_x',r['unique_one_sided_x'],'rels',r['valid_relations'],'factor',r['factor_recovered'],'median_us',r['time_us_median'],flush=True)
    out=ROOT/'results'/'one_sided_benchmark.json'
    out.write_text(json.dumps({'protocol':{'K':32768,'B':1021,'repetitions':3,'same_smoothness_base':True,'methods_have_different_unique_x_coverage':True},'rows':rows},indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
