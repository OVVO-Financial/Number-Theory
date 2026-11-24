#!/usr/bin/env python3
"""Positive/negative controls; generated data, not real sensor measurements."""
import json,csv
from pathlib import Path
import numpy as np
from algebraic_codec import as_fixed
from benchmark_real import per_signal_scores

def run(out='results/local_controls'):
    output=Path(out);output.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261118)
    n=900
    # Meaningful smooth rotations with a time-varying axis; values are true
    # unit quaternions, but are generated numerically for a positive control.
    t=np.arange(n)
    theta=.018*t +.04*np.sin(.017*t)
    axis=np.stack([np.sin(.005*t),np.cos(.007*t),.15*np.ones(n)],axis=1)
    axis/=np.linalg.norm(axis,axis=1)[:,None]
    w=np.cos(theta/2)
    xyz=axis*np.sin(theta/2)[:,None]
    q=as_fixed(np.column_stack((w,xyz)))
    broken=q.copy();broken[:,0]=broken[rng.permutation(n),0]
    independent=as_fixed(np.c_[rng.uniform(-.5,.5,size=(n,4))])
    outputs=[]
    for name,data in [('constructed_unit_quaternion',q),('shuffled_w_control',broken),('independent_channels',independent)]:
        scores,_=per_signal_scores(data,name,None)
        outputs+=scores
        best={s['mode']:s['encoded_bytes'] for s in scores}
        print(name,'raw',best.get('raw'),'ar2',best.get('ar2'),'norm_ar2',best.get('norm_ar2'),'mul',best.get('mul'))
    with (output/'measurements.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(outputs[0]));w.writeheader();w.writerows(outputs)
    (output/'README.md').write_text('This data is constructed, not a real-world result. See --download for NanoBench.\n')
if __name__=='__main__':run()
