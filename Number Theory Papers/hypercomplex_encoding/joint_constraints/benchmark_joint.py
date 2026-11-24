#!/usr/bin/env python3
"""Held-out, exact joint compression ablations for public EuRoC state trajectories.

Calibrate ordinary temporal predictor choices and kinematic direction ONLY on first
chronological quarter of each dataset. Never tune on held-out suffix.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,platform,statistics,time,urllib.request
from dataclasses import dataclass
from pathlib import Path
import brotli,numpy as np
from joint_codec import (GROUPS,MODE_INDEX,SCALE,all_group_config,build_residuals,
                         decode,encode,parse_euroc_csv,serialize_residuals,synthetic)

REPO='https://raw.githubusercontent.com/rpng/open_vins/master/ov_data/euroc_mav/'
DATASET=['MH_01_easy','MH_03_medium','MH_05_difficult']


def ensure_source(root,name,allow_download):
    f=Path(root)/(name+'.csv')
    if f.exists():return f
    if not allow_download:raise FileNotFoundError(f)
    f.parent.mkdir(parents=True,exist_ok=True)
    with urllib.request.urlopen(urllib.request.Request(REPO+name+'.csv',headers={'User-Agent':'number-theory-joint-constraints/1.0'}),timeout=90) as r:
        f.write_bytes(r.read())
    return f


def choose_baseline(t,arr):
    # Four DISJOINT groups: position, orientation, velocity and biases.
    cfg={};detail={}
    for name,idx in GROUPS.items():
        vals=[]
        for mode in ('raw','delta','ar2'):
            tmp=build_residuals(t,arr,base=all_group_config(mode))[:,list(idx)]
            # Column-major varint, then common Brotli backend; model choice trained
            # solely on the calibration prefix, not the evaluation suffix.
            size=len(brotli.compress(serialize_residuals(tmp),quality=9))
            vals.append((size,mode))
        size,selected=min(vals,key=lambda x:(x[0],MODE_INDEX[x[1]]))
        cfg[name]=selected;detail[name]={k:s for s,k in vals}
    return cfg,detail


def choose_kinematic(t,arr,cfg):
    vals=[]
    for kin in ('p_from_v','v_from_p'):
        z=encode(t,arr,base=cfg,kinematic=kin)
        vals.append((len(z),kin))
    return min(vals)[1],{k:int(size) for size,k in vals}


def configurations(kin):
    return {
      'temporal':(False,False,'none',False),
      'N_norm':(True,False,'none',False),
      'K_position_from_velocity':(False,False,'p_from_v',False),
      'K_velocity_from_position':(False,False,'v_from_p',False),
      'K_selected':(False,False,kin,False),
      'R_quat_motion':(False,True,'none',False),
      'B_bias_static':(False,False,'none',True),
      'N+K':(True,False,kin,False),
      'N+R':(True,True,'none',False),
      'N+B':(True,False,'none',True),
      'K+R':(False,True,kin,False),
      'K+B':(False,False,kin,True),
      'N+K+R':(True,True,kin,False),
      'N+K+B':(True,False,kin,True),
      'N+R+B':(True,True,'none',True),
      'K+R+B':(False,True,kin,True),
      'N+K+R+B':(True,True,kin,True),
    }


def compare_one(name,t,arr,meta,outdir):
    # Fixed training-only tuning; test data compressed as independent stream.
    div=len(t)//4
    train_t,train_x=t[:div],arr[:div]
    test_t,test_x=t[div:],arr[div:]
    cfg,baselines=choose_baseline(train_t,train_x)
    direction,direction_sizes=choose_kinematic(train_t,train_x,cfg)
    modes=configurations(direction)
    stats=[];timings=[]
    for label,(n,r,k,b) in modes.items():
        start=time.perf_counter()
        z=encode(test_t,test_x,base=cfg,norm=n,rotation=r,kinematic=k,bias_delta=b)
        et=time.perf_counter()-start
        start=time.perf_counter()
        rt,rv=decode(z)
        dt=time.perf_counter()-start
        if not np.array_equal(rt,test_t) or not np.array_equal(rv,test_x):raise AssertionError((name,label))
        xres=build_residuals(test_t,test_x,cfg,n,r,k,b)
        # Cost for each group without cross-channel correlation: exploratory only.
        group_bytes={g:len(brotli.compress(serialize_residuals(xres[:,list(idx)]),quality=9)) for g,idx in GROUPS.items()}
        stats.append({'dataset':name,'mode':label,'test_rows':len(test_x),'bytes':len(z),
                      'encode_s':round(et,4),'decode_s':round(dt,4),
                      'position_bytes':group_bytes['position'],'quaternion_bytes':group_bytes['quaternion'],
                      'velocity_bytes':group_bytes['velocity'],'bias_bytes':group_bytes['bias'],
                      'exact':True,'codec_backend':'brotli9','quantized_sha256':hashlib.sha256(test_x.tobytes()).hexdigest(),
                      'source_sha256':meta['source_sha256']})
    # Brought-in timestamp bytes count unchanged across methods.
    settings={'dataset':name,'original_rows':len(arr),'train_rows':div,'test_rows':len(test_x),
              'quantized_scale':SCALE,'stride':meta.get('stride',1),
              'source_sha256':meta['source_sha256'],'trained_temporal_modes':cfg,
              'training_group_scores':baselines,'selected_kinematic':direction,
              'kinematic_training_bytes':direction_sizes,'reference_url':REPO+name+'.csv',
              'test_quantized_sha256':hashlib.sha256(test_x.tobytes()).hexdigest(),
              'median_timestamp_step_ns':int(np.median(np.diff(test_t))),
              'baseline_train_only':True}
    return stats,settings


def append_controls(t,x,cfg,kin,seed=1776):
    # Shuffle the norm-dependent q_w independently, and mismatch velocity and position.
    rng=np.random.default_rng(seed)
    c={}
    q=x.copy();q[:,3]=q[rng.permutation(len(q)),3]
    c['q_w_disconnected']=q
    v=x.copy();v[:,7:10]=v[rng.permutation(len(v)),7:10]
    c['velocity_disconnected']=v
    b=x.copy();b[:,10:16]=b[rng.permutation(len(b)),10:16]
    c['bias_disconnected']=b
    out=[]
    for label,s in c.items():
        for name,n,r,k,bd in [
           ('temporal',False,False,'none',False),
           ('N',True,False,'none',False),
           ('K',False,False,kin,False),
           ('B',False,False,'none',True),
           ('N+K+B',True,False,kin,True)]:
            try:
                z=encode(t,s,cfg,n,r,k,bd);rt,rv=decode(z)
                assert np.array_equal(rv,s) and np.array_equal(rt,t)
                out.append({'control':label,'mode':name,'bytes':len(z),'exact':True})
            except (OverflowError,ValueError) as ex:
                out.append({'control':label,'mode':name,'bytes':None,'exact':False,'error':str(ex)})
    return out


def write_report(scores,meta,controls,out):
    by={};
    for row in scores:by.setdefault(row['dataset'],{})[row['mode']]=row['bytes']
    lines=['# Joint algebraic-constraint discovery on EuRoC estimated motion states',
        '', 'All reported sizes are standalone **quantized-integer** lossless frames: 16 state channels at scale 1e6, timestamp varints, model flags, CRC32, and Brotli9 residual payload. Source CSV decimal text/float64 byte patterns are **not** preserved.',
        '', 'Each sequence is split chronologically: first 25% trains temporal modes and the preferred p/v kinematic direction; last 75% is evaluated without tuning.',
        '', '| Sequence | Samples | Temporal baseline | Quaternion norm N | Kinematic K | N+K | All N+K+R+B | Combined vs temporal |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for ds,d in by.items():
        s=next(r for r in scores if r['dataset']==ds)
        base=d['temporal'];joint=d['N+K+R+B']
        lines.append(f"| {ds} | {s['test_rows']:,} | {base:,} | {d['N_norm']:,} | {d['K_selected']:,} | {d['N+K']:,} | {joint:,} | {(100*(base-joint)/base):+.2f}% |")
    lines+=['','## Independent constraint family and incremental contribution','','* **N**: quaternion unit norm, reconstruct scalar coordinate from current vector part, then transmit exact correction.','* **K**: translation–velocity kinematics, using measured contemporaneous velocity to predict position, or measured position to predict velocity; choose direction exclusively on calibration prefix.','* **R**: temporal quaternion Hamilton-product recurrence, a separate motion assumption from the same-row norm.','* **B**: near-static gyro/accelerometer biases, predicting bias from its prior decoded value. This is a temporal prior, not an exact physical equality or a hypercomplex-specific result.','',
        'Rotation-matrix orthogonality and unit-quaternion conjugation are intentionally **not** counted as independent equations: they follow from the same orientation norm.',
        '', '## Per-sequence incremental net compressed bytes',
        '', '| Dataset | N gain | K gain after N | R gain after N+K | B gain after N+K+R |',
        '|---|---:|---:|---:|---:|']
    for ds,d in by.items():
        lines.append(f"| {ds} | {d['temporal']-d['N_norm']:+,} | {d['N_norm']-d['N+K']:+,} | {d['N+K']-d['N+K+R']:+,} | {d['N+K+R']-d['N+K+R+B']:+,} |")
    lines+=['', 'A positive gain means fewer bytes; a negative gain means the proposed constraint harms compression. **Independence of physical assumptions does not guarantee additive information savings**, especially when baseline temporal prediction already exploits the same regularity.', '', '## Negative controls','', '| Control | Temporal | N | K | B |', '|---|---:|---:|---:|---:|']
    controls_by={}
    for c in controls:controls_by.setdefault(c['control'],{})[c['mode']]=c['bytes']
    for n,c in controls_by.items():lines.append(f"| {n} | {c.get('temporal')} | {c.get('N')} | {c.get('K')} | {c.get('B')} |")
    lines+=['', '## Critical scope','', '- This is EuRoC *estimated ground-truth states*, not raw inertial observation channels. Estimates may be normalized and smoothed by the estimator.', '- Kinematics are approximate in sampled/estimated data; corrections are always transmitted, so no quantized coefficient is literally omitted unless its correction is zero.', '- The 16D state vector is ordinary position/quaternion/velocity/bias data; higher Cayley–Dickson dimension contributes no new saving.', '- The comparison is against a strong temporal predictor with the same bit framing/backend, not an exhaustive competition with PCL or specialized state-estimate codecs.', '- Dataset originals are downloaded by the reproducible script, with immutable SHA256s recorded in source_manifests.json.', '', '## Reproduction','', '```sh','python -m pip install -r requirements.txt','python -m unittest -v test_joint_codec','python benchmark_joint.py --download --stride 10 --datasets MH_01_easy MH_03_medium MH_05_difficult','```','']
    Path(out,'REPORT.md').write_text('\n'.join(lines))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data-dir',default='data')
    ap.add_argument('--out',default='results')
    ap.add_argument('--download',action='store_true')
    ap.add_argument('--stride',type=int,default=10)
    ap.add_argument('--max-rows',type=int,default=None)
    ap.add_argument('--datasets',nargs='+',default=DATASET)
    ap.add_argument('--synthetic',action='store_true')
    args=ap.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    rows=[];meta=[];controls=[]
    for name in (['synthetic'] if args.synthetic else args.datasets):
        if args.synthetic:
            ts,arr=synthetic(700)
            src={'source_sha256':hashlib.sha256(arr.tobytes()).hexdigest(),'rows':len(arr)}
        else:
            f=ensure_source(args.data_dir,name,args.download)
            ts,arr,src=parse_euroc_csv(f)
            ts,arr=ts[::args.stride],arr[::args.stride]
            if args.max_rows:ts,arr=ts[:args.max_rows],arr[:args.max_rows]
        if len(arr)<150:raise ValueError('Not enough samples for '+name)
        src['stride']=args.stride
        z,man=compare_one(name,ts,arr,src,out)
        rows+=z;meta.append(man)
        if name==(args.datasets[0] if not args.synthetic else 'synthetic'):
            a=len(ts)//4;tt=ts[a:];xx=arr[a:];
            controls=append_controls(tt,xx,man['trained_temporal_modes'],man['selected_kinematic'])
        d={r['mode']:r['bytes'] for r in z}
        print(name,'rows',len(arr),'base',d['temporal'],'N',d['N_norm'],'K',d['K_selected'],'joint',d['N+K+R+B'],flush=True)
    with (out/'results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    Path(out,'source_manifests.json').write_text(json.dumps(meta,indent=2)+'\n')
    Path(out,'negative_controls.json').write_text(json.dumps(controls,indent=2)+'\n')
    write_report(rows,meta,controls,out)
    print('saved',out.resolve(),flush=True)

if __name__=='__main__':main()
