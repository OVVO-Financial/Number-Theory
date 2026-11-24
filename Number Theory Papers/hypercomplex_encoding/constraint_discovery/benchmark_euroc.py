#!/usr/bin/env python3
"""Reproduce held-out EuRoC estimated-trajectory tests and compressor comparisons.

Public source: https://github.com/uzh-rpg/rpg_trajectory_evaluation
Format: timestamp x y z qx qy qz qw. Estimates, NOT sensor ground truth.
Only 1e-6 quantized integer coordinates are preserved bit-for-bit.
"""
import argparse,csv,hashlib,json,math,struct,urllib.request,random
from pathlib import Path
import numpy as np
from algebraic_codec import SCALE,as_fixed,residuals,encode,decode

FILES={
 'MH_01_mono':'results/euroc_mono_stereo/laptop/vio_mono/laptop_vio_mono_MH_01/stamped_traj_estimate.txt',
 'MH_03_mono':'results/euroc_mono_stereo/laptop/vio_mono/laptop_vio_mono_MH_03/stamped_traj_estimate.txt',
 'MH_05_mono':'results/euroc_mono_stereo/laptop/vio_mono/laptop_vio_mono_MH_05/stamped_traj_estimate.txt',
 'MH_05_stereo':'results/euroc_mono_stereo/laptop/vio_stereo/laptop_vio_stereo_MH_05/stamped_traj_estimate.txt',
}
ROOT='https://raw.githubusercontent.com/uzh-rpg/rpg_trajectory_evaluation/master/'


def read_data(path):
    data=Path(path).read_bytes()
    rows=[]
    for line in data.decode('utf-8').splitlines():
        if not line.strip() or line.startswith('#'):continue
        parts=line.split()
        if len(parts)<8:continue
        try:items=list(map(float,parts[:8]))
        except ValueError:continue
        if all(np.isfinite(items)):
            rows.append(items)
    if len(rows)<100:raise ValueError('Insufficient trajectory rows: '+str(path))
    a=np.asarray(rows,dtype=np.float64)
    q=as_fixed(a[:,[7,4,5,6]])
    nonq=as_fixed(a[:,[1,2,3,7]])
    return q,nonq,hashlib.sha256(data).hexdigest()


def find_source(name,folder,download):
    file=Path(folder)/f'{name}.txt'
    if file.exists():return file
    if not download:raise FileNotFoundError(file)
    file.parent.mkdir(parents=True,exist_ok=True)
    url=ROOT+FILES[name]
    req=urllib.request.Request(url,headers={'User-Agent':'OVVO-reproducible-math-experiment/1.0'})
    with urllib.request.urlopen(req,timeout=60) as r:file.write_bytes(r.read())
    return file


def choose_codec(train,mode):
    candidates=[]
    for packing in ('i32_row','i32_col','varint_col'):
        for backend in ('zlib','brotli','zstd'):
            try:
                z=encode(train,mode,packing,backend)
                candidates.append((len(z),packing,backend))
            except ImportError:continue
            except OverflowError:continue
    if not candidates:raise RuntimeError(mode)
    return min(candidates,key=lambda x:x[0])


def empirical_h(x):
    _,ct=np.unique(x,return_counts=True)
    p=ct/ct.sum()
    return float(-(p*np.log2(p)).sum())


def score(q,name,source_sha):
    idx=len(q)//4
    train=q[:idx]
    test=q[idx:]
    out=[]
    for mode in ('raw','delta','ar2','mul','norm_delta','norm_ar2','pointwise'):
        _,pack,backend=choose_codec(train,mode)
        packet=encode(test,mode,pack,backend)
        if not np.array_equal(decode(packet),test):raise AssertionError(name+' '+mode)
        err=residuals(test,mode)
        out.append({'dataset':name,'rows':len(test),'mode':mode,'packet_bytes':len(packet),
         'raw_i32_bytes':len(test)*16,'packing':pack,'backend':backend,
         'w_entropy_bits':round(empirical_h(err[:,0]),3),
         'all_entropy_bits_per_quad':round(sum(empirical_h(err[:,i]) for i in range(4)),3),
         'quantized_sha256':hashlib.sha256(test.tobytes()).hexdigest(),
         'csv_sha256':source_sha,'roundtrip_exact':True})
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--download',action='store_true')
    p.add_argument('--folder',default='results/euroc_source')
    p.add_argument('--out',default='results/euroc')
    args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    scores=[];manifest=[]
    for name in FILES:
        path=find_source(name,args.folder,args.download)
        q,nonq,sha=read_data(path)
        rows={'name':name,'path':FILES[name],'csv_sha256':sha,'total_rows':len(q),'heldout_rows':len(q)-len(q)//4}
        manifest.append(rows)
        n=len(q);split=n//4
        tmp=q.copy()
        shuffled=q[split:].copy();size=len(shuffled)
        for i in range(size):shuffled[i,0]=q[split+(i*643+113)%size,0]
        scoremap=[('recorded_quaternion',q),
         ('nonquaternion_position_w_control',nonq),
         # place shuffled 75% after copy of 25% training prefix
         ('independently_shuffled_w_control',np.vstack([q[:split],shuffled])),
         ('time_shuffled_quaternion',np.vstack([q[:split],q[split:][np.random.default_rng(20261118).permutation(size)]]))]
        for label,mat in scoremap:
            res=score(mat,name+'_'+label,sha)
            scores.extend(res)
        dd={r['mode']:r['packet_bytes'] for r in scores if r['dataset']==name+'_recorded_quaternion'}
        print(name, 'heldout',len(q)-split,'raw',dd.get('raw'),'AR2',dd.get('ar2'),
              'mul',dd.get('mul'),'norm_AR2',dd.get('norm_ar2'),flush=True)
    with (out/'results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(scores[0]));writer.writeheader();writer.writerows(scores)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(out/'results.csv')

if __name__=='__main__':main()
