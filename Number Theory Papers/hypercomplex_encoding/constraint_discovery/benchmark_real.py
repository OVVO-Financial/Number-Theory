#!/usr/bin/env python3
"""Real-world algebraic constraint discovery on recorded drone quaternion telemetry.

The benchmark is reproducible in a networked environment; optionally give local
CSV files with --csv instead of downloading third-party data.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,sys,time,urllib.request,statistics,random
from pathlib import Path
import numpy as np
from algebraic_codec import (
 SCALE,as_fixed,encode,decode,residuals,characteristics,MODES,mul,conjugate)

# Precommitted, deterministic trajectory families; select first filename
# alphabetically in each class without inspecting file values. Public source:
# https://github.com/syediu/nanobench (BSD-3-Clause)
SOURCE_REPO='syediu/nanobench'
SELECTED_CLASSES=['B2_circle_slow','B3_','B11_','C4_']
AGENT='HypercomplexConstraintResearch/1.0 (public research; reproducible test)'


def retrieve(url, dest):
    req=urllib.request.Request(url,headers={'User-Agent':AGENT,'Accept':'application/vnd.github+json'})
    with urllib.request.urlopen(req,timeout=60) as response:
        data=response.read()
    dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_bytes(data)
    return data


def find_public_flows(dest, limit=4):
    base=Path(dest);base.mkdir(parents=True,exist_ok=True)
    tree=json.loads(retrieve('https://api.github.com/repos/syediu/nanobench/git/trees/main?recursive=1',base/'github_tree.json'))
    allpaths=sorted(t['path'] for t in tree['tree'] if t.get('type')=='blob' and
                    t['path'].startswith('datasets/dataset/') and t['path'].endswith('.csv'))
    picks=[]
    for prefix in SELECTED_CLASSES:
        chosen=next((s for s in allpaths if Path(s).name.startswith(prefix)),None)
        if chosen is not None and chosen not in picks:picks.append(chosen)
    if len(picks)<min(limit,3):raise RuntimeError(f'Unexpected dataset layout; found {len(picks)} filenames')
    picks=picks[:limit]
    out=[]
    for p in picks:
        target=base/Path(p).name
        url='https://raw.githubusercontent.com/'+SOURCE_REPO+'/main/'+p
        if not target.exists():retrieve(url,target)
        out.append((target,url))
    return out


def load_series(path,max_rows=1600):
    path=Path(path)
    # The Vicon quaternions are scalar-last; internally use (w,x,y,z).
    # Preserve source CSV bytes and hash for provenance. Fixed-point stream is
    # a deliberately different object; encoding *that* stream is lossless.
    with path.open('r',newline='',encoding='utf-8-sig') as f:
        rows=csv.DictReader(f)
        fieldnames=rows.fieldnames or []
        qkeys=('qw','qx','qy','qz')
        controls=('px','py','pz','imu_gyro_x')
        if not all(k in fieldnames for k in qkeys):
            raise ValueError(f'Missing quaternion columns. Available: {fieldnames}')
        has_controls=all(k in fieldnames for k in controls)
        allq=[];allcontrol=[];row_total=0;nanrows=0
        for row in rows:
            row_total+=1
            try:
                q=np.array([float(row[k]) for k in qkeys],dtype=np.float64)
                if not np.all(np.isfinite(q)):
                    nanrows+=1;continue
                allq.append(q)
                if has_controls:
                    nonq=np.array([float(row[k]) for k in controls],dtype=np.float64)
                    allcontrol.append(nonq if np.all(np.isfinite(nonq)) else np.zeros(4))
            except (ValueError,TypeError,KeyError):
                nanrows+=1
        full=np.asarray(allq,dtype=np.float64)
        ctrl=np.asarray(allcontrol,dtype=np.float64) if has_controls else None
        if len(full)<120:raise ValueError('Need >=120 finite quaternion observations')
        # Equidistant sampling; no train/test switching after results; preserve
        # temporal order and include warmup in each split.
        if max_rows and len(full)>max_rows:
            start=min(max(0,len(full)//10),len(full)-max_rows)
            full=full[start:start+max_rows]
            if ctrl is not None:ctrl=ctrl[start:start+max_rows]
        fixed=as_fixed(full)
        control=None
        if ctrl is not None:
            try:control=as_fixed(ctrl)
            except ValueError:control=None
        return fixed,control,{'csv_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                              'source_rows':row_total,'skipped_rows':nanrows,'selected_rows':len(fixed),
                              'columns':fieldnames,'exact_target':'int32 round-to-nearest of float quaternions times 1e6'}


def empirical_entropy(v):
    _,cnt=np.unique(v,return_counts=True)
    p=cnt/cnt.sum()
    return float(-(p*np.log2(p)).sum())


def per_signal_scores(arr, names, out, warmup_fraction=.25):
    """One held-out suffix per signal; calibrate compressor/packing on prefix."""
    arr=np.asarray(arr,dtype='<i4')
    border=max(16,int(len(arr)*warmup_fraction))
    training=arr[:border]
    test=arr[border:]
    mode_list=('raw','delta','ar2','mul','norm_delta','norm_ar2','pointwise')
    encodings=('i32_row','i32_col','varint_col')
    compressors=('zlib','brotli')
    try:
        import zstandard  # noqa: F401
        compressors=('zlib','brotli','zstd')
    except ImportError:pass
    result=[]
    for mode in mode_list:
        if not names.endswith('_quaternion') and mode in ('mul','norm_delta','norm_ar2'):
            # Negative controls are measured under norm/multiplication too.
            pass
        calib=[]
        for packing in encodings:
            for compressor in compressors:
                try:
                    blob=encode(training,mode,packing,compressor)
                    calib.append((len(blob),packing,compressor))
                except (OverflowError,ValueError):pass
        if not calib:continue
        calib.sort(key=lambda t:t[0]);_,chosen_packing,chosen_compressor=calib[0]
        t0=time.perf_counter()
        try:blob=encode(test,mode,chosen_packing,chosen_compressor)
        except (OverflowError,ValueError):continue
        e_ms=(time.perf_counter()-t0)*1000
        t0=time.perf_counter()
        try:rec=decode(blob)
        except Exception as ex:raise AssertionError(f'Failed decoding {names}/{mode}: {ex}')
        d_ms=(time.perf_counter()-t0)*1000
        if not np.array_equal(rec,test):raise AssertionError(f'Lossy: {names}/{mode}')
        residual=residuals(test,mode)
        wres=residual[:,0]
        result.append({'dataset':names,'mode':mode,'rows_train':border,
                       'rows_heldout':len(test), 'original_bytes':len(test)*16,
                       'encoded_bytes':len(blob),'packing':chosen_packing,
                       'compressor':chosen_compressor,'calib_bytes':calib[0][0],
                       'encode_ms':round(e_ms,3),'decode_ms':round(d_ms,3),
                       'w_median_abs_error':float(np.median(np.abs(wres.astype(np.int64)))),
                       'w_entropy_bits_per_value':round(empirical_entropy(wres),4),
                       'all_median_abs_error':float(np.median(np.abs(residual.astype(np.int64)))),
                       'exact':True,
                       'sha256_heldout':hashlib.sha256(test.tobytes()).hexdigest()})
    return result,{'signal':names,'train_rows':border,'heldout_rows':len(test),
                    'quaternion_like_norm':characteristics(test) if names.endswith('_quaternion') else None}


def selected_signals(path,seed=20261118,max_rows=1600):
    q,ctrl,meta=load_series(path,max_rows=max_rows)
    rng=np.random.default_rng(seed)
    shuffled_w=q.copy()
    shuffled_w[:,0]=q[rng.permutation(len(q)),0]
    temporal_shuffle=q[rng.permutation(len(q))].copy()
    series=[('recorded_quaternion',q),('shuffled_w_control',shuffled_w),('time_shuffled_quaternion',temporal_shuffle)]
    if ctrl is not None:series.append(('nonquaternion_4d_control',ctrl))
    return series,meta


def compile_report(rows,details,manifest,output):
    lines=['# Real-data arithmetic-constraint experiment',
           '', 'Source: NanoBench Vicon ground-truth quaternion trajectories; https://github.com/syediu/nanobench',
           '', '**Exactness:** only fixed-point int32 quaternion samples are reconstructed bit-for-bit; the original CSV decimals are not preserved.',
           '**Holdout:** first 25% of each selected signal used only to choose packing/compression, last 75% independently measured; modes were specified before viewing holdout.',
           '**Controls:** time-shuffled quaternion rows, independently shuffled w, and non-quaternion position/gyro 4D stream.',
           '','## Held-out compressed bytes (includes 18-byte frame)', '']
    methods=('raw','delta','ar2','mul','norm_delta','norm_ar2','pointwise')
    lines.append('| Signal | '+' | '.join(methods)+' |')
    lines.append('|'+'---|'*(len(methods)+1))
    datasets=sorted({x['dataset'] for x in rows})
    for ds in datasets:
        lookup={x['mode']:x for x in rows if x['dataset']==ds}
        lines.append('| '+ds+' | '+' | '.join(str(lookup[m]['encoded_bytes']) if m in lookup else '—' for m in methods)+' |')
    lines+=['','## Source files','']
    for source in manifest:
        lines.append(f"- `{source['filename']}`: {source['rows']} selected rows; original CSV SHA256 `{source['sha256']}`; source URL {source['url']}")
    lines+=['','## Interpretation rules',
            '- A **norm residual** is not a zero-bit coordinate: its correction stream and sign are accounted for by the codec.',
            '- Even if a norm restriction helps, this does not show a novel advantage: unit-quaternion parameterizations (3-axis rotation vectors, axis-angle, etc.) are established alternatives.',
            '- A universal identity of a transformed coordinate system does not lower source entropy; only source-specific concentration can.',
            '- Relative multiplication is causal and uses only previously decoded samples; no original future data or free model tables are available to the decoder.',
            '- Results are from 1e-6 fixed-point samples, not verbatim CSV bytes or float64 bit patterns.',
            '- All competitor codecs use the same frame format, fixed scale, and calibrated packing/entropy-codec selection. Different models have different Python inference speed.',
            '', '## Experiment reproducibility',
            'Run `python benchmark_real.py --download --max-rows 1600 --out results/real`.',
            'Unit tests: `python -m unittest -v test_algebraic_codec.py`.',
            'The downloaded public CSV files are not redistributed by this PR.']
    Path(output).write_text('\n'.join(lines)+'\n')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--download',action='store_true',help='fetch public NanoBench files from GitHub')
    p.add_argument('--csv',nargs='*',default=[],help='local NanoBench-compatible CSVs')
    p.add_argument('--max-rows',type=int,default=1600)
    p.add_argument('--out',default='results/real')
    p.add_argument('--max-flights',type=int,default=4)
    args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    source=[]
    if args.download:source=find_public_flows(out/'source_data',limit=args.max_flights)
    for file in args.csv:source.append((Path(file), 'local-file:'+str(file)))
    if not source:raise SystemExit('Supply --download or --csv PATH')
    rows=[];details=[];manifest=[]
    for filepath,url in source:
        print('Processing:',filepath.name,flush=True)
        signals,meta=selected_signals(filepath,max_rows=args.max_rows)
        manifest.append({'filename':filepath.name,'url':url,'sha256':meta['csv_sha256'],'rows':meta['selected_rows'],
                         'total_rows':meta['source_rows'],'skipped_rows':meta['skipped_rows']})
        for short,mat in signals:
            res,stats=per_signal_scores(mat,filepath.stem+'_'+short,out)
            rows.extend(res);details.append(stats)
            score={r['mode']:r['encoded_bytes'] for r in res}
            print('   ',short,'rows',len(mat),'raw',score.get('raw'),'ar2',score.get('ar2'),'mul',score.get('mul'),
                  'norm-ar2',score.get('norm_ar2'),flush=True)
    with (out/'measurements.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/'manifest.json').write_text(json.dumps({'source_repo':SOURCE_REPO,'seed':20261118,'selected':manifest,
                   'details':details,'version':1},indent=2))
    compile_report(rows,details,manifest,out/'REPORT.md')
    print('Complete:',out/'REPORT.md',flush=True)

if __name__=='__main__':main()
