"""Reproducible SOTA compression evaluation: train on dev files, freeze, test heldout.

Usage: python -m sota_benchmark.benchmark calibrate --manifest ... --config ...
       python -m sota_benchmark.benchmark evaluate --manifest ... --config ... --output ...
Do not tune configurations on the test split or select the per-file oracle as a deployable codec.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import os
import platform
import statistics
import threading
import time
import tracemalloc
from dataclasses import asdict
from pathlib import Path
from typing import Callable

import numpy as np
import psutil

from hypercomplex_cascade_test import cd_cascade_inverse_checked
from .codec import Zstd, HybridConfig, VARIANT, compress_named, decompress_named, hybrid_encode, hybrid_decode, basis_codebook_encode, basis_codebook_decode

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST=Path(__file__).with_name('corpus_manifest.json')
STANDARD=['zstd1','zstd3','zstd9','zstd19','brotli5','brotli9','brotli11',
          'lzma6','lzma9','zlib6','zlib9','zstd_dict3','zstd_dict9']
VARIANTS=list(VARIANT)
BLOCK_SIZES=[32,128,512]  # in physical bytes; D=8*block_bytes, D-1 data bits/block


def digest(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def get_manifest(path:Path)->dict:
    m=json.loads(path.read_text())
    ids=[];files=[]
    for split in ('calibration','heldout'):
        for row in m[split]:
            ids.append(row['id'])
            files.append(str((PROJECT/row['path']).resolve()))
    if len(set(ids))!=len(ids) or len(set(files))!=len(files):
        raise ValueError('manifest leaks files across calibration and heldout')
    return m


def read_split(manifest:dict,split:str,smoke:bool=False,max_bytes:int=0,
               permit_missing:bool=False)->list[dict]:
    output=[]
    for entry in manifest[split]:
        path=(PROJECT/entry['path']).resolve()
        if not path.is_file():
            if permit_missing:
                continue
            raise FileNotFoundError(f'missing {split} file: {path}')
        full=path.read_bytes()
        data=full[:max_bytes] if max_bytes else full
        if smoke:data=data[:min(len(data),3072)]
        if not data:raise ValueError('empty input: '+str(path))
        output.append({'id':entry['id'],'family':entry['family'],'source':entry['path'],
                       'full_sha256':digest(full),'sha256':digest(data),'bytes':len(data),
                       'data':data,'label':'heldout_natural' if split=='heldout' else 'calibration'})
    if not output:raise ValueError('no '+split+' inputs')
    return output


def get_training_samples(calibration:list[dict])->list[bytes]:
    # Data and dictionary samples are exclusively from calibration entries.
    chunks=[]
    for obj in calibration:
        data=obj['data']
        for offset in range(0,len(data),256):
            piece=data[offset:offset+256]
            if len(piece)>16:chunks.append(piece)
    return chunks


def config_integrity(calibration:list[dict])->list[dict]:
    return [{'id':x['id'],'sha256':x['sha256'],'bytes':x['bytes']} for x in calibration]


def calibrate(manifest_path:Path,config_path:Path,smoke=False,allow_missing=False)->dict:
    manifest=get_manifest(manifest_path)
    calibration=read_split(manifest,'calibration',smoke=smoke,permit_missing=allow_missing)
    zstd=Zstd()
    samples=get_training_samples(calibration)
    start=time.perf_counter_ns()
    capacity=1024 if smoke else 4096
    dictionary=zstd.train(samples,capacity)
    dict_ms=(time.perf_counter_ns()-start)/1e6
    chosen={};scores={}
    calibration_start=time.perf_counter_ns()
    for variant in VARIANTS:
        scorecard=[]
        for block_size in ([32,128] if smoke else BLOCK_SIZES):
            cfg=HybridConfig(block_size,variant)
            total_bytes=0;times=[]
            for source in calibration:
                start=time.perf_counter_ns()
                coded=hybrid_encode(source['data'],cfg,zstd)
                times.append((time.perf_counter_ns()-start)/1e6)
                # Round-trip check is part of calibration (not included in encode times).
                assert hybrid_decode(coded,zstd)==source['data']
                total_bytes+=len(coded)
            scorecard.append({'block_bytes':block_size,'sum_encoded_bytes':total_bytes,
                              'sum_encode_ms':round(sum(times),3)})
        best=min(scorecard,key=lambda x:(x['sum_encoded_bytes'],x['sum_encode_ms'],x['block_bytes']))
        chosen[variant]={'block_bytes':best['block_bytes'],'variant':variant}
        scores[variant]=scorecard
    search_ms=(time.perf_counter_ns()-calibration_start)/1e6
    result={'schema':'hypercomplex_sota_benchmark_config_v1','created_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
            'protocol':manifest['protocol'],'smoke':smoke,'max_bytes':3072 if smoke else 0,
            'training_files':config_integrity(calibration),'heldout_files_not_accessed':True,
            'zstd_version':zstd.version,'dictionary_size':len(dictionary),
            'dictionary_sha256':digest(dictionary),'dictionary_base64':base64.b64encode(dictionary).decode(),
            'dictionary_training_ms':round(dict_ms,3), 'calibration_search_ms':round(search_ms,3),
            'block_candidates':[32,128] if smoke else BLOCK_SIZES,
            'selected':chosen,'candidate_scores':scores}
    config_path.parent.mkdir(parents=True,exist_ok=True)
    config_path.write_text(json.dumps(result,indent=2)+'\n')
    return result


def memory_observation(fn:Callable)->dict:
    # Separate replay, so the RSS/tracemalloc instrumentation does not affect throughput timings.
    process=psutil.Process()
    observed=[];halt=threading.Event()
    def watch():
        while not halt.wait(0.001):
            try:observed.append(process.memory_info().rss)
            except psutil.Error:break
    baseline=process.memory_info().rss
    watcher=threading.Thread(target=watch,daemon=True)
    tracemalloc.start()
    watcher.start()
    try:
        fn()
        _,peak=tracemalloc.get_traced_memory()
    finally:
        halt.set();watcher.join(timeout=2)
        tracemalloc.stop()
    return {'rss_before_bytes':baseline,'rss_peak_sampled_bytes':max([baseline,process.memory_info().rss,*observed]),
            'python_alloc_peak_bytes':peak}


def benchmark_one(data:bytes,encode:Callable[[],bytes],decode:Callable[[bytes],bytes],
                  repeats:int=2,include_memory:bool=True)->tuple[dict,bytes]:
    # Warmup excludes compilation, then full call repeated. All byte-for-byte verified.
    coded=encode()
    if decode(coded)!=data:raise AssertionError('codec round-trip failed')
    tenc=[];tdec=[]
    for _ in range(repeats):
        st=time.perf_counter_ns();c=encode();tenc.append((time.perf_counter_ns()-st)/1e6)
        if c!=coded:raise AssertionError('unstable encoded payload')
        st=time.perf_counter_ns();out=decode(c);tdec.append((time.perf_counter_ns()-st)/1e6)
        if out!=data:raise AssertionError('decode mismatch')
    e=statistics.median(tenc);d=statistics.median(tdec)
    memory=memory_observation(lambda:(encode(),decode(coded))) if include_memory else {}
    return ({'encoded_bytes':len(coded),'encode_ms':round(e,4),'decode_ms':round(d,4),
             'encode_MBps':round((len(data)/1e6)/(e/1000),3) if e else None,
             'decode_MBps':round((len(data)/1e6)/(d/1000),3) if d else None,
             'exact':True,'encode_trials_ms':[round(v,4) for v in tenc],
             'decode_trials_ms':[round(v,4) for v in tdec],**memory},coded)


def synthetic_basis_signal(D:int=1024,blocks:int=4,seed:int=20073)->bytes:
    rng=np.random.default_rng(seed)
    raw_bits=[]
    for _ in range(blocks):
        idx=int(rng.integers(0,D))
        y=np.zeros(D,dtype=np.int64);y[idx]=D
        basis=cd_cascade_inverse_checked(y)
        if basis[0]==-1:basis=-basis
        if basis[0]!=1 or np.any(np.abs(basis[1:])!=1):raise AssertionError('not a bipolar basis codeword')
        raw_bits.append(((basis[1:]+1)//2).astype(np.uint8))
    return np.packbits(np.concatenate(raw_bits)).tobytes()


def evaluate(manifest_path:Path,config_path:Path,output_dir:Path,smoke=False,
             allow_missing=False,synthetic=False,repeats:int=2,memory=True,max_bytes=0)->dict:
    manifest=get_manifest(manifest_path)
    cfg=json.loads(config_path.read_text())
    if cfg['schema']!='hypercomplex_sota_benchmark_config_v1':raise ValueError('wrong config schema')
    if bool(cfg['smoke'])!=bool(smoke):raise ValueError('smoke/calibration mode mismatch')
    calibration=read_split(manifest,'calibration',smoke=smoke,permit_missing=allow_missing)
    if cfg['training_files']!=config_integrity(calibration):
        raise ValueError('calibration sources changed since configuration was frozen')
    test=read_split(manifest,'heldout',smoke=smoke,max_bytes=max_bytes,permit_missing=allow_missing)
    if synthetic:
        # One cross-dimension diagnostic and one aligned to the frozen calibration choice.
        frozen_dim=8*cfg['selected']['hypercomplex']['block_bytes']
        for label,D,blocks in [('basis_cross_dimension_1024',1024,4 if smoke else 24),
                               ('basis_frozen_dimension',frozen_dim,4 if smoke else 12)]:
            data=synthetic_basis_signal(D=D,blocks=blocks)
            test.append({'id':label,'family':'synthetic_basis','source':f'in-code generator seed 20073 / D={D}',
                         'full_sha256':digest(data),'sha256':digest(data),'bytes':len(data),
                         'data':data,'label':'synthetic_diagnostic','basis_dimension':D})
    zstd=Zstd()
    dictionary=base64.b64decode(cfg['dictionary_base64'])
    if digest(dictionary)!=cfg['dictionary_sha256']:raise ValueError('dictionary/config mismatch')
    names=STANDARD if not smoke else ['zstd1','zstd3','zstd9','zstd19','brotli5','brotli9','brotli11','lzma6','lzma9','zstd_dict3','zstd_dict9']
    results=[]; files=[]
    output_dir.mkdir(parents=True,exist_ok=True)
    for source in test:
        data=source['data']; entries=[]
        methods=[(name,lambda name=name:compress_named(data,name,zstd,dictionary),
                        lambda out,name=name:decompress_named(out,name,zstd,len(data)),None) for name in names]
        for variant in VARIANTS:
            sel=cfg['selected'][variant]
            hcfg=HybridConfig(**sel)
            stats={}
            methods.append(('hybrid_'+variant,lambda conf=hcfg,s=stats:hybrid_encode(data,conf,zstd,s),
                           lambda out:hybrid_decode(out,zstd),stats))
        # Causal ablation: fix the calibrated hypercomplex block size and change ONLY the transform.
        # Unlike the separately tuned candidates above, these are controlled comparisons.
        common_block=cfg['selected']['hypercomplex']['block_bytes']
        for variant in ('walsh','xor','no_spectral','basis_index'):
            hcfg=HybridConfig(common_block,variant)
            stats={}
            methods.append(('fixed_ablation_'+variant,lambda conf=hcfg,s=stats:hybrid_encode(data,conf,zstd,s),
                           lambda out:hybrid_decode(out,zstd),stats))
        print(f'[{source["label"]}] {source["id"]}: {len(data)} bytes / {len(methods)} codecs',flush=True)
        for name,enc,dec,mode_stats in methods:
            metrics,coded=benchmark_one(data,enc,dec,repeats,include_memory=memory)
            row={'dataset':source['id'],'label':source['label'],'family':source['family'],'source':source['source'],
                 'sha256':source['sha256'],'input_bytes':len(data),'codec':name,
                 'ratio':round(metrics['encoded_bytes']/len(data),6),
                 'mode_counts':dict(mode_stats.get('modes',{})) if mode_stats is not None else {},
                 **metrics}
            results.append(row);entries.append(row)
        # Specialized codebook-index comparator. Undefined on non-family inputs, not an oracle fallback.
        basis_dim=source.get('basis_dimension',1024)
        candidate=basis_codebook_encode(data,basis_dim)
        if candidate is not None:
            bmetrics,_=benchmark_one(data,lambda:basis_codebook_encode(data,basis_dim),basis_codebook_decode,repeats,include_memory=memory)
            brow={'dataset':source['id'],'label':source['label'],'family':source['family'],
                  'source':source['source'],'sha256':source['sha256'],'input_bytes':len(data),
                  'codec':f'basis_codebook_{basis_dim}','ratio':round(bmetrics['encoded_bytes']/len(data),6),
                  'mode_counts':{'basis_codewords':(len(data)*8)//1023},**bmetrics}
            results.append(brow);entries.append(brow)
        baseline=min((r for r in entries if not r['codec'].startswith(('hybrid_','fixed_ablation_','basis_codebook_'))),key=lambda x:x['encoded_bytes'])
        hyper=next(r for r in entries if r['codec']=='hybrid_hypercomplex')
        files.append({'id':source['id'],'label':source['label'],'input_bytes':len(data),
                      'best_baseline_oracle':baseline['codec'],'best_baseline_bytes':baseline['encoded_bytes'],
                      'frozen_hypercomplex_bytes':hyper['encoded_bytes'],
                      'change_vs_oracle_pct':round(100*(hyper['encoded_bytes']/baseline['encoded_bytes']-1),3),
                      'source_sha256':source['sha256'],
                      'basis_codebook_bytes':next((r['encoded_bytes'] for r in entries if r['codec'].startswith('basis_codebook_')),None)})
    environment={'python':platform.python_version(),'platform':platform.platform(),'CPU':platform.processor(),
                 'zstd':zstd.version,'brotli':getattr(__import__('brotli'),'__version__','unknown'),
                 'numpy':np.__version__,'memory_measurement':'sampled process RSS + peak Python allocations; not a hardware memory ceiling'}
    output={'schema':'hypercomplex_sota_benchmark_results_v1','protocol':cfg['protocol'],
            'calibration_config_sha256':digest(config_path.read_bytes()),
            'calibration_training_files':cfg['training_files'],
            'dictionary_bytes_charged_per_file':len(dictionary),'dictionary_training_ms':cfg['dictionary_training_ms'],
            'calibration_search_ms':cfg['calibration_search_ms'],
            'selected_configurations':cfg['selected'], 'environment':environment,
            'smoke':smoke,'fairness':'frozen hyperparameters on heldout; per-file conventional minimum only oracle diagnostic',
            'datasets':files,'measurements':results}
    (output_dir/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    with (output_dir/'measurements.csv').open('w',newline='') as f:
        fields=['dataset','label','family','sha256','input_bytes','codec','encoded_bytes','ratio','encode_ms','decode_ms',
                'encode_MBps','decode_MBps','rss_before_bytes','rss_peak_sampled_bytes','python_alloc_peak_bytes','exact','mode_counts']
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader()
        for row in results:
            w.writerow({k:(json.dumps(row[k]) if k=='mode_counts' else row.get(k)) for k in fields})
    lines=['# Held-out compression benchmark',
           '',f'Calibration: {len(calibration)} files; held-out: {sum(v["label"]=="heldout_natural" for v in files)} files.',
           '', 'Configurations frozen on calibration files before evaluating held-out files. The best conventional compressor shown is a **post-hoc per-file oracle**, not an operational method.',
           '', '| Dataset | Input (B) | Best conventional (B) | Frozen hypercomplex (B) | Change |',
           '|---|---:|---:|---:|---:|']
    for f in files:
        lines.append(f'| {f["id"]} ({f["label"]}) | {f["input_bytes"]} | {f["best_baseline_bytes"]} ({f["best_baseline_oracle"]}) | {f["frozen_hypercomplex_bytes"]} | {f["change_vs_oracle_pct"]:+.2f}% |')
    lines += ['', f'Calibration-search CPU time: {cfg["calibration_search_ms"]} ms; dictionary train CPU time: {cfg["dictionary_training_ms"]} ms.',
              '', 'Detailed per-method compression/decompression throughput, peak sampled RSS, and Python allocation measurements are in `measurements.csv`.',
              'Zstandard dictionary frames include the entire calibration dictionary once **per encoded file**; no hidden side channel.',
              'Synthetic basis-family examples are explicitly labeled and excluded from natural-data conclusions.',
              'All held-out source SHA-256 hashes, file paths, library versions, and bit-exact verifications are in results.json.',
              'Limitations: one machine, small suite, per-method RSS sampling misses short allocation spikes, and no GPU kernels.']
    (output_dir/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    return output


def main():
    p=argparse.ArgumentParser(description='Five-part held-out SOTA hypercomplex compression benchmark')
    p.add_argument('phase',choices=('calibrate','evaluate','all'))
    p.add_argument('--manifest',type=Path,default=DEFAULT_MANIFEST)
    p.add_argument('--config',type=Path,default=Path('/tmp/hypercomplex_sota_calibration.json'))
    p.add_argument('--output',type=Path,default=Path('/tmp/hypercomplex_sota_results'))
    p.add_argument('--smoke',action='store_true',help='small inputs, bounded calibration candidates; still tests all codecs')
    p.add_argument('--allow-missing',action='store_true',help='skip absent optional repository files (local snapshots only)')
    p.add_argument('--max-bytes',type=int,default=0,help='explicit heldout byte truncation, 0 means whole files')
    p.add_argument('--synthetic',action='store_true',help='evaluate constructed basis signals separately from heldout')
    p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--no-memory',action='store_true',help='skip expensive memory replay; size/speed still measured')
    args=p.parse_args()
    if args.repeats<1:raise ValueError('repeats must be positive')
    if args.phase in ('calibrate','all'):
        c=calibrate(args.manifest,args.config,args.smoke,args.allow_missing)
        print('Frozen calibration:',args.config,'training files',len(c['training_files']),
              'dictionary',c['dictionary_size'],'bytes, search',c['calibration_search_ms'],'ms')
    if args.phase in ('evaluate','all'):
        out=evaluate(args.manifest,args.config,args.output,args.smoke,args.allow_missing,
                     args.synthetic,args.repeats,not args.no_memory,args.max_bytes)
        print('Wrote',args.output,'with',len(out['datasets']),'datasets and',len(out['measurements']),'rows')


if __name__=='__main__':main()
