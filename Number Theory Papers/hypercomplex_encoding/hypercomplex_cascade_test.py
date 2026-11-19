#!/usr/bin/env python3
"""Full Cayley-Dickson imaginary-dimension *mixing* transform and codec.

Generalizes Viole phi by multiplying successively by 1+e_(2**r).
The first stage is exactly phi on every coordinate pair; subsequent
stages are actual hypercomplex multiplications, not a tensor definition.
Although the Cayley-Dickson algebra is nonassociative from dimension 8,
the sequence of *linear maps* is invertible because each L_(e_j)^2=-I.
"""
from __future__ import annotations
import csv, json, os, struct, statistics, time, zlib, lzma, platform
from pathlib import Path
import numpy as np
from numba import njit
import brotli
from hypercomplex_highdim_test import (ROOT, SEED, cd_sign, sign_cd_ref, corpus, repeated_time)

MAGIC=b'HCSP'
HEADER='<4sIQI'
HSZ=struct.calcsize(HEADER)

@njit(cache=True)
def cd_cascade(coeff):
    n=len(coeff)
    x=np.empty(n,dtype=np.int64)
    for i in range(n): x[i]=coeff[i]
    stride=1
    while stride<n:
        y=np.empty(n,dtype=np.int64)
        for i in range(n):
            j=i^stride
            y[i]=x[i]+cd_sign(stride,j,n)*x[j]
        x=y
        stride*=2
    return x

@njit(cache=True)
def cd_cascade_inverse(transformed):
    n=len(transformed)
    x=transformed.copy()
    stride=n//2
    exact=True
    while stride>=1:
        y=np.empty(n,dtype=np.int64)
        for i in range(n):
            j=i^stride
            v=x[i]-cd_sign(stride,j,n)*x[j]
            if v % 2 != 0:
                exact=False
            y[i]=v//2
        x=y
        stride//=2
    return x,exact


def cd_cascade_inverse_checked(y):
    x,exact=cd_cascade_inverse(y)
    if not exact: raise ValueError('nonintegral inverse: modified or invalid spectrum')
    return x


def varint(v):
    assert v>=0
    b=bytearray()
    while v>=128:
        b.append((v&127)|128)
        v>>=7
    b.append(v)
    return bytes(b)


def unvarint(data,pos):
    x=0;shift=0
    for _ in range(11):
        if pos>=len(data): raise ValueError('truncated varint')
        v=data[pos];pos+=1
        x|=(v&127)<<shift
        if v<128:return x,pos
        shift+=7
    raise ValueError('varint too long')


def zigzag(x):return x*2 if x>=0 else -2*x-1

def unzigzag(x):return x//2 if x%2==0 else -(x+1)//2


def spectrum_counts(seed_rng):
    results=[]
    for n in (2,4,8,16,32,64,256,1024,4096,65536,1048576):
        vecs={}
        if n<=4096:
            vecs['ones']=np.ones(n,dtype=np.int8)
            vecs['alternating']=np.tile([0,1],n//2).astype(np.int8)
            vecs['repeat1011']=np.tile([1,0,1,1],(n+3)//4)[:n].astype(np.int8)
        vecs['random']=seed_rng.integers(0,2,size=n,dtype=np.int8)
        vecs['single_coefficient']=np.eye(1,n,0,dtype=np.int8)[0]
        for name, x in vecs.items():
            y,enc_ms=repeated_time(lambda:cd_cascade(x), 2)
            recovered,dec_ms=repeated_time(lambda:cd_cascade_inverse_checked(y), 2)
            assert np.array_equal(recovered,x)
            norm_before=int(np.dot(x.astype(np.int64),x.astype(np.int64)))
            norm_after=int(np.dot(y,y))
            assert norm_after==n*norm_before,(n,name,norm_before,norm_after)
            results.append({'dimension':n,'signal':name,'nonzero_input':int(np.count_nonzero(x)),
                            'nonzero_spectrum':int(np.count_nonzero(y)),
                            'spectrum_maxabs':int(np.max(np.abs(y))),
                            'spectrum_l2_squared':norm_after,
                            'encode_ms_median':enc_ms,'decode_ms_median':dec_ms,
                            'exact_inverse':True})
    return results


def check_imaginary_mixing(seed_rng):
    # Verify actual paper phi (2D) and quaternion psi (4D) via left (1+i) action.
    for a in range(-4,5):
        for b in range(-4,5):
            v=np.array([a,b],dtype=np.int64)
            stage_v=np.array([v[i]+cd_sign(1,i^1,2)*v[i^1] for i in range(2)],dtype=np.int64)
            assert np.array_equal(stage_v,[a-b,a+b])
    for _ in range(300):
        a,b,c,d=seed_rng.integers(-20,21,size=4)
        v=np.array([a,b,c,d],dtype=np.int64)
        stage_v=np.array([v[i]+cd_sign(1,i^1,4)*v[i^1] for i in range(4)],dtype=np.int64)
        assert np.array_equal(stage_v,[a-b,a+b,c-d,c+d])
    # Validate higher basis operators L_j^2=-I (for every basis j up to 256)
    all_basis=0
    for n in (2,4,8,16,32,64,128,256):
        for j in range(1,n):
            for i in range(n):
                assert cd_sign(j,i,n)*cd_sign(j,i^j,n)==-1,(n,j,i)
                all_basis+=1
    # Random checks at a million dimensions, including e_j^2=-1 and anticommutativity.
    huge=1<<20
    for _ in range(4000):
        i=int(seed_rng.integers(0,huge))
        j=int(seed_rng.integers(1,huge))
        assert cd_sign(j,i,huge)==sign_cd_ref(j,i,huge)
        assert cd_sign(j,i,huge)*cd_sign(j,i^j,huge)==-1
    # Full matrix exact orthogonality and coefficients ±1, not just norm on random signals.
    gram={}
    for n in (2,4,8,16,32,64):
        M=np.empty((n,n),dtype=np.int64)
        for i in range(n):
            ei=np.zeros(n,dtype=np.int8); ei[i]=1
            M[:,i]=cd_cascade(ei)
        assert np.all(np.abs(M)==1),('not Hadamard',n)
        assert np.array_equal(M.T@M,n*np.eye(n,dtype=np.int64))
        gram[n]={'all_entries_plusminus_one':True,'gram_nI':True,
                 'row_all_ones_exists':bool(np.any(np.all(M==1,axis=1) | np.all(M==-1,axis=1))),
                 'distinct_from_tensor_map':None}
    return {'paper_phi_and_psi_verified':True,'basis_operator_pairs':all_basis,
            'million_dim_basis_checks':4000,'small_gram_checks':gram}


def encode_sparse_cascade(payload,n,bipolar=False):
    assert n>=2 and n&(n-1)==0
    bits=np.unpackbits(np.frombuffer(payload,dtype=np.uint8))
    per=n-1;blocks=(len(bits)+per-1)//per
    out=bytearray(struct.pack(HEADER,b'HCSB' if bipolar else MAGIC,n,len(payload),zlib.crc32(payload)))
    specs=0;fills=0;nonzero=[]
    for block in range(blocks):
        chunk=bits[block*per:(block+1)*per]
        vec=np.ones(n,dtype=np.int8) if bipolar else np.zeros(n,dtype=np.int8)
        if bipolar:
            vec[1:]= -1
            vec[1:len(chunk)+1]=2*chunk.astype(np.int8)-1
        else:
            vec[1:len(chunk)+1]=chunk
        if np.all(chunk==0):
            out.append(2);fills+=1;continue
        if np.all(chunk==1) and len(chunk)==per:
            out.append(3);fills+=1;continue
        raw_bits=np.zeros(per,dtype=np.uint8)
        raw_bits[:len(chunk)]=chunk
        raw=np.packbits(raw_bits).tobytes()
        y=cd_cascade(vec)
        idx=np.flatnonzero(y)
        nonzero.append(len(idx))
        # Build sparse only if below reasonable best-case threshold; then true byte count chooses.
        if len(idx)<len(raw)//2:
            sp=bytearray(varint(len(idx)))
            prev=0
            for k,coordinate in enumerate(idx):
                coordinate=int(coordinate)
                sp.extend(varint(coordinate-prev));sp.extend(varint(zigzag(int(y[coordinate]))))
                prev=coordinate
            if len(sp)<len(raw):
                out.append(1);out.extend(sp);specs+=1;continue
        out.append(0);out.extend(raw)
    return bytes(out),{'spectral':specs,'fill':fills,'raw':blocks-specs-fills,
                       'blocks':blocks,'avg_spectrum_nonzero':round(statistics.mean(nonzero),1) if nonzero else 0}


def decode_sparse_cascade(blob):
    if len(blob)<HSZ:raise ValueError('incomplete header')
    magic,n,rawbytes,crc=struct.unpack_from(HEADER,blob,0)
    if magic not in (MAGIC,b'HCSB') or n<2 or n&(n-1):raise ValueError('not hypercomplex spectral codec')
    bipolar = magic==b'HCSB'
    per=n-1;blocks=(rawbytes*8+per-1)//per
    bits_out=np.zeros(blocks*per,dtype=np.uint8)
    pos=HSZ; per_raw=(per+7)//8
    for k in range(blocks):
        if pos>=len(blob):raise ValueError('truncated input')
        mode=blob[pos];pos+=1
        if mode==0:
            if pos+per_raw>len(blob):raise ValueError('truncated raw block')
            bits_out[k*per:(k+1)*per]=np.unpackbits(np.frombuffer(blob[pos:pos+per_raw],dtype=np.uint8))[:per]
            pos+=per_raw
        elif mode==1:
            count,pos=unvarint(blob,pos)
            if count>n:raise ValueError('invalid sparse count')
            spectrum=np.zeros(n,dtype=np.int64);idx=0
            for i in range(count):
                delta,pos=unvarint(blob,pos)
                v,pos=unvarint(blob,pos)
                idx+=delta
                if idx>=n or (i>0 and delta==0):raise ValueError('invalid sparse index')
                spectrum[idx]=unzigzag(v)
            coeff=cd_cascade_inverse_checked(spectrum)
            if bipolar:
                if coeff[0]!=1 or np.any((coeff[1:]!=-1)&(coeff[1:]!=1)):
                    raise ValueError('invalid bipolar spectrum')
                bits_out[k*per:(k+1)*per]=((coeff[1:]+1)//2).astype(np.uint8)
            else:
                if coeff[0]!=0 or np.any((coeff!=0)&(coeff!=1)):
                    raise ValueError('invalid binary spectrum')
                bits_out[k*per:(k+1)*per]=coeff[1:].astype(np.uint8)
        elif mode==2:
            pass
        elif mode==3:
            bits_out[k*per:(k+1)*per]=1
        else:raise ValueError('bad mode')
    if pos!=len(blob):raise ValueError('unexpected bytes')
    original=np.packbits(bits_out[:rawbytes*8]).tobytes()
    if zlib.crc32(original)!=crc:raise ValueError('CRC mismatch')
    return original


def spectral_codec_comparison(rng):
    datasets=corpus(rng)
    # Synthetic signals aligned to the *actual* Cayley-Dickson mixing basis.
    for n,rows in ((1024,512),(4096,128)):
        generated=[]
        for i in range(rows):
            ei=np.zeros(n,dtype=np.int8)
            ei[int(rng.integers(0,n))]=1
            # Columns of the transpose M are ±1. Encode 0/1 from those signs.
            # Strict imaginary input: drop real coefficient then pack n-1 per block.
            # This construction can be highly compressible using a basis-specific codec.
            # Derive ±1 pattern using inverse transform of n*e_j: M^{-1}(n e_j)=M^T e_j.
            ei=int(n)*ei.astype(np.int64)
            raw=cd_cascade_inverse_checked(ei)
            assert np.all(np.abs(raw)==1)
            if raw[0] == -1:raw=-raw
            assert raw[0]==1
            bitvector=((raw+1)//2).astype(np.uint8)
            generated.append(bitvector[1:])
        merged=np.concatenate(generated)
        # Our block n-1 doesn't align exactly with byte boundary; packing pads end only.
        datasets[f'CD_basis_binary_{n}']=np.packbits(merged).tobytes()
    results=[];details=[]
    for name,payload in datasets.items():
        z=zlib.compress(payload,9)
        bro=brotli.compress(payload,quality=9)
        lz=lzma.compress(payload,preset=6)
        baseline=min([('zlib9',len(z)),('brotli9',len(bro)),('lzma6',len(lz))],key=lambda x:x[1])
        per=[]
        for n in (256,1024,4096,65536):
            if n==65536 and len(payload)<8192:continue
            for bipolar in (False,True):
                (encoded,info), tenc=repeated_time(lambda:encode_sparse_cascade(payload,n,bipolar),2)
                recovered,tdec=repeated_time(lambda:decode_sparse_cascade(encoded),2)
                assert recovered==payload,(name,n,bipolar)
                cz=len(zlib.compress(encoded,9))
                cb=len(brotli.compress(encoded,quality=9))
                compressed=min([('zlib9',cz),('brotli9',cb)],key=lambda x:x[1])
                per.append({'dimension':n,'bipolar':bipolar,'raw_file_bytes':len(encoded),
                            'zlib9_bytes':cz,'brotli9_bytes':cb,'compressed_bytes':compressed[1],
                            'compressor':compressed[0],'encode_ms_median':tenc,'decode_ms_median':tdec,
                            'mode_counts':info,'exact':True})
        best=min(per,key=lambda x:x['compressed_bytes'])
        results.append({'dataset':name,'input_bytes':len(payload),
                        'zlib9_bytes':len(z),'brotli9_bytes':len(bro),'lzma6_bytes':len(lz),
                        'best_conventional_name':baseline[0], 'best_conventional_bytes':baseline[1],
                        'best_cd_cascade_dim':best['dimension'],'bipolar':best['bipolar'],'best_cd_cascade_bytes':best['compressed_bytes'],
                        'best_cd_mode_counts':best['mode_counts'],
                        'pct_larger_than_baseline':round(100*(best['compressed_bytes']/baseline[1]-1),3),
                        'exact':True})
        details.append({'dataset':name,'variants':per,'baseline':baseline})
        print(f"[cascade codec] {name:22s} baseline={baseline[1]:7d} cascade={best['compressed_bytes']:7d} dim={best['dimension']:,} modes={best['mode_counts']['spectral']}",flush=True)
    return results,details


def run():
    rng=np.random.default_rng(SEED+7)
    t0=time.perf_counter()
    cd_cascade(np.array([1,0],dtype=np.int8));cd_cascade_inverse(np.array([1,1],dtype=np.int64))
    codec,details=spectral_codec_comparison(np.random.default_rng(SEED+701))
    algebra=check_imaginary_mixing(rng)
    print('[cascade] phi psi verified; gram sizes',list(algebra['small_gram_checks']), flush=True)
    spectrum=spectrum_counts(rng)
    for x in spectrum:
        if x['dimension'] in (1024,4096,1048576):
            print('[spectrum]',x['dimension'],x['signal'],'input',x['nonzero_input'],'output',x['nonzero_spectrum'],'ms',x['encode_ms_median'], flush=True)
    result={'seed':SEED+7,'codec_seed':SEED+701,'algebra':algebra,'spectrum_scaling':spectrum,
            'compression':codec,'compression_details':details,
            'runtime_seconds':round(time.perf_counter()-t0,3)}
    (ROOT/'hypercomplex_cascade_results.json').write_text(json.dumps(result,indent=2)+'\n')
    for name,rows in [('hypercomplex_cascade_spectrum.csv',spectrum),('hypercomplex_cascade_compression.csv',codec)]:
        with (ROOT/name).open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
    print('[cascade] complete in',result['runtime_seconds'],'s')
    return result

if __name__=='__main__':
    import argparse, subprocess, sys
    arg=argparse.ArgumentParser(description='Test actual Cayley-Dickson high-dimensional cascade encoding')
    arg.add_argument('--phase',choices=('all','codec','spectrum'),default='all')
    args=arg.parse_args()
    if args.phase=='codec':
        rows,details=spectral_codec_comparison(np.random.default_rng(SEED+701))
        output={'seed':SEED+701,'compression':rows,'compression_details':details}
        (ROOT/'hypercomplex_cascade_bipolar_results.json').write_text(json.dumps(output,indent=2)+'\n')
        with (ROOT/'hypercomplex_cascade_compression_full.csv').open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
    elif args.phase=='spectrum':
        rng=np.random.default_rng(SEED+7)
        algebra=check_imaginary_mixing(rng)
        rows=spectrum_counts(rng)
        out={'seed':SEED+7,'algebra':algebra,'spectrum_scaling':rows}
        (ROOT/'hypercomplex_cascade_spectrum_only.json').write_text(json.dumps(out,indent=2)+'\n')
        with (ROOT/'hypercomplex_cascade_spectrum.csv').open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
        print('[cascade spectrum] validated',len(rows),'examples, max dimension',max(x['dimension'] for x in rows),flush=True)
    else:
        for phase in ('codec','spectrum'):
            subprocess.run([sys.executable,str(Path(__file__).resolve()),'--phase',phase],check=True)
        print('[cascade] complete in isolated phases',flush=True)
