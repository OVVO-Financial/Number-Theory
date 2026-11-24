#!/usr/bin/env python3
"""High-dimensional *imaginary-algebra* bit encoding test.

One bit per independent imaginary basis coefficient in Cayley-Dickson
algebras C,H,O,sedenions,..., through 1,048,576 total real dimensions.
Distinguishes an actual Cayley-Dickson multiplication from the Walsh
coordinate transform. Also tests associative Clifford algebra vector encoding.

Requires: numpy, numba, brotli, Pillow for the page bitmap case.
All data deterministic; output JSON, CSV, report created by run().
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import struct
import time
import zlib
import lzma
import brotli
import numpy as np
from numba import njit
from PIL import Image

ROOT = Path(__file__).resolve().parent
SEED = 20261118
MAGIC = b'HCD1'
HEADER = '<4sIIQI'
H_SIZE = struct.calcsize(HEADER)


def sign_cd_ref(p: int, q: int, n: int) -> int:
    """Basis product: e_p * e_q = sign * e_(p XOR q).

    Cayley-Dickson convention: (a,b)(c,d)=(ac-conj(d)b, da+b conj(c)).
    i*j=k, j*i=-k in dimension four. For n>=16 this algebra is
    nonassociative and has zero divisors; that is part of the test.
    """
    if n == 1:
        assert p == q == 0
        return 1
    half = n // 2
    if p < half and q < half:
        return sign_cd_ref(p, q, half)
    if p < half:
        return sign_cd_ref(q-half, p, half)
    if q < half:
        return (-1 if q else 1) * sign_cd_ref(p-half, q, half)
    return (1 if q > half else -1) * sign_cd_ref(q-half, p-half, half)


@njit(cache=True)
def cd_sign(p: int, q: int, n: int) -> int:
    sg = 1
    while n > 1:
        h = n // 2
        if p < h:
            if q >= h:
                tmp = p
                p = q - h
                q = tmp
        else:
            if q < h:
                if q != 0:
                    sg = -sg
                p -= h
            else:
                if q == h:
                    sg = -sg
                tmp = p - h
                p = q - h
                q = tmp
        n = h
    return sg


@njit(cache=True)
def cd_product_dense(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = len(a)
    result = np.zeros(n, dtype=np.int64)
    for i in range(n):
        ai = int(a[i])
        if not ai:
            continue
        for j in range(n):
            bj = int(b[j])
            if bj:
                result[i ^ j] += ai * bj * cd_sign(i, j, n)
    return result


@njit(cache=True)
def cd_product_sparse_indices(ia: np.ndarray, va: np.ndarray,
                              ib: np.ndarray, vb: np.ndarray, n: int) -> np.ndarray:
    result = np.zeros(n, dtype=np.int64)
    for ai in range(len(ia)):
        for bj in range(len(ib)):
            i = ia[ai]
            j = ib[bj]
            result[i ^ j] += int(va[ai]) * int(vb[bj]) * cd_sign(i, j, n)
    return result


@njit(cache=True)
def cd_unit_product(bits: np.ndarray, key: int) -> np.ndarray:
    """Actual algebra product e_key * sum_i bits[i] e_i.

    The output has signed coefficients. Values in {-1,0,1} for input bits.
    """
    n = len(bits)
    result = np.zeros(n, dtype=np.int8)
    for i in range(n):
        if bits[i]:
            result[i ^ key] = np.int8(cd_sign(key, i, n))
    return result


def repeated_time(fn, min_trials=3):
    timings = []
    last = None
    for _ in range(min_trials):
        t0 = time.perf_counter_ns()
        last = fn()
        timings.append((time.perf_counter_ns() - t0)/1e6)
    return last, round(statistics.median(timings), 4)


def tests_algebra(seed_rng):
    base_checks = 0
    # Independent recursive vs JIT definition for ALL basis pairs <= 128
    for n in (2,4,8,16,32,64,128):
        for i in range(n):
            assert cd_sign(i, i, n) == (-1 if i else 1)
            for j in range(n):
                sg = cd_sign(i, j, n)
                assert sg == sign_cd_ref(i, j, n), (i,j,n,sg)
                if i and j and i!=j:
                    assert sg == -cd_sign(j, i, n)
                base_checks += 1
    assert cd_sign(1,2,4)==1 and cd_sign(2,1,4)==-1
    assert cd_sign(1,3,4)==-1 and cd_sign(3,1,4)==1
    assert cd_sign(2,3,4)==1 and cd_sign(3,2,4)==-1
    # octonion associator; absent for quaternions
    assoc = {}
    for n in (4,8,16):
        first = None
        for i in range(1,n):
            for j in range(1,n):
                for k in range(1,n):
                    sleft = cd_sign(i,j,n)*cd_sign(i^j,k,n)
                    sright = cd_sign(j,k,n)*cd_sign(i,j^k,n)
                    if sleft != sright:
                        first = [i,j,k,sleft,sright,i^j^k]
                        break
                if first: break
            if first: break
        assoc[n] = first
    assert assoc[4] is None and assoc[8] is not None

    # Zero divisor by explicit two-sparse combinations. Search entire n=16
    zero_pair = None
    for n in (8,16):
        found = None
        for i in range(1,n):
            for j in range(i+1,n):
                for k in range(1,n):
                    for l in range(k+1,n):
                        for s in (-1,1):
                            terms = {}
                            for u,v,w in ((i,k,1),(i,l,s),(j,k,1),(j,l,s)):
                                idx=u^v
                                terms[idx]=terms.get(idx,0)+w*cd_sign(u,v,n)
                            if all(v==0 for v in terms.values()):
                                found=[n,[i,j],[k,l],s]
                                break
                        if found: break
                    if found: break
                if found: break
            if found: break
        if n==8:
            assert found is None, 'Octonion unexpected zero divisor'
        else:
            assert found is not None, 'Sedenion zero divisor not found'
            zero_pair = found
    n, left, right, s = zero_pair
    a=np.zeros(n,dtype=np.int64); b=np.zeros(n,dtype=np.int64)
    a[left]=1; b[right[0]]=1; b[right[1]]=s
    assert not np.any(cd_product_dense(a,b))
    assert np.any(a) and np.any(b)

    # Norm composition: exact for R,C,H,O; counterexample in sedenions
    norm_trials = []
    for n in (2,4,8,16,32,64,128):
        failures = 0
        ratios=[]
        trials=100 if n<=64 else 40
        example=None
        for t in range(trials):
            a=seed_rng.integers(-1,2,size=n,dtype=np.int8)
            b=seed_rng.integers(-1,2,size=n,dtype=np.int8)
            p=cd_product_dense(a,b)
            na=int(a@a); nb=int(b@b)
            np_=int(p@p)
            if na*nb != np_:
                failures+=1
                if example is None: example={'Na':na,'Nb':nb,'Nab':np_}
            ratios.append(np_/(na*nb) if na*nb else 1.)
        assert failures==0 if n<=8 else True
        norm_trials.append({'dimension':n,'trials':trials,'norm_multiplicativity_failures':failures,
                            'avg_ratio':round(statistics.mean(ratios),5),'first_failure':example})
    return {'basis_pairs_checked':base_checks, 'first_associativity_failure':assoc,
            'sedenion_zero_divisor':{'a':f'e{left[0]}+e{left[1]}',
                                     'b':f'e{right[0]}' + ('+' if s==1 else '-')+f'e{right[1]}',
                                     'product':'zero'},'norm_trials':norm_trials}


def exact_dim_scaling(rng):
    results=[]
    # Note: D dimensions, 1 real + D-1 independently addressable imaginary coefficients.
    for n in (2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576):
        bits = np.zeros(n,dtype=np.uint8)
        bits[1:] = rng.integers(0,2,size=n-1,dtype=np.uint8)
        key = min(7,n-1)
        signed, t_unit = repeated_time(lambda: cd_unit_product(bits,key),3)
        # The product e_key*b is a signed coordinate permutation.
        assert np.array_equal(np.abs(signed), bits[np.arange(n,dtype=np.int64)^key])
        recovered = np.abs(signed)[np.arange(n,dtype=np.int64)^key]
        assert np.array_equal(recovered,bits)
        assert signed.min()>=-1 and signed.max()<=1
        raw_packed=np.packbits(bits[1:]).tobytes()
        packed_dec=np.unpackbits(np.frombuffer(raw_packed,dtype=np.uint8))[:n-1]
        assert np.array_equal(packed_dec,bits[1:])
        # Low-level dimension scaling, exclude Python loop allocation from decode metric.
        _, t_de = repeated_time(lambda: np.abs(signed)[np.arange(n,dtype=np.int64)^key],3)
        results.append({
            'dimension_D':n,'imaginary_axes':n-1,'message_bits':n-1,
            'raw_packed_bytes':len(raw_packed),
            'dense_int8_coefficient_bytes':int(bits.nbytes),
            'dense_float64_coefficient_bytes':int(bits.size*8),
            'unit_product_signed_int8_bytes':int(signed.nbytes),
            'nonzero_coefficients':int(np.count_nonzero(bits)),
            'unit_product_nonzero_coefficients':int(np.count_nonzero(signed)),
            'unit_product_ms_median':t_unit, 'unit_inverse_ms_median':t_de,
            'exact_roundtrip':True
        })
    return results


def dense_product_scaling(rng):
    results=[]
    for n in (2,4,8,16,32,64,128,256,512,1024,2048,4096):
        a=rng.integers(0,2,size=n,dtype=np.int8)
        b=rng.integers(0,2,size=n,dtype=np.int8)
        product,ms=repeated_time(lambda: cd_product_dense(a,b),3 if n<=1024 else 2)
        assert len(product)==n
        assert np.all(np.abs(product) <= n), 'Bounds on product coefficient violated'
        nnza=int(np.count_nonzero(a)); nnzb=int(np.count_nonzero(b))
        results.append({'dimension':n,'nnz_left':nnza,'nnz_right':nnzb,
                        'pair_multiply_contributions':nnza*nnzb,
                        'nnz_output':int(np.count_nonzero(product)),
                        'output_absmax':int(np.max(np.abs(product))),
                        'negative_output_coefficients':int(np.sum(product<0)),
                        'output_int16_bytes':int(product.size*2),
                        'cd_product_ms_median':ms,
                        'norm_ratio':round(float(product@product)/(nnza*nnzb),5) if nnza*nnzb else None})
    # Larger dimensions via sparse input, an explicit scaling distinction
    sparse=[]
    for n in (4096,65536,262144,1048576):
        count=64
        ia=np.sort(rng.choice(n,size=count,replace=False).astype(np.int64))
        ib=np.sort(rng.choice(n,size=count,replace=False).astype(np.int64))
        va=np.ones(count,dtype=np.int64)
        vb=np.ones(count,dtype=np.int64)
        prod,ms=repeated_time(lambda: cd_product_sparse_indices(ia,va,ib,vb,n),3)
        sparse.append({'dimension':n,'nnz_left':count,'nnz_right':count,
                       'pairs':count*count,'output_nonzero':int(np.count_nonzero(prod)),
                       'output_nnz_rate_percent':round(100*np.count_nonzero(prod)/n,4),
                       'sparse_product_ms_median':ms,'exact_integer':bool(np.all(np.equal(prod,prod.astype(np.int64))))})
    return results, sparse


def clifford_vector_product(a,b):
    """Associative Clifford Cl_{0,n} product of two grade-1 vectors.

    Scalar: -(a.b); grade-2: (a_i b_j - a_j b_i) e_i e_j.
    Does not allocate the full 2^n-dimensional algebra.
    """
    n = len(a)
    scalar = -int(np.dot(a.astype(np.int64),b.astype(np.int64)))
    A=a.astype(np.int16)
    B=b.astype(np.int16)
    wedge=A[:,None]*B[None,:]-B[:,None]*A[None,:]
    up=wedge[np.triu_indices(n,1)]
    return scalar,up


def clifford_test(rng):
    results=[]
    for n in (4,8,16,32,64,128,256,512,1024):
        a = rng.integers(0,2,size=n,dtype=np.uint8)
        b = rng.integers(0,2,size=n,dtype=np.uint8)
        (scalar,wedge),ms=repeated_time(lambda:clifford_vector_product(a,b),3)
        assert scalar == -int(a.astype(np.int64)@b.astype(np.int64))
        # Multiplying a vector by itself is pure scalar (all bivectors 0).
        same_scalar, same_wedge=clifford_vector_product(a,a)
        assert same_scalar == -int(a.astype(np.int64)@a.astype(np.int64)), (n,same_scalar)
        assert not np.any(same_wedge)
        results.append({'generators_bits':n,'full_clifford_dimension_exact':'2^'+str(n),
                        'vector_coefficients':n,'bivector_possible':n*(n-1)//2,
                        'nonzero_bivectors':int(np.count_nonzero(wedge)),
                        'scalar_coefficient':scalar,'vector_product_ms_median':ms,
                        'array_bivector_i16_bytes':int(wedge.nbytes),
                        'roundtrip_source_vectors':True})
    return results


def corpus(rng):
    files={
        'paper_text':ROOT/'realpair_text.txt',
        'python_source':ROOT/'real_pair_factor_test.py',
        'paper_PDF':ROOT.parent/'Real_Pair_Representation_of_Complex_Numbers_and_Quaternions.pdf',
    }
    out={}
    for k,p in files.items():
        if p.exists():out[k]=p.read_bytes()
    rendered=ROOT/'realpair_renders/page-3.png'
    if rendered.exists():
        im=Image.open(rendered).convert('L').resize((384,496),Image.Resampling.LANCZOS)
        bw=(np.asarray(im)<150).astype(np.uint8)
        out['document_bitmap']=np.packbits(bw,axis=1).tobytes()
    # add held-out, naturally structured but not purpose-built for CD
    out['repeated_pattern']=bytes([0xB5])*65536
    out['random_control']=rng.integers(0,256,size=65536,dtype=np.uint8).tobytes()
    out['sparse_binary']=bytes(bytearray(8192-1)+bytearray([1])) # highly sparse, already trivial for RLE
    return out


def encode_unit(payload, n, key):
    """Lossless codec using e_key multiplication, sign-normalized to bits.

    1 scalar and n-1 imaginary axes, each holding one binary coefficient.
    Crucial: sign stripping is reversible because sign depends on key and axis.
    """
    assert n>=2 and n&(n-1)==0 and 1<=key<n
    raw_bits=np.unpackbits(np.frombuffer(payload,dtype=np.uint8))
    per=n-1
    blocks=(len(raw_bits)+per-1)//per
    transformed=np.zeros(blocks*n,dtype=np.uint8)
    for b in range(blocks):
        chunk=raw_bits[b*per:(b+1)*per]
        coeff=np.zeros(n,dtype=np.uint8)
        coeff[1:1+len(chunk)]=chunk
        signed=cd_unit_product(coeff,key)
        transformed[b*n:(b+1)*n]=np.abs(signed).astype(np.uint8)
    packed=np.packbits(transformed).tobytes()
    return struct.pack(HEADER,MAGIC,n,key,len(payload),zlib.crc32(payload))+packed


def decode_unit(serialized):
    if len(serialized)<H_SIZE:raise ValueError('missing header')
    magic,n,key,original_bytes,crc=struct.unpack(HEADER,serialized[:H_SIZE])
    if magic!=MAGIC or n<2 or (n&(n-1)) or key<=0 or key>=n:
        raise ValueError('bad codec header')
    per=n-1
    blocks=(original_bytes*8+per-1)//per
    if len(serialized[H_SIZE:])!=(blocks*n+7)//8:
        raise ValueError('invalid payload length')
    enc=np.unpackbits(np.frombuffer(serialized[H_SIZE:],dtype=np.uint8))[:blocks*n]
    output=np.zeros(blocks*per,dtype=np.uint8)
    inverse=np.arange(n,dtype=np.int64)^key
    for b in range(blocks):
        block=enc[b*n:(b+1)*n]
        output[b*per:(b+1)*per]=block[inverse][1:]
    recovered=np.packbits(output[:original_bytes*8]).tobytes()
    if zlib.crc32(recovered)!=crc:
        raise ValueError('checksum failed')
    return recovered


def compression_test(rng):
    tests=corpus(rng)
    results=[]
    examples=[]
    for name,data in tests.items():
        z=zlib.compress(data,9)
        bro=brotli.compress(data,quality=9)
        zl=lzma.compress(data,preset=6)
        best=sorted([('zlib9',len(z)),('brotli9',len(bro)),('lzma6',len(zl))],key=lambda x:x[1])[0]
        choices=[]
        for n in (256,4096,65536,1048576):
            if n==1048576 and len(data)<65536:
                continue # giant block padding would dominate shorter files
            key=min(7,n-1)
            coded,enc_ms=repeated_time(lambda:encode_unit(data,n,key),2)
            recovered,dec_ms=repeated_time(lambda:decode_unit(coded),2)
            assert recovered==data and len(coded)==H_SIZE+((len(data)*8+(n-2))//(n-1)*n+7)//8
            cz=zlib.compress(coded,9)
            cb=brotli.compress(coded,quality=9)
            ch=min([('zlib9',len(cz)),('brotli9',len(cb))],key=lambda x:x[1])
            choices.append({'dimension':n,'uncompressed_signed_normalized_file_bytes':len(coded),
                            'zlib9_bytes':len(cz),'brotli9_bytes':len(cb),
                            'best_transformed_name':ch[0], 'best_transformed_bytes':ch[1],
                            'encode_ms_median':enc_ms,'decode_ms_median':dec_ms,
                            'roundtrip':True})
        selected=min(choices,key=lambda x:x['best_transformed_bytes'])
        results.append({'dataset':name,'original_bytes':len(data),
                        'zlib9_bytes':len(z),'brotli9_bytes':len(bro),'lzma6_bytes':len(zl),
                        'best_conventional':best[0],'best_conventional_bytes':best[1],
                        'best_cd_dimension':selected['dimension'],
                        'best_cd_bytes':selected['best_transformed_bytes'],
                        'pct_cd_minus_conventional':round(100*(selected['best_transformed_bytes']/best[1]-1),3),
                        'exact_reconstruction':True})
        examples.append({'dataset':name, 'baseline':best, 'variants':choices})
        print(f"[codec] {name:20s} bytes={len(data):7d} baseline={best[0]}:{best[1]:7d} imaginary={selected['dimension']:,}:{selected['best_transformed_bytes']:7d}", flush=True)
    return results,examples


def run():
    rng=np.random.default_rng(SEED)
    print('Cayley-Dickson / Clifford actual imaginary-algebra test; seed',SEED, flush=True)
    t0=time.perf_counter()
    # Warm up explicit compilation (excluded from reported timings)
    cd_sign(1,2,4)
    cd_unit_product(np.array([0,1,1,0],dtype=np.uint8),1)
    cd_product_dense(np.array([1,1,0,1],dtype=np.int8),np.array([0,1,1,0],dtype=np.int8))
    cd_product_sparse_indices(np.array([1],dtype=np.int64),np.array([1],dtype=np.int64),
                              np.array([2],dtype=np.int64),np.array([1],dtype=np.int64),4)
    algebra=tests_algebra(rng)
    print('[algebra] basis checks',algebra['basis_pairs_checked'],'zero divisor',algebra['sedenion_zero_divisor'], flush=True)
    dimensions=exact_dim_scaling(rng)
    print('[dimensions] exact to',dimensions[-1]['dimension_D'], 'median unit time', dimensions[-1]['unit_product_ms_median'], flush=True)
    dense,sparse=dense_product_scaling(rng)
    print('[products] dense to',dense[-1]['dimension'],'sparse to',sparse[-1]['dimension'], flush=True)
    clifford=clifford_test(rng)
    print('[clifford] through',clifford[-1]['generators_bits'],'imaginary generators',flush=True)
    compression, detailed = compression_test(rng)
    result={
        'date':'2026-11-18','source_paper':'Real-Pair Representations of Complex Numbers and Quaternions, November 4 2026',
        'scope':'Actual Cayley-Dickson imaginary coefficient encoding, multiplication and Clifford alternatives',
        'seed':SEED,'Python':platform.python_version(), 'numpy':np.__version__,'numba':__import__('numba').__version__,
        'algebra':algebra,'dimensions':dimensions,'dense_products':dense,'sparse_products':sparse,
        'clifford':clifford,'compression':compression,'compression_details':detailed,
        'runtime_sec':round(time.perf_counter()-t0,3)}
    (ROOT/'hypercomplex_highdim_results.json').write_text(json.dumps(result,indent=2)+'\n')
    for name,rows in [('hypercomplex_highdim_dimensions.csv',dimensions),
                      ('hypercomplex_highdim_dense_multiply.csv',dense),
                      ('hypercomplex_highdim_sparse_multiply.csv',sparse),
                      ('hypercomplex_highdim_compression.csv',compression),
                      ('hypercomplex_highdim_clifford.csv',clifford)]:
        with (ROOT/name).open('w',newline='') as f:
            dw=csv.DictWriter(f,fieldnames=list(rows[0]));dw.writeheader();dw.writerows(rows)
    print('[all] completed in',result['runtime_sec'],'s')
    return result

if __name__=='__main__':
    run()