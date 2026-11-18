#!/usr/bin/env python3
"""Multi-bit imaginary coefficients and channel-error tests for high-dimensional encoding."""
from __future__ import annotations
import csv, json, sys, time, statistics
from pathlib import Path
import numpy as np
from numba import njit
from hypercomplex_highdim_test import ROOT, SEED, cd_sign, repeated_time
from hypercomplex_cascade_test import cd_cascade, cd_cascade_inverse_checked

@njit(cache=True)
def unit_wide(a,key):
    n=len(a)
    out=np.zeros(n,dtype=np.int64)
    for i in range(n):
        if a[i]!=0:
            out[i^key]=int(a[i])*cd_sign(key,i,n)
    return out

@njit(cache=True)
def unit_wide_inverse(b,key):
    n=len(b)
    out=np.zeros(n,dtype=np.int64)
    for i in range(n):
        out[i]=int(b[i^key])*cd_sign(key,i,n)
    return out

def test():
    rng=np.random.default_rng(SEED+900)
    unit_wide(np.array([0,1],dtype=np.int32),1)
    unit_wide_inverse(np.array([0,1],dtype=np.int64),1)
    rows=[]
    for n in (256,65536,1048576):
        for w in (1,2,4,8,16):
            c=np.zeros(n,dtype=np.int32)
            c[1:]=rng.integers(0,1<<w,size=n-1,dtype=np.int32)
            signed,enc_ms=repeated_time(lambda:unit_wide(c,7 if n>7 else 1),3)
            recovered,dec_ms=repeated_time(lambda:unit_wide_inverse(signed,7 if n>7 else 1),3)
            assert np.array_equal(c,recovered)
            assert np.count_nonzero(c)==np.count_nonzero(signed)
            # Integer pipeline; no floating point means exact beyond 2**53 if needed.
            rows.append({'dimensions':n,'imaginary_axes':n-1,'bits_per_coefficient':w,
                         'input_bits':(n-1)*w,'tight_storage_bytes':((n-1)*w+7)//8,
                         'dense_int32_bytes':int(c.nbytes),
                         'signed_int64_product_bytes':int(signed.nbytes),
                         'nonzero_coefficients':int(np.count_nonzero(c)),
                         'multiplier_basis_key':7,'roundtrip_exact':True,
                         'unit_encode_ms':enc_ms,'unit_decode_ms':dec_ms})
        print('[widths] n',n,'done',flush=True)
    signed_mix=[]
    for n,w in ((256,1),(65536,8),(65536,16),(1048576,16)):
        c=np.zeros(n,dtype=np.int32)
        c[1:]=rng.integers(0,1<<w,size=n-1,dtype=np.int32)
        mixed,enc_ms=repeated_time(lambda:cd_cascade(c),2)
        recovered,dec_ms=repeated_time(lambda:cd_cascade_inverse_checked(mixed),2)
        assert np.array_equal(recovered,c)
        norm_in=int(np.dot(c.astype(np.int64),c.astype(np.int64)))
        norm_out=sum(int(z)*int(z) for z in mixed)
        assert norm_out==n*norm_in
        signed_mix.append({'dimensions':n,'width_bits':w,
                           'max_input_coeff':int(c.max()),
                           'max_output_abs_coeff':int(np.max(np.abs(mixed))),
                           'nonzero_input_coeff':int(np.count_nonzero(c)),
                           'nonzero_mixed_coeff':int(np.count_nonzero(mixed)),
                           'norm_multiplier':n,'encode_ms':enc_ms,'decode_ms':dec_ms,
                           'roundtrip_exact':True})
        print('[mixed width] n',n,'bits',w,'done',flush=True)
    # One physical bit flip in a sign-normalized unit-encoded stream -> exactly one data-bit flip.
    n=4096;key=7;trials=1000
    msgs=rng.integers(0,2,size=(trials,n),dtype=np.uint8)
    msgs[:,0]=0
    hamming=[]
    for i in range(trials):
        coded=np.abs(unit_wide(msgs[i].astype(np.int32),key)).astype(np.uint8)
        flip=int(rng.integers(0,n))
        # Exclude real-axis index of original message; choose valid image coordinate.
        while (flip^key)==0:flip=int(rng.integers(0,n))
        coded[flip]^=1
        corrupt=coded[np.arange(n)^key]
        hamming.append(int(np.sum(corrupt!=msgs[i])))
    assert set(hamming)=={1}
    report={'seed':SEED+900,'multi_bit_widths':rows,'mixed_widths':signed_mix,
            'unit_channel_bit_flip':{'trials':trials,'input_dimension':n,
                                     'single_encoded_bit_flips':trials,
                                     'data_bits_flipped_per_trial_unique':list(set(hamming)),
                                     'native_error_correction':False}}
    (ROOT/'hypercomplex_widths_results.json').write_text(json.dumps(report,indent=2)+'\n')
    with (ROOT/'hypercomplex_widths.csv').open('w',newline='') as f:
        dw=csv.DictWriter(f,fieldnames=rows[0]);dw.writeheader();dw.writerows(rows)
    print('[widths] success',len(rows),'configurations',flush=True)
    return report

if __name__=='__main__':test()