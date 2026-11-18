#!/usr/bin/env python3
"""Reproduce 1,048,576-dimensional signed-imaginary codeword result."""
from pathlib import Path
import hashlib,json,time,zlib,lzma,brotli
import numpy as np
from hypercomplex_highdim_test import ROOT
from hypercomplex_cascade_test import cd_cascade_inverse_checked,encode_sparse_cascade,decode_sparse_cascade

def run():
    n=1048576; axis=123456
    spectrum=np.zeros(n,dtype=np.int64);spectrum[axis]=n
    beginning=time.perf_counter()
    bipolar=cd_cascade_inverse_checked(spectrum)
    if bipolar[0]<0:bipolar=-bipolar
    assert bipolar[0]==1 and np.all(np.abs(bipolar)==1)
    bitdata=np.packbits(((bipolar[1:]+1)//2).astype(np.uint8)).tobytes()
    coded,stats=encode_sparse_cascade(bitdata,n,True)
    decoded=decode_sparse_cascade(coded)
    assert decoded==bitdata and stats['spectral']==1 and stats['fill']==1
    sizes={
        'original':len(bitdata),'hypercomplex_raw':len(coded),
        'hypercomplex_zlib9':len(zlib.compress(coded,9)),
        'hypercomplex_brotli9':len(brotli.compress(coded,quality=9)),
        'zlib9':len(zlib.compress(bitdata,9)),
        'brotli9':len(brotli.compress(bitdata,quality=9)),
        'lzma6':len(lzma.compress(bitdata,preset=6)),
    }
    result={'dimension':n,'imaginary_axes':n-1,'axis':axis,
            'input_bits_including_byte_padding':len(bitdata)*8,
            'message_bits':n-1,'sha256':hashlib.sha256(bitdata).hexdigest(),
            'sizes_bytes':sizes,'codec_modes':stats,'exact_roundtrip':True,
            'wall_sec':round(time.perf_counter()-beginning,3)}
    (ROOT/'hypercomplex_million_dim_demo.hcsp').write_bytes(coded)
    (ROOT/'hypercomplex_million_dim_codeword_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result,flush=True)
    return result

if __name__=='__main__':run()