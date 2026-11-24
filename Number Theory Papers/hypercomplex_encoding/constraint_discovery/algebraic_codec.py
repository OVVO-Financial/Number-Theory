#!/usr/bin/env python3
"""Exactly reversible fixed-point quaternion constraint probes.

Only integer-valued, quantized input streams are encoded losslessly. The source
CSV numeric text itself is NOT reconstructed: see benchmark_real.py.
"""
from __future__ import annotations
import math
import struct
import zlib
from dataclasses import dataclass

import numpy as np

SCALE = 1_000_000
MODES = ('raw','delta','ar2','mul','norm','norm_ar2','norm_delta','pointwise')
HEADER = '<4sHIII'
MAGIC = b'ACQ1'
HEADER_SIZE = struct.calcsize(HEADER)


def mul(q, r):
    """Hamilton product, scalar-first (w,x,y,z)."""
    a,b,c,d = [float(v) for v in q]
    e,f,g,h = [float(v) for v in r]
    return np.array([a*e-b*f-c*g-d*h,
                     a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f,
                     a*h+b*g-c*f+d*e], dtype=np.float64)


def conjugate(q):
    a,b,c,d = q
    return np.array([a,-b,-c,-d], dtype=np.float64)


def as_fixed(q, scale=SCALE):
    q = np.asarray(q, dtype=np.float64)
    if q.ndim != 2 or q.shape[1] != 4 or not np.all(np.isfinite(q)):
        raise ValueError('require N x 4 finite array')
    if np.max(np.abs(q), initial=0) * scale >= (2**31 - 1):
        raise ValueError('quantized values overflow i32')
    return np.rint(q * scale).astype('<i4')


def predict(history, mode, channel, current=None, scale=SCALE):
    """Predict a single fixed-point coefficient from reconstructed predecessors.

    In norm modes xyz are encoded first in the current row, allowing a
    deterministic norm prediction for w. w is channel 0 in scalar-first data.
    """
    idx=len(history)
    if mode=='raw' or idx==0:
        return 0
    prev = history[-1]
    if mode in ('norm', 'norm_ar2','norm_delta') and channel == 0:
        if current is None: raise ValueError('norm prediction needs current xyz')
        xyz = np.asarray(current[1:4],dtype=np.float64) / scale
        sq = 1.0 - float(np.dot(xyz,xyz))
        w_hat = math.sqrt(max(0.0, sq))
        sign = -1.0 if prev[0] < 0 else 1.0
        return int(round(sign*w_hat*scale))
    if mode in ('delta','norm','norm_delta') or idx<2:
        return int(prev[channel])
    earlier = history[-2]
    if mode in ('ar2','norm_ar2'):
        return int(2*int(prev[channel])-int(earlier[channel]))
    if mode == 'pointwise':
        p=float(prev[channel]);pp=float(earlier[channel]);
        if abs(pp)<500: return int(2*p-pp)
        z = p * p / pp
        return int(round(np.clip(z, -2*scale,2*scale)))
    if mode=='mul':
        q1=np.asarray(prev,dtype=np.float64)/scale
        q0=np.asarray(earlier,dtype=np.float64)/scale
        denominator=float(np.dot(q0,q0))
        if denominator < 1e-9: return int(2*int(prev[channel])-int(earlier[channel]))
        speed=mul(q1,conjugate(q0))/denominator
        p=mul(speed,q1)
        return int(round(np.clip(p[channel]*scale, -2*scale,2*scale)))
    raise ValueError(mode)


def _predict_row(hist, mode, now, scale=SCALE):
    out=np.zeros(4,dtype=np.int64)
    # Channels 1-3 must precede channel 0 for norm models.
    order=(1,2,3,0) if mode.startswith('norm') else (0,1,2,3)
    for c in order:
        out[c]=predict(hist,mode,c,now,scale)
    return out


def residuals(arr, mode, scale=SCALE):
    arr=np.asarray(arr,dtype=np.int32)
    if arr.ndim!=2 or arr.shape[1]!=4: raise ValueError('must be N x 4')
    if mode not in MODES:raise ValueError(mode)
    residual=np.empty(arr.shape,dtype=np.int32)
    # Compression is a causal transform; encode uses original already seen rows,
    # decoder uses reconstructed rows. An explicit decoder verifies equality.
    hist=[]
    for t in range(len(arr)):
        p=_predict_row(hist,mode,arr[t],scale)
        error=arr[t].astype(np.int64)-p
        if np.any(error>2**31-1) or np.any(error<-(2**31)):
            raise OverflowError('residual overflow')
        residual[t]=error.astype(np.int32)
        hist.append(arr[t])
    return residual


def recover(res, mode, scale=SCALE):
    res=np.asarray(res,dtype=np.int32)
    if res.ndim!=2 or res.shape[1]!=4:raise ValueError('must be N x 4')
    hist=[]
    for row in res:
        # Recover xyz first in norm mode. Other predictors do not depend on
        # current row; all channels can be reconstructed independently.
        if mode.startswith('norm'):
            cur=np.empty(4,dtype=np.int64)
            for ch in (1,2,3):
                cur[ch]=int(row[ch])+predict(hist,mode,ch,cur,scale)
            cur[0]=int(row[0])+predict(hist,mode,0,cur,scale)
        else:
            cur=np.empty(4,dtype=np.int64)
            for ch in (0,1,2,3):
                cur[ch]=int(row[ch])+predict(hist,mode,ch,cur,scale)
        if np.any(cur>2**31-1) or np.any(cur<-(2**31)):
            raise OverflowError('decoded coefficient overflow')
        hist.append(cur.astype(np.int32))
    return np.asarray(hist,dtype='<i4').reshape((-1,4))


def pack_uvarint(n):
    if n<0:raise ValueError(n)
    output=bytearray()
    while n>127:
        output.append((n&127)|128)
        n>>=7
    output.append(n)
    return bytes(output)


def read_uvarint(data, offset):
    out=0
    for j in range(10):
        if offset>=len(data):raise ValueError('truncated varint')
        b=data[offset];offset+=1
        out|=(b&127)<<(7*j)
        if b<128:return out,offset
    raise ValueError('overlong varint')


def zigzag(v):return (v<<1) ^ (v>>31)
def unzigzag(v):return (v>>1) ^ -(v&1)


def pack_payload(res, encoding):
    if encoding=='i32_row':return res.astype('<i4').tobytes(order='C')
    if encoding=='i32_col':return res.T.copy().astype('<i4').tobytes(order='C')
    if encoding=='varint_col':
        bs=bytearray()
        for ch in res.T:
            for v in ch:
                bs+=pack_uvarint(zigzag(int(v)))
        return bytes(bs)
    raise ValueError(encoding)


def unpack_payload(data, rows, encoding):
    if encoding=='i32_row':
        if len(data)!=rows*16:raise ValueError('invalid i32 row length')
        return np.frombuffer(data,dtype='<i4').reshape(rows,4).copy()
    if encoding=='i32_col':
        if len(data)!=rows*16:raise ValueError('invalid i32 column length')
        return np.frombuffer(data,dtype='<i4').reshape(4,rows).T.copy()
    if encoding=='varint_col':
        out=np.zeros((4,rows),dtype=np.int32)
        pos=0
        for ch in range(4):
            for i in range(rows):
                val,pos=read_uvarint(data,pos)
                a=unzigzag(val)
                if val>2**32-1 or a>=2**31 or a<-(2**31):raise ValueError('invalid varint coefficient')
                out[ch,i]=a
        if pos!=len(data):raise ValueError('trailing varint data')
        return out.T.copy()
    raise ValueError(encoding)


def encode(arr, mode, encoding, compressor='zlib'):
    arr=np.asarray(arr,dtype='<i4')
    rr=residuals(arr,mode)
    data=pack_payload(rr,encoding)
    if compressor=='zlib':payload=zlib.compress(data,9)
    elif compressor=='brotli':
        import brotli
        payload=brotli.compress(data,quality=11)
    elif compressor=='zstd':
        import zstandard as zstd
        payload=zstd.ZstdCompressor(level=9).compress(data)
    elif compressor=='none':payload=data
    else:raise ValueError(compressor)
    mid=MODES.index(mode)
    eid=('i32_row','i32_col','varint_col').index(encoding)
    cid=('none','zlib','brotli','zstd').index(compressor)
    meta=(mid<<16)|(eid<<8)|cid
    return struct.pack(HEADER,MAGIC,1,len(arr),meta,zlib.crc32(arr.tobytes()))+payload


def decode(blob):
    if len(blob)<HEADER_SIZE:raise ValueError('missing header')
    magic,version,rows,meta,checksum=struct.unpack_from(HEADER,blob,0)
    if magic!=MAGIC or version!=1 or rows>10_000_000:raise ValueError('invalid metadata')
    mid=(meta>>16)&255;eid=(meta>>8)&255;cid=meta&255
    if mid>=len(MODES) or eid>2 or cid>3:raise ValueError('invalid codec settings')
    data=blob[HEADER_SIZE:]
    if cid==1:data=zlib.decompress(data)
    elif cid==2:
        import brotli
        data=brotli.decompress(data)
    elif cid==3:
        import zstandard as zstd
        data=zstd.ZstdDecompressor().decompress(data)
    res=unpack_payload(data,rows,('i32_row','i32_col','varint_col')[eid])
    recovered=recover(res,MODES[mid])
    if zlib.crc32(recovered.tobytes())!=checksum:raise ValueError('checksum mismatch')
    return recovered


def characteristics(a,scale=SCALE):
    arr=np.asarray(a,dtype='<i4')
    q=arr.astype(np.float64)/scale
    norms=np.sum(q*q,axis=1)
    xyz_norm=np.linalg.norm(q[:,1:],axis=1)
    err_norm=np.abs(norms-1)
    metrics={
      'norm_abs_max':float(np.max(err_norm)),
      'norm_abs_median':float(np.median(err_norm)),
      'norm_within_1e-3_frac':float(np.mean(err_norm<1e-3)),
      'norm_within_1e-4_frac':float(np.mean(err_norm<1e-4)),
      'normalized_w_residual_zero_pct':float(np.mean(residuals(arr,'norm')[:,0]==0)),
      'normalized_w_residual_median_abs':float(np.median(np.abs(residuals(arr,'norm')[:,0].astype(np.int64)))),
    }
    return metrics
