"""Experimental lossless *causal* hypercomplex-multiplication predictor.

A predictor only sees two PREVIOUS byte vectors at each time t. The trained
model is fitted to an initial prefix (up to 25% of the blocks, capped at 256) of the SAME file; the prefix is sent in the
compressed frame, and frozen, quantized model coefficients are included explicitly.

Cayley-Dickson features are genuine bilinear products using the source
paper's CD algebra, not a linear coordinate transform or a future-data oracle.
Residuals are serialized with an invertible modulo-256 zigzag map.
"""
from __future__ import annotations
import struct
import zlib
import brotli
import numpy as np
from numba import njit
from hypercomplex_highdim_test import cd_sign
from sota_benchmark.codec import Zstd

MAGIC = b'HCP2'
VERSION = 1
Q = 1024
HEADER = struct.Struct('<4sBBBBQIII')  # magic, version, log2 D, model, compressor, length, prefix frames, crc, coefficients
LENGTHS = struct.Struct('<II')
MODELS = ('ar2', 'ar2_cd', 'ar2_xorconv', 'ar2_pointwise', 'ar2_circular', 'ar2_shifted', 'lag1', 'delta2', 'xor_lag')
COMPRESSORS = ('brotli11', 'zstd9', 'zstd3')
MODEL_IDS = {x: i for i, x in enumerate(MODELS)}
COMP_IDS = {x: i for i, x in enumerate(COMPRESSORS)}
TRAINED = frozenset(('ar2', 'ar2_cd', 'ar2_xorconv', 'ar2_pointwise', 'ar2_circular', 'ar2_shifted'))
EXTRA = {'ar2_cd':1,'ar2_xorconv':2,'ar2_pointwise':3,'ar2_circular':4,'ar2_shifted':5}


def _zstd():
    return Zstd()


def pack_stream(data: bytes, backend: str) -> bytes:
    if backend=='brotli11':return brotli.compress(data, quality=11)
    if backend=='zstd9':return _zstd().compress(data,9)
    if backend=='zstd3':return _zstd().compress(data,3)
    raise ValueError(backend)


def unpack_stream(data: bytes, backend: str, length: int) -> bytes:
    if backend=='brotli11':
        out=brotli.decompress(data)
        if len(out)!=length:raise ValueError('incorrect Brotli output length')
        return out
    if backend=='zstd9' or backend=='zstd3':return _zstd().decompress(data,length)
    raise ValueError(backend)


@njit(cache=True)
def pair_feature(x, y, which):
    """Int64 channelwise interaction from two *previous* centered byte blocks.

    which: 1 actual Cayley-Dickson product, 2 unsigned XOR convolution,
           3 coordinatewise multiplication, 4 ordinary cyclic convolution,
           5 single shifted signed product. Each normalized to byte scale.
    """
    D=len(x)
    result=np.zeros(D,dtype=np.int64)
    if which==3:
        for k in range(D):result[k]=x[k]*y[k]
        den=128
    elif which==5:
        for k in range(D):result[k]=x[k]*y[k^1]*cd_sign(k, k^1, D)
        den=128
    else:
        for i in range(D):
            a=x[i]
            for j in range(D):
                if which==1:
                    result[i^j] += a*y[j]*cd_sign(i,j,D)
                elif which==2:
                    result[i^j] += a*y[j]
                else:
                    result[(i+j)%D] += a*y[j]
        den=128*int(np.sqrt(D))
    for k in range(D):
        z=result[k]
        if z>=0: result[k]=(z+den//2)//den
        else: result[k]=-((-z+den//2)//den)
    return result


def feature_matrix(prev1: np.ndarray, prev2: np.ndarray, mode: str) -> np.ndarray:
    D=len(prev1)
    p1=prev1.astype(np.int64)-128
    p2=prev2.astype(np.int64)-128
    intercept=np.ones(D,dtype=np.int64)
    if mode=='ar2': return np.column_stack((intercept,p1,p2))
    if mode in EXTRA:
        bilinear=pair_feature(p1,p2,EXTRA[mode])
        return np.column_stack((intercept,p1,p2,bilinear))
    raise ValueError('unsupported fitted model '+mode)


def fit_coefficients(data: bytes, D: int, prefix_frames: int, mode: str) -> np.ndarray:
    if mode not in TRAINED:return np.zeros(0,dtype=np.int32)
    ncol=4 if mode in EXTRA else 3
    x=np.frombuffer(data[:prefix_frames*D], dtype=np.uint8).reshape((-1,D))
    F=[]; targets=[]
    for t in range(2,len(x)):
        F.append(feature_matrix(x[t-1],x[t-2],mode))
        targets.append(x[t].astype(np.int64)-128)
    if not F: return np.zeros(ncol,dtype=np.int32)
    mat=np.concatenate(F).astype(np.float64)
    labels=np.concatenate(targets).astype(np.float64)
    # Identical ridge penalty for matched-capacity feature models.
    gram=mat.T@mat
    gram+=np.eye(ncol)*0.1
    try:weights=np.linalg.solve(gram,mat.T@labels)
    except np.linalg.LinAlgError:weights=np.zeros(ncol)
    weights=np.clip(weights,-127.,127.)
    return np.rint(weights*Q).astype(np.int32)


def div_round(x:np.ndarray,denom:int=Q)->np.ndarray:
    return np.where(x>=0,(x+denom//2)//denom,-((-x+denom//2)//denom))


def predict(prev1:np.ndarray,prev2:np.ndarray,mode:str,coef:np.ndarray)->np.ndarray:
    if mode=='lag1':return prev1.copy()
    if mode=='delta2':return np.clip(2*prev1.astype(np.int64)-prev2.astype(np.int64),0,255).astype(np.uint8)
    if mode=='xor_lag':return np.bitwise_xor(prev1,prev2)
    if mode not in TRAINED:raise ValueError(mode)
    F=feature_matrix(prev1,prev2,mode)
    out=128+div_round(F@coef.astype(np.int64))
    return np.clip(out,0,255).astype(np.uint8)


def encode_residual(actual:np.ndarray, predicted:np.ndarray) -> np.ndarray:
    d=(actual.astype(np.int64)-predicted.astype(np.int64))%256
    signed=np.where(d<128,d,d-256)
    return np.where(signed>=0,2*signed,(-signed)*2-1).astype(np.uint8)


def decode_residual(encoded:np.ndarray, predicted:np.ndarray)->np.ndarray:
    e=encoded.astype(np.int64)
    signed=np.where(e%2==0,e//2,-((e+1)//2))
    return ((predicted.astype(np.int64)+signed)%256).astype(np.uint8)


@njit(cache=True)
def _predict_jit(p1,p2,model_id,coeff):
    D=len(p1)
    out=np.empty(D,dtype=np.uint8)
    if model_id==6:
        return p1.copy()
    if model_id==7:
        for i in range(D):
            v=2*int(p1[i])-int(p2[i])
            out[i]=np.uint8(min(255,max(0,v)))
        return out
    if model_id==8:
        for i in range(D):out[i]=p1[i]^p2[i]
        return out
    a=np.empty(D,dtype=np.int64);b=np.empty(D,dtype=np.int64)
    for i in range(D):
        a[i]=int(p1[i])-128
        b[i]=int(p2[i])-128
    extra=np.zeros(D,dtype=np.int64)
    if model_id in (1,2,3,4,5):extra=pair_feature(a,b,model_id)
    for i in range(D):
        v=int(coeff[0])+int(coeff[1])*a[i]+int(coeff[2])*b[i]
        if model_id in (1,2,3,4,5):v+=int(coeff[3])*extra[i]
        if v>=0:v=(v+Q//2)//Q
        else:v=-((-v+Q//2)//Q)
        out[i]=np.uint8(min(255,max(0,128+v)))
    return out


@njit(cache=True)
def _encode_residual_jit(frames,train,model_id,coeff):
    n,D=frames.shape
    out=np.empty((n-train)*D,dtype=np.uint8)
    pos=0
    for t in range(train,n):
        pred=_predict_jit(frames[t-1],frames[t-2],model_id,coeff)
        for i in range(D):
            delta=(int(frames[t,i])-int(pred[i]))%256
            signed=delta if delta<128 else delta-256
            v=2*signed if signed>=0 else -2*signed-1
            out[pos]=np.uint8(v);pos+=1
    return out


@njit(cache=True)
def _decode_residual_jit(prefix,residual,nframes,D,train,model_id,coeff):
    output=np.empty(nframes*D,dtype=np.uint8)
    for j in range(len(prefix)):output[j]=prefix[j]
    for t in range(train,nframes):
        p1=output[(t-1)*D:t*D]
        p2=output[(t-2)*D:(t-1)*D]
        pred=_predict_jit(p1,p2,model_id,coeff)
        for i in range(D):
            e=int(residual[(t-train)*D+i])
            signed=e//2 if e%2==0 else -((e+1)//2)
            output[t*D+i]=np.uint8((int(pred[i])+signed)%256)
    return output


def residual_stream(data:bytes,mode:str,D:int,prefix_frames:int,coef:np.ndarray) -> tuple[bytes,bytes]:
    nframes=len(data)//D
    x=np.frombuffer(data[:nframes*D],dtype=np.uint8).reshape((nframes,D))
    if nframes>prefix_frames and prefix_frames<2:
        raise ValueError('two previous blocks required')
    encoded=_encode_residual_jit(x,prefix_frames,MODEL_IDS[mode],coef)
    return data[:prefix_frames*D],encoded.tobytes()+data[nframes*D:]


def rebuild(prefix:bytes,resid:bytes,mode:str,D:int,prefix_frames:int,coef:np.ndarray,total_length:int)->bytes:
    if len(prefix)!=prefix_frames*D:raise ValueError('invalid training prefix length')
    nframes=total_length//D
    if len(resid)!=total_length-prefix_frames*D:raise ValueError('residual length mismatch')
    if nframes>prefix_frames and prefix_frames<2:raise ValueError('missing two previous blocks')
    p=np.frombuffer(prefix,dtype=np.uint8)
    r=np.frombuffer(resid[:(nframes-prefix_frames)*D],dtype=np.uint8)
    raw=_decode_residual_jit(p,r,nframes,D,prefix_frames,MODEL_IDS[mode],coef)
    return raw.tobytes()+resid[(nframes-prefix_frames)*D:]


def default_prefix_frames(length:int,D:int)->int:
    n=length//D
    if n<3:return n
    return min(n,max(3,min(256,n//4)))


def encode(data:bytes,mode='ar2_cd',D=16,backend='zstd9')->tuple[bytes,dict]:
    if D not in (2,4,8,16,32,64):raise ValueError('unsupported D')
    if mode not in MODEL_IDS or backend not in COMP_IDS:raise ValueError('bad predictor or compressor')
    train=default_prefix_frames(len(data),D)
    if train<2 and len(data)>=D*2:raise ValueError('insufficient prefix')
    coef=fit_coefficients(data,D,train,mode)
    prefix, residual=residual_stream(data,mode,D,train,coef)
    cp=pack_stream(prefix,backend)
    cr=pack_stream(residual,backend)
    frame=(HEADER.pack(MAGIC,VERSION,D.bit_length()-1,MODEL_IDS[mode],COMP_IDS[backend],len(data),train,zlib.crc32(data),len(coef))+
           coef.astype('<i4').tobytes()+LENGTHS.pack(len(cp),len(cr))+cp+cr)
    return frame, {'prefix_frames':train,'prefix_bytes':len(prefix),'residual_bytes':len(residual),
                   'compressed_prefix_bytes':len(cp),'compressed_residual_bytes':len(cr),
                   'coefficient_bytes':4*len(coef),'model':mode,'D':D,'backend':backend}


def decode(frame:bytes)->bytes:
    if len(frame)<HEADER.size+LENGTHS.size:raise ValueError('truncated prediction frame')
    magic,v,shift,mode_id,backend_id,total,train,crc,ncoef=HEADER.unpack_from(frame)
    if magic!=MAGIC or v!=VERSION or shift>6 or shift<1:raise ValueError('wrong predictive frame')
    if mode_id>=len(MODELS) or backend_id>=len(COMPRESSORS):raise ValueError('unknown model or compressor')
    if total>1<<30 or ncoef>4:raise ValueError('unsafe length / coefficient count')
    D=1<<shift
    mode=MODELS[mode_id]; backend=COMPRESSORS[backend_id]
    if ncoef != (4 if mode in EXTRA else 3 if mode=='ar2' else 0):raise ValueError('invalid coefficient count')
    if train!=default_prefix_frames(total,D):raise ValueError('invalid training prefix')
    pos=HEADER.size
    if pos+4*ncoef+LENGTHS.size>len(frame):raise ValueError('missing coefficients')
    coef=np.frombuffer(frame[pos:pos+4*ncoef],dtype='<i4').copy()
    pos+=4*ncoef
    a,b=LENGTHS.unpack_from(frame,pos);pos+=LENGTHS.size
    if pos+a+b!=len(frame):raise ValueError('invalid frame lengths')
    prefix=unpack_stream(frame[pos:pos+a],backend,train*D)
    resid=unpack_stream(frame[pos+a:pos+a+b],backend,total-train*D)
    raw=rebuild(prefix,resid,mode,D,train,coef,total)
    if len(raw)!=total or zlib.crc32(raw)!=crc:raise ValueError('CRC mismatch')
    return raw


def select_prefix_model(data:bytes,D:int=16,candidates:tuple[str,...]=MODELS)->tuple[str,dict]:
    """Choose mode ONLY from an initial within-file prefix, never suffix."""
    prefix_frames=default_prefix_frames(len(data),D)
    if prefix_frames<12:return 'lag1',{'decision':'short_prefix_fallback'}
    split=max(4,(prefix_frames*3)//4)
    blocks=np.frombuffer(data[:prefix_frames*D],dtype=np.uint8).reshape((-1,D))
    scores={}
    for mode in candidates:
        coef=fit_coefficients(data,D,split,mode)
        buf=bytearray()
        for t in range(split,prefix_frames):
            p=predict(blocks[t-1],blocks[t-2],mode,coef)
            buf.extend(encode_residual(blocks[t],p).tobytes())
        # Same compressed validation bytes, plus charge serialized coefficients.
        scores[mode]=len(pack_stream(bytes(buf),'zstd3'))+4*len(coef)
    return min(candidates,key=lambda m:(scores[m],m!='ar2',m)),scores


def encode_adaptive(data:bytes,D:int=16,backend:str='zstd9')->tuple[bytes,dict]:
    mode,scores=select_prefix_model(data,D)
    framed,info=encode(data,mode,D,backend)
    info['selection']='prefix_only'
    info['validation_scores']=scores
    return framed,info
