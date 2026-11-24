#!/usr/bin/env python3
"""Exact fixed-point, multi-constraint joint encoding of EuRoC 17-column state trajectories.

Input is first quantized to signed integer micro-units; this is NOT a byte-identical
representation of the source CSV. Position p, orientation q, velocity v, gyro and
accel biases b are preserved identically at that integer precision.

Independent hypotheses:
 N: quaternion scalar coordinate from unit norm and current decoded vector part.
 K: translation kinematics p_t from velocity v_t and the prior state, OR inverse.
 R: near-constant quaternion relative-rotation product, for q vector components.
 B: near-constant gyro/accelerometer biases, before residual correction.

N and R share orientation but are distinct static and temporal restrictions; K concerns
translation/velocity, B concerns sensor bias. Every prediction has an exact residual.
"""
from __future__ import annotations
import csv,hashlib,json,math,struct,zlib
from pathlib import Path
import numpy as np
import brotli

SCALE=1_000_000
GROUPS={"position":(0,1,2),"quaternion":(3,4,5,6),"velocity":(7,8,9),"bias":(10,11,12,13,14,15)}
MODE_INDEX={"raw":0,"delta":1,"ar2":2}
INDEX_MODE={v:k for k,v in MODE_INDEX.items()}
HEADER='<4sBIIIII'
HEADER_SIZE=struct.calcsize(HEADER)
MAGIC=b'JAC1'


def parse_euroc_csv(path,limit=None):
    """EuRoC groundtruth (time ns,p3,q4,v3,bias_g3,bias_a3), preserving only 1e-6 integer grid."""
    rows=[]; t=[]; raw=Path(path).read_bytes(); skipped=0
    for line in raw.decode('utf-8').splitlines():
        if not line.strip() or line.lstrip().startswith('#'): continue
        p=line.strip().split(',')
        if len(p)<17: skipped+=1;continue
        try:
            # direct timestamps as int64 nanoseconds, no floating timestamp precision loss
            ti=int(float(p[0])) if ('e' in p[0].lower() or '.' in p[0]) else int(p[0])
            x=np.asarray([float(c) for c in p[1:17]],dtype=float)
            if not np.all(np.isfinite(x)):skipped+=1;continue
            vals=np.rint(x*SCALE)
            if np.any(np.abs(vals)>2**31-1):skipped+=1;continue
        except (ValueError,OverflowError):skipped+=1;continue
        t.append(ti);rows.append(vals.astype(np.int32))
        if limit and len(rows)>=limit:break
    if len(rows)<100:raise ValueError('Fewer than 100 valid 17-channel EuRoC samples from '+str(path)+', skipped '+str(skipped))
    t=np.asarray(t,dtype=np.int64)
    if not np.all(np.diff(t)>0):
        mask=np.r_[True,np.diff(t)>0]
        t=t[mask]; rows=np.asarray(rows,dtype=np.int32)[mask]
    else:rows=np.asarray(rows,dtype=np.int32)
    return t,rows,{"source_sha256":hashlib.sha256(raw).hexdigest(),"rows":len(t),"skipped":skipped,
                   "timestamp_median_ns":int(np.median(np.diff(t)))}


def parse_tum_pose(path,limit=None):
    """TUM estimate (time seconds,p3,qx,qy,qz,qw) to 7 fixed columns; no v/bias channels."""
    raw=Path(path).read_bytes(); t=[]; rows=[]
    for line in raw.decode().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):continue
        p=line.split()
        if len(p)<8:continue
        try:
            v=list(map(float,p[:8])); assert all(math.isfinite(x) for x in v)
            ti=int(round(v[0]*1e9)); x=np.array([*v[1:4],v[7],*v[4:7]])
            y=np.rint(x*SCALE).astype(np.int32)
        except (ValueError,OverflowError,AssertionError):continue
        if t and ti<=t[-1]:continue
        t.append(ti);rows.append(y)
        if limit and len(t)>=limit:break
    if len(t)<100:raise ValueError('Insufficient TUM pose observations')
    return np.array(t,dtype=np.int64),np.array(rows,dtype=np.int32),{"source_sha256":hashlib.sha256(raw).hexdigest(),"rows":len(t)}


def _qmul(q,r):
    a,b,c,d=q; e,f,g,h=r
    return np.array([a*e-b*f-c*g-d*h,a*f+b*e+c*h-d*g,a*g-b*h+c*e+d*f,a*h+b*g-c*f+d*e])


def qrot_pred(prev,last):
    # q_{t-1} * conjugate(q_{t-2}) * q_{t-1} / ||q_{t-2}||^2
    u=prev.astype(float)/SCALE;v=last.astype(float)/SCALE
    den=float(np.dot(v,v))
    if den<1e-9:return np.rint(2*prev.astype(float)-last).astype(np.int64)
    qq=_qmul(_qmul(u,np.array([v[0],-v[1],-v[2],-v[3]])),u)/den
    return np.rint(np.clip(qq,-2.0,2.0)*SCALE).astype(np.int64)


def norm_pred(xyz,prev_w):
    vv=sum(int(v)*int(v) for v in xyz)
    rem=max(0,SCALE*SCALE-vv)
    a=math.isqrt(rem)
    if (a+1)*(a+1)-rem<rem-a*a:a+=1
    return -a if prev_w<0 else a


def round_div_signed(num,den):
    """Symmetric half-away-from-zero, deterministic for integer arithmetic."""
    if den<=0:raise ValueError('denominator')
    return (num+den//2)//den if num>=0 else -((-num+den//2)//den)


def temporal_pred(hist,ch,mode):
    if not hist or mode=='raw':return 0
    if mode=='ar2' and len(hist)>=2:return 2*int(hist[-1][ch])-int(hist[-2][ch])
    return int(hist[-1][ch])


def all_group_config(mode='ar2'):
    return {g:mode for g in GROUPS}


def _index_cfg(cfg):
    return sum(MODE_INDEX[cfg[g]]<<(2*i) for i,g in enumerate(GROUPS))


def _decode_cfg(v):
    return {g:INDEX_MODE[(v>>(2*i))&3] for i,g in enumerate(GROUPS)}


def build_residuals(t,source,base=None,norm=False,rotation=False,kinematic='none',bias_delta=False):
    """True sequential encoder, each current-row dependency decoded before use."""
    if base is None:base=all_group_config()
    if len(t)!=len(source):raise ValueError('timestamps')
    data=np.asarray(source,dtype=np.int64)
    channels=data.shape[1]
    if channels not in (7,16):raise ValueError('7- or 16-column source required')
    if kinematic!='none' and channels!=16:raise ValueError('kinematic needs measured v')
    errors=np.zeros_like(data,dtype=np.int64)
    past=[]
    for i,cur in enumerate(data):
        prev=past[-1] if past else None
        # Source-specific rotation hypothesis predicts current quaternion from two prior orientations.
        rp=qrot_pred(past[-1][3:7],past[-2][3:7]) if rotation and len(past)>1 else None
        order=[7,8,9,4,5,6,3,0,1,2,10,11,12,13,14,15] if channels==16 else [4,5,6,3,0,1,2]
        if kinematic=='v_from_p':order=[0,1,2,4,5,6,3,7,8,9,10,11,12,13,14,15]
        for ch in order:
            if ch<3: grp='position'
            elif ch<7:grp='quaternion'
            elif ch<10:grp='velocity'
            else:grp='bias'
            pred=temporal_pred(past,ch,'delta' if bias_delta and grp=='bias' else base[grp])
            if rp is not None and 3<=ch<7:pred=int(rp[ch-3])
            if norm and ch==3 and prev is not None:
                pred=norm_pred(cur[4:7],int(prev[3]))
            if prev is not None and kinematic=='p_from_v' and ch<3:
                delta=int(t[i]-t[i-1]);pred=int(prev[ch])+round_div_signed(delta*(int(prev[ch+7])+int(cur[ch+7])),2_000_000_000)
            if prev is not None and kinematic=='v_from_p' and 7<=ch<10:
                delta=int(t[i]-t[i-1]);pred=round_div_signed((int(cur[ch-7])-int(prev[ch-7]))*1_000_000_000,delta)
            errors[i,ch]=int(cur[ch])-pred
        past.append(cur)
    return errors


def reconstruct(t,err,base=None,norm=False,rotation=False,kinematic='none',bias_delta=False):
    if base is None:base=all_group_config()
    err=np.asarray(err,dtype=np.int64);channels=err.shape[1];past=[];recovered=np.zeros_like(err)
    for i,row in enumerate(err):
        prev=past[-1] if past else None
        rp=qrot_pred(past[-1][3:7],past[-2][3:7]) if rotation and len(past)>1 else None
        order=[7,8,9,4,5,6,3,0,1,2,10,11,12,13,14,15] if channels==16 else [4,5,6,3,0,1,2]
        if kinematic=='v_from_p':order=[0,1,2,4,5,6,3,7,8,9,10,11,12,13,14,15]
        cur=np.zeros(channels,dtype=np.int64)
        for ch in order:
            if ch<3:grp='position'
            elif ch<7:grp='quaternion'
            elif ch<10:grp='velocity'
            else:grp='bias'
            pred=temporal_pred(past,ch,'delta' if bias_delta and grp=='bias' else base[grp])
            if rp is not None and 3<=ch<7:pred=int(rp[ch-3])
            if norm and ch==3 and prev is not None:pred=norm_pred(cur[4:7],int(prev[3]))
            if prev is not None and kinematic=='p_from_v' and ch<3:
                delta=int(t[i]-t[i-1]);pred=int(prev[ch])+round_div_signed(delta*(int(prev[ch+7])+int(cur[ch+7])),2_000_000_000)
            if prev is not None and kinematic=='v_from_p' and 7<=ch<10:
                delta=int(t[i]-t[i-1]);pred=round_div_signed((int(cur[ch-7])-int(prev[ch-7]))*1_000_000_000,delta)
            cur[ch]=int(row[ch])+pred
        if np.any(cur>2**31-1) or np.any(cur<-(2**31)):
            raise OverflowError('decoded i32 overflow')
        recovered[i]=cur;past.append(cur)
    return recovered.astype('<i4')


def uvarint(v):
    if v<0:raise ValueError('negative')
    b=bytearray()
    while v>127:b.append((v&127)|128);v>>=7
    b.append(v);return b


def read_uvarint(data,pos):
    v=0
    for j in range(11):
        if pos>=len(data):raise ValueError('truncated varint')
        k=data[pos];pos+=1;v|=(k&127)<<(7*j)
        if k<128:return v,pos
    raise ValueError('overlong varint')


def zz(x):return (2*x if x>=0 else -2*x-1)
def unzz(x):return x//2 if x%2==0 else -(x+1)//2


def serialize_times(t):
    out=bytearray(struct.pack('<q',int(t[0])))
    for dt in np.diff(t):out.extend(uvarint(int(dt)))
    return bytes(out)


def deserialize_times(data,n):
    if len(data)<8:raise ValueError('timestamps truncated')
    t=np.empty(n,dtype=np.int64);t[0]=struct.unpack('<q',data[:8])[0];pos=8
    for i in range(1,n):
        dt,pos=read_uvarint(data,pos)
        if dt<=0:raise ValueError('nonmonotonic timestamps')
        t[i]=int(t[i-1])+int(dt)
    if pos!=len(data):raise ValueError('trailing time bytes')
    return t


def serialize_residuals(arr):
    """Cross-channel predictor residuals, column major zigzag varint."""
    out=bytearray()
    for col in arr.T:
        for x in col:out.extend(uvarint(zz(int(x))))
    return bytes(out)


def deserialize_residuals(data,n,c):
    arr=np.zeros((n,c),dtype=np.int64);pos=0
    for col in range(c):
        for i in range(n):
            v,pos=read_uvarint(data,pos);arr[i,col]=unzz(v)
    if pos!=len(data):raise ValueError('trailing residual payload')
    return arr


def encode(t,source,base=None,norm=False,rotation=False,kinematic='none',bias_delta=False,backend='brotli'):
    if base is None:base=all_group_config()
    if len(t)<1:raise ValueError('empty input')
    t=np.asarray(t,dtype=np.int64);source=np.asarray(source,dtype='<i4')
    if not np.all(np.diff(t)>0):raise ValueError('timestamps not increasing')
    res=build_residuals(t,source,base,norm,rotation,kinematic,bias_delta)
    dat=serialize_residuals(res);times=serialize_times(t)
    if backend=='brotli':payload=brotli.compress(dat,quality=9);be=0
    elif backend=='zlib':payload=zlib.compress(dat,9);be=1
    else:raise ValueError(backend)
    kin={'none':0,'p_from_v':1,'v_from_p':2}[kinematic]
    flags=int(norm)+(int(rotation)<<1)+(int(bias_delta)<<2)+(kin<<3)+(be<<5)
    meta=(flags<<8)|_index_cfg(base)
    crc=zlib.crc32(source.tobytes())
    return struct.pack(HEADER,MAGIC,1,len(source),source.shape[1],meta,len(times),crc)+times+payload


def decode(blob):
    if len(blob)<HEADER_SIZE:raise ValueError('short input')
    magic,version,n,c,meta,tlen,crc=struct.unpack_from(HEADER,blob,0)
    if magic!=MAGIC or version!=1 or not (0<n<=2_000_000) or c not in (7,16):raise ValueError('format/version')
    flags=meta>>8;cfg=_decode_cfg(meta&255)
    norm=bool(flags&1);rotation=bool(flags&2);bias_delta=bool(flags&4)
    kin={0:'none',1:'p_from_v',2:'v_from_p'}.get((flags>>3)&3)
    if kin is None:raise ValueError('kinematic metadata')
    if HEADER_SIZE+tlen>len(blob):raise ValueError('short timestamp payload')
    t=deserialize_times(blob[HEADER_SIZE:HEADER_SIZE+tlen],n)
    res_comp=blob[HEADER_SIZE+tlen:]
    res=deserialize_residuals(zlib.decompress(res_comp) if flags&(1<<5) else brotli.decompress(res_comp),n,c)
    arr=reconstruct(t,res,cfg,norm,rotation,kin,bias_delta)
    if zlib.crc32(arr.tobytes())!=crc:raise ValueError('CRC failure')
    return t,arr


def synthetic(n=320,seed=20261118):
    rng=np.random.default_rng(seed)
    t=(1400000000000000000+np.arange(n,dtype=np.int64)*5_000_000)
    out=np.zeros((n,16),dtype=np.int32)
    r=np.array([1.,0.,0.,0.]);p=np.array([0.,0.,0.]);v=np.array([1.,0.,0.])
    for i in range(n):
        w=np.array([.003,.006,-.001])+rng.normal(0,.0002,3)
        dq=np.r_[np.sqrt(1-np.dot(w,w)),w]
        r=_qmul(r,dq);r/=np.linalg.norm(r)
        v=v+rng.normal(0,.0001,3)
        p=p+v*.005
        out[i,0:3]=np.rint(p*SCALE)
        out[i,3:7]=np.rint(r*SCALE)
        out[i,7:10]=np.rint(v*SCALE)
        out[i,10:16]=np.rint(np.array([0.02,0.01,-.02,.002,.001,-.003])*SCALE+rng.normal(0,5,6)).astype(np.int32)
    return t,out


if __name__=='__main__':
    t,arr=synthetic(1000)
    for norm in (False,True):
      for rotation in (False,True):
        for kin in ('none','p_from_v','v_from_p'):
          for bias in (False,True):
            b=encode(t,arr,norm=norm,rotation=rotation,kinematic=kin,bias_delta=bias)
            tt,rec=decode(b)
            assert np.array_equal(tt,t) and np.array_equal(rec,arr)
            print(norm,rotation,kin,bias,len(b),flush=True)
