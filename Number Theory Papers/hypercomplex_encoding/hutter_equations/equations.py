"""Causal integer equations for byte prediction; all operations are mod 256.

A model explicitly predicts the next *byte*, not a label for the previous context.
All formula inputs are past decoded bytes, and the residual is always reversible.
No equations here should be confused with algebraic identities universally obeyed by
English text. Validity is measured, not assumed.
"""
from __future__ import annotations
import struct
import zlib
import numpy as np
import brotli

SEED_LENGTH = 16
HEADER = struct.Struct('<4sBBII')
MAGIC = b'EQD2'

# All formula names and families are declared before measuring validation/test.
FORMULAS = (
    ('literal_space', 'literal'), ('literal_e', 'literal'),
    ('copy_1', 'copy'), ('copy_2', 'copy'), ('copy_4', 'copy'),
    ('copy_8', 'copy'), ('copy_16', 'copy'),
    ('affine_2x_minus_y', 'affine'), ('affine_xyz', 'affine'),
    ('affine_x4_y5', 'affine'), ('affine_x8_y9', 'affine'),
    ('affine_x_minus_y_plus_z', 'affine'),
    ('xor_1_2', 'xor'), ('xor_1_2_3', 'xor'), ('xor_1_4', 'xor'),
    ('quadratic_12', 'quadratic'), ('quadratic_12_plus_3', 'quadratic'),
    ('quadratic_12_minus_34', 'quadratic'),
    ('complex_real', 'complex'), ('complex_imag', 'complex'),
    ('complex_real_hi', 'complex'), ('complex_imag_hi', 'complex'),
    ('quaternion_scalar', 'quaternion'), ('quaternion_i', 'quaternion'),
    ('quaternion_j', 'quaternion'), ('quaternion_k', 'quaternion'),
    ('quaternion_scalar_hi', 'quaternion'), ('quaternion_i_hi', 'quaternion'),
    ('quaternion_j_hi', 'quaternion'), ('quaternion_k_hi', 'quaternion'),
)
NAMES = tuple(x[0] for x in FORMULAS)
FAMILIES = tuple(x[1] for x in FORMULAS)
IDX = {name: i for i, name in enumerate(NAMES)}


def prev_class(p: int) -> int:
    """Eight frozen byte categories, determined from the previous decoded byte."""
    if 97 <= p <= 122: return 0  # lowercase
    if 65 <= p <= 90: return 1   # uppercase
    if 48 <= p <= 57: return 2   # digit
    if p in (9, 10, 13, 32): return 3
    if p in (60, 62, 38, 47, 61): return 4  # XML punctuation
    if p < 32 or p == 127: return 5
    if p >= 128: return 6
    return 7


def _formula(n: str, d: list[int]) -> int:
    """d[0] is preceding byte, d[15] is 16th preceding byte."""
    x, y, z, u = d[0:4]
    if n == 'literal_space': return 32
    if n == 'literal_e': return 101
    if n.startswith('copy_'): return d[int(n.split('_')[1])-1]
    if n == 'affine_2x_minus_y': return 2*x-y
    if n == 'affine_xyz': return x+y-z
    if n == 'affine_x4_y5': return x+d[3]-d[4]
    if n == 'affine_x8_y9': return x+d[7]-d[8]
    if n == 'affine_x_minus_y_plus_z': return x-y+u
    if n == 'xor_1_2': return x ^ y
    if n == 'xor_1_2_3': return x ^ y ^ z
    if n == 'xor_1_4': return x ^ u
    if n == 'quadratic_12': return x*y
    if n == 'quadratic_12_plus_3': return x*y+z
    if n == 'quadratic_12_minus_34': return x*y-z*u
    if n.startswith('complex_'):
        a, b, c, dd = (d[j]-128 for j in range(4))
        real, imag = a*c-b*dd, a*dd+b*c
        v = real if 'real' in n else imag
        return v >> 8 if n.endswith('_hi') else v
    if n.startswith('quaternion_'):
        a, b, c, e = (d[j]-128 for j in range(4))
        f, g, h, k = (d[j]-128 for j in range(4,8))
        products = (a*f-b*g-c*h-e*k,
                    a*g+b*f+c*k-e*h,
                    a*h-b*k+c*f+e*g,
                    a*k+b*h-c*g+e*f)
        axis = n[len('quaternion_'):].removesuffix('_hi')
        v = products[('scalar', 'i', 'j', 'k').index(axis)]
        return v >> 8 if n.endswith('_hi') else v
    raise ValueError('unknown equation: '+n)


def predict_one(model_id: int, history: bytes | bytearray) -> int:
    if model_id < 0 or model_id >= len(NAMES):
        raise ValueError('invalid equation ID')
    d = [history[-1-j] if j < len(history) else 0 for j in range(SEED_LENGTH)]
    return _formula(NAMES[model_id], d) & 255


def predictions(data: bytes) -> dict[str, np.ndarray]:
    """Vectorized predictions. No data[i:] contributes to any prediction[i]."""
    s = np.frombuffer(data, dtype=np.uint8).astype(np.int64)
    def lag(k):
        a = np.zeros(s.size, dtype=np.int64)
        a[k:] = s[:-k]
        return a
    d = [lag(k) for k in range(1, SEED_LENGTH+1)]
    x, y, z, u = d[:4]
    a, b, c, e, f, g, h, k = (v-128 for v in d[:8])
    q = (a*f-b*g-c*h-e*k,
         a*g+b*f+c*k-e*h,
         a*h-b*k+c*f+e*g,
         a*k+b*h-c*g+e*f)
    re = (x-128)*(z-128) - (y-128)*(u-128)
    im = (x-128)*(u-128) + (y-128)*(z-128)
    formulas = {
        'literal_space': np.full(s.size,32,dtype=np.int64),
        'literal_e':np.full(s.size,101,dtype=np.int64),
        'copy_1':x,'copy_2':y,'copy_4':u,'copy_8':d[7],'copy_16':d[15],
        'affine_2x_minus_y':2*x-y,
        'affine_xyz':x+y-z,
        'affine_x4_y5':x+d[3]-d[4],
        'affine_x8_y9':x+d[7]-d[8],
        'affine_x_minus_y_plus_z':x-y+u,
        'xor_1_2':np.bitwise_xor(x,y),
        'xor_1_2_3':np.bitwise_xor(np.bitwise_xor(x,y),z),
        'xor_1_4':np.bitwise_xor(x,u),
        'quadratic_12':x*y,
        'quadratic_12_plus_3':x*y+z,
        'quadratic_12_minus_34':x*y-z*u,
        'complex_real':re,'complex_imag':im,
        'complex_real_hi':re >>8,'complex_imag_hi':im>>8,
    }
    for ax, arr in zip(('scalar','i','j','k'),q):
        formulas['quaternion_'+ax]=arr
        formulas['quaternion_'+ax+'_hi']=arr>>8
    assert set(formulas)==set(NAMES)
    return {name:(arr & 255).astype(np.uint8) for name,arr in formulas.items()}


def _raw_residual(data: bytes, equation_id: int, class_mapping: bytes | None = None) -> bytes:
    """Remainder byte stream; seed held literally, residuals use causal predictions."""
    if len(data) < SEED_LENGTH: return b''
    p = predictions(data)
    orig = np.frombuffer(data,dtype=np.uint8)
    if class_mapping is None:
        pred = p[NAMES[equation_id]]
    else:
        if len(class_mapping)!=8: raise ValueError('8 class rules required')
        prev=np.empty(len(orig),dtype=np.uint8);prev[0]=0;prev[1:]=orig[:-1]
        cls=np.fromiter((prev_class(int(v)) for v in prev),dtype=np.uint8,count=len(orig))
        pred = np.zeros(len(orig), dtype=np.uint8)
        for klass, model in enumerate(class_mapping):
            pred[cls==klass] = p[NAMES[model]][cls==klass]
    residual=((orig.astype(np.int16)-pred.astype(np.int16))&255).astype(np.uint8)
    return residual[SEED_LENGTH:].tobytes()


def encode(data:bytes, equation_id:int=0, class_mapping:bytes | None=None, quality:int=11)->bytes:
    """Fully standalone lossless framed compressed residual; decoder needs no source text."""
    if not 0<=equation_id<len(NAMES): raise ValueError('bad equation index')
    if class_mapping is not None and (len(class_mapping)!=8 or any(i>=len(NAMES) for i in class_mapping)):
        raise ValueError('bad model map')
    seed=data[:min(len(data),SEED_LENGTH)]
    flags=1 if class_mapping is not None else 0
    header=HEADER.pack(MAGIC, flags, equation_id, len(data),zlib.crc32(data))
    meta=bytes(class_mapping) if class_mapping else b''
    return header+meta+seed+brotli.compress(_raw_residual(data,equation_id,class_mapping),quality=quality)


def decode(frame:bytes)->bytes:
    if len(frame)<HEADER.size: raise ValueError('truncated header')
    magic,flags,mid,n,crc=HEADER.unpack_from(frame)
    if magic!=MAGIC or flags not in (0,1) or mid>=len(NAMES):raise ValueError('not equation-coded')
    if n>50_000_000:raise ValueError('oversized frame')
    p=HEADER.size
    clmap=None
    if flags:
        if p+8>len(frame): raise ValueError('truncated map')
        clmap=frame[p:p+8];p+=8
        if any(i>=len(NAMES) for i in clmap):raise ValueError('invalid class model')
    seedsize=min(n,SEED_LENGTH)
    if p+seedsize>len(frame):raise ValueError('truncated seed')
    out=bytearray(frame[p:p+seedsize]); p+=seedsize
    residual=brotli.decompress(frame[p:])
    if len(residual)!=n-seedsize:raise ValueError('wrong residual length')
    for r in residual:
        eqid=clmap[prev_class(out[-1])] if clmap is not None else mid
        v=(predict_one(eqid,out)+r)&255
        out.append(v)
    if len(out)!=n or zlib.crc32(out)!=crc:raise ValueError('CRC mismatch')
    return bytes(out)
