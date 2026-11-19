"""Byte-exact standalone hybrid and ablation codecs for the SOTA comparison.

Only the hypercomplex variant invokes Cayley-Dickson basis multiplication.
All framed sizes include headers, block modes, and lengths; no external side data.
"""
from __future__ import annotations

import ctypes
import ctypes.util
import lzma
import struct
import zlib
from collections import Counter
from dataclasses import dataclass

import brotli
import numpy as np
from numba import njit

from hypercomplex_cascade_test import cd_cascade, cd_cascade_inverse_checked, varint, unvarint, zigzag, unzigzag

_HEADER = struct.Struct('<4sBBIIQI')  # magic, version, variant, chunk bytes, D, original bytes, CRC32
MAGIC = b'HYHC'
VERSION = 1
VARIANT = {'hypercomplex': 0, 'walsh': 1, 'xor': 2, 'no_spectral': 3, 'basis_index': 4}
REVERSE = {v: k for k, v in VARIANT.items()}
RAW, ZSTD3, BROTLI5, SPECTRAL, FILL0, FILL1, XOR_ZSTD, XOR_BROTLI, BASIS_INDEX = range(9)
MODE_NAMES = {RAW:'raw', ZSTD3:'zstd3', BROTLI5:'brotli5', SPECTRAL:'spectral', FILL0:'fill0', FILL1:'fill1', XOR_ZSTD:'xor_zstd3', XOR_BROTLI:'xor_brotli5', BASIS_INDEX:'basis_index'}


class Zstd:
    """libzstd C-ABI wrapper: levels, trained dictionaries, and exact decoding.

    Avoids subprocess launch overhead when benchmarking many short blocks. Requires
    libzstd 1.4+ in the environment (Ubuntu's libzstd1 satisfies this).
    """
    def __init__(self):
        path = ctypes.util.find_library('zstd')
        if not path:
            raise RuntimeError('libzstd unavailable; install libzstd1 to run benchmark')
        self.lib = ctypes.CDLL(path)
        L = self.lib
        v = ctypes.c_void_p
        z = ctypes.c_size_t
        L.ZSTD_versionString.argtypes = []
        L.ZSTD_versionString.restype = ctypes.c_char_p
        L.ZSTD_compressBound.argtypes = [z]
        L.ZSTD_compressBound.restype = z
        L.ZSTD_compress.argtypes = [v,z,v,z,ctypes.c_int]
        L.ZSTD_compress.restype = z
        L.ZSTD_decompress.argtypes = [v,z,v,z]
        L.ZSTD_decompress.restype = z
        L.ZSTD_isError.argtypes = [z]
        L.ZSTD_isError.restype = ctypes.c_uint
        L.ZSTD_getErrorName.argtypes = [z]
        L.ZSTD_getErrorName.restype = ctypes.c_char_p
        L.ZSTD_createCCtx.argtypes = []
        L.ZSTD_createCCtx.restype = v
        L.ZSTD_freeCCtx.argtypes = [v]
        L.ZSTD_freeCCtx.restype = z
        L.ZSTD_createDCtx.argtypes = []
        L.ZSTD_createDCtx.restype = v
        L.ZSTD_freeDCtx.argtypes = [v]
        L.ZSTD_freeDCtx.restype = z
        L.ZSTD_compress_usingDict.argtypes = [v,v,z,v,z,v,z,ctypes.c_int]
        L.ZSTD_compress_usingDict.restype = z
        L.ZSTD_decompress_usingDict.argtypes = [v,v,z,v,z,v,z]
        L.ZSTD_decompress_usingDict.restype = z
        L.ZDICT_trainFromBuffer.argtypes = [v,z,v,ctypes.POINTER(z),ctypes.c_uint]
        L.ZDICT_trainFromBuffer.restype = z
        L.ZDICT_isError.argtypes = [z]
        L.ZDICT_isError.restype = ctypes.c_uint
        L.ZDICT_getErrorName.argtypes = [z]
        L.ZDICT_getErrorName.restype = ctypes.c_char_p
        self.version = L.ZSTD_versionString().decode()

    def _check(self, n, dictionary=False):
        L = self.lib
        if (L.ZDICT_isError(n) if dictionary else L.ZSTD_isError(n)):
            fn = L.ZDICT_getErrorName if dictionary else L.ZSTD_getErrorName
            raise RuntimeError(fn(n).decode())
        return n

    def compress(self, data: bytes, level: int=3, dictionary: bytes|None=None) -> bytes:
        L = self.lib
        src = ctypes.create_string_buffer(data)
        dest = ctypes.create_string_buffer(L.ZSTD_compressBound(len(data)))
        if dictionary is None:
            n = L.ZSTD_compress(dest,len(dest),src,len(data),int(level))
        else:
            dic = ctypes.create_string_buffer(dictionary)
            ctx = L.ZSTD_createCCtx()
            if not ctx: raise MemoryError('ZSTD_createCCtx')
            try: n = L.ZSTD_compress_usingDict(ctx,dest,len(dest),src,len(data),dic,len(dictionary),int(level))
            finally: L.ZSTD_freeCCtx(ctx)
        return dest.raw[:self._check(n)]

    def decompress(self, blob: bytes, expected: int, dictionary: bytes|None=None) -> bytes:
        if not 0 <= expected <= 1<<31:
            raise ValueError('unreasonable decompressed length')
        L = self.lib
        src = ctypes.create_string_buffer(blob)
        dest = ctypes.create_string_buffer(max(1,expected))
        if dictionary is None:
            n = L.ZSTD_decompress(dest,expected,src,len(blob))
        else:
            dic = ctypes.create_string_buffer(dictionary)
            ctx = L.ZSTD_createDCtx()
            if not ctx: raise MemoryError('ZSTD_createDCtx')
            try: n = L.ZSTD_decompress_usingDict(ctx,dest,expected,src,len(blob),dic,len(dictionary))
            finally: L.ZSTD_freeDCtx(ctx)
        self._check(n)
        if n != expected:
            raise ValueError(f'zstd wrong length: {n} != {expected}')
        return dest.raw[:n]

    def train(self, samples: list[bytes], capacity: int=4096) -> bytes:
        # Training data never includes held-out evaluation bytes.
        samples = [x for x in samples if len(x)>=16]
        if len(samples)<8: raise ValueError('need >=8 calibration samples for zstd dictionary')
        data = b''.join(samples)
        lengths = (ctypes.c_size_t * len(samples))(*(len(x) for x in samples))
        dest = ctypes.create_string_buffer(capacity)
        src = ctypes.create_string_buffer(data)
        n = self.lib.ZDICT_trainFromBuffer(dest,capacity,src,lengths,len(samples))
        return dest.raw[:self._check(n,dictionary=True)]


@njit(cache=True)
def fwht(vec):
    x = np.empty(len(vec),np.int64)
    for i in range(len(vec)):
        x[i] = vec[i]
    stride = 1
    while stride < len(x):
        for base in range(0,len(x),2*stride):
            for j in range(stride):
                a=x[base+j]; b=x[base+j+stride]
                x[base+j]=a+b
                x[base+j+stride]=a-b
        stride *= 2
    return x


def inverse_fwht(vec):
    x=fwht(vec)
    if np.any(x % len(vec)):
        raise ValueError('invalid Walsh spectrum')
    return x//len(vec)


def _xor_predict(data: bytes) -> bytes:
    if not data: return b''
    return bytes([data[0]]+[a^b for a,b in zip(data[1:],data[:-1])])


def _xor_inverse(data: bytes) -> bytes:
    out = bytearray()
    last = 0
    for x in data:
        last ^= x
        out.append(last)
    return bytes(out)


def _spectral_encode(bits: np.ndarray, D: int, variant: str) -> tuple[bytes|None,bytes|None,int]:
    # Exactly D-1 imaginary bits per block, including unaligned byte boundaries.
    if len(bits)>D-1:raise ValueError('block exceeds its hypercomplex dimension')
    coeff=np.full(D,-1,dtype=np.int8)
    coeff[0]=1
    coeff[1:len(bits)+1]=2*bits.astype(np.int8)-1
    spectrum=fwht(coeff) if variant=='walsh' else cd_cascade(coeff)
    idx=np.flatnonzero(spectrum)
    # Do not serialize a dense spectrum, because raw and conventional modes remain available.
    if len(idx)>max(1,D//16):
        return None,None,len(idx)
    encoded=bytearray(varint(len(idx)))
    prev=0
    for coordinate in idx:
        i=int(coordinate)
        encoded+=varint(i-prev)
        encoded+=varint(zigzag(int(spectrum[i])))
        prev=i
    basis=None
    if len(idx)==1 and abs(int(spectrum[idx[0]]))==D:
        # A specialized family-specific basis-index code, not a new entropy bound.
        basis=varint(int(idx[0]))+bytes([int(spectrum[idx[0]])<0])
    return bytes(encoded),basis,len(idx)


def _spectral_decode(blob:bytes,D:int,variant:str,basis=False)->np.ndarray:
    spectrum=np.zeros(D,dtype=np.int64)
    pos=0
    if basis:
        idx,pos=unvarint(blob,pos)
        if pos+1!=len(blob) or idx>=D or blob[pos] not in (0,1):
            raise ValueError('invalid basis-index block')
        spectrum[idx]= -D if blob[pos] else D
    else:
        count,pos=unvarint(blob,pos)
        if count>D:raise ValueError('invalid spectrum count')
        idx=0
        for k in range(count):
            delta,pos=unvarint(blob,pos)
            val,pos=unvarint(blob,pos)
            idx+=delta
            if idx>=D or (k>0 and delta==0):raise ValueError('invalid spectrum coordinate')
            spectrum[idx]=unzigzag(val)
        if pos!=len(blob):raise ValueError('extra spectral bytes')
    coeff=inverse_fwht(spectrum) if variant=='walsh' else cd_cascade_inverse_checked(spectrum)
    if coeff[0]!=1 or np.any((coeff[1:]!=-1)&(coeff[1:]!=1)):
        raise ValueError('spectrum not a valid bit encoding')
    return ((coeff[1:]+1)//2).astype(np.uint8)


@dataclass(frozen=True)
class HybridConfig:
    block_bytes:int=127
    variant:str='hypercomplex'
    def __post_init__(self):
        if not (self.block_bytes in (8,16,32,64,128,256,512,1024)):
            raise ValueError('block_bytes must be a power-of-two multiple of 8 (D=8*bytes)')
        if self.variant not in VARIANT:raise ValueError('unknown variant')
    @property
    def dimension(self):
        return self.block_bytes*8


def hybrid_encode(data:bytes,cfg:HybridConfig,zstd:Zstd,stats:dict|None=None)->bytes:
    """Encode the *continuous* file bitstream in exactly D-1-bit chunks.

    This is important: rounding each chunk down to whole bytes would destroy
    the complete basis-index spectral codewords used in earlier research.
    """
    D=cfg.dimension
    out=bytearray(_HEADER.pack(MAGIC,VERSION,VARIANT[cfg.variant],cfg.block_bytes,D,len(data),zlib.crc32(data)))
    counts=Counter()
    bits=np.unpackbits(np.frombuffer(data,dtype=np.uint8))
    per=D-1
    for offset in range(0,len(bits),per):
        chunk=bits[offset:offset+per]
        aligned=np.zeros(D,dtype=np.uint8)  # adds exactly one pad bit per full block
        aligned[:len(chunk)]=chunk
        block=np.packbits(aligned).tobytes()
        options=[(RAW,block)]
        if not np.any(chunk):options.append((FILL0,b''))
        if bool(np.all(chunk==1)):options.append((FILL1,b''))
        options.append((ZSTD3,zstd.compress(block,3)))
        options.append((BROTLI5,brotli.compress(block,quality=5)))
        if cfg.variant=='xor':
            predicted=_xor_predict(block)
            options.extend([(XOR_ZSTD,zstd.compress(predicted,3)),(XOR_BROTLI,brotli.compress(predicted,quality=5))])
        if cfg.variant in ('hypercomplex','walsh','basis_index'):
            transform='walsh' if cfg.variant=='walsh' else 'hypercomplex'
            spectral,basis,_=_spectral_encode(chunk,D,transform)
            if spectral is not None and cfg.variant!='basis_index':options.append((SPECTRAL,spectral))
            if basis is not None:options.append((BASIS_INDEX,basis))
        mode,payload=min(options,key=lambda item:1+len(varint(len(item[1])))+len(item[1]))
        counts[MODE_NAMES[mode]]+=1
        out.append(mode)
        out.extend(varint(len(payload)))
        out.extend(payload)
    if stats is not None:
        stats.update({'modes':dict(counts),'blocks':sum(counts.values()),'dimension':D,
                      'variant':cfg.variant,'block_bytes':cfg.block_bytes})
    return bytes(out)


def hybrid_decode(blob:bytes,zstd:Zstd)->bytes:
    if len(blob)<_HEADER.size:raise ValueError('truncated hybrid header')
    magic,version,variant,block_bytes,D,original_len,crc=_HEADER.unpack_from(blob)
    if magic!=MAGIC or version!=VERSION or variant not in REVERSE:
        raise ValueError('wrong hybrid format')
    cfg=HybridConfig(block_bytes,REVERSE[variant])
    if cfg.dimension!=D:raise ValueError('dimension does not match block size')
    if original_len>1<<30:raise ValueError('input too long for benchmark decoder')
    pos=_HEADER.size
    bits_out=np.empty(original_len*8,dtype=np.uint8)
    per=D-1
    blocks=(original_len*8+per-1)//per
    for block_idx in range(blocks):
        if pos>=len(blob):raise ValueError('truncated hybrid stream')
        mode=blob[pos];pos+=1
        count,pos=unvarint(blob,pos)
        if count>len(blob)-pos:raise ValueError('invalid compressed block length')
        payload=blob[pos:pos+count];pos+=count
        offset=block_idx*per
        want=min(per,original_len*8-offset)
        if mode==FILL0:
            if payload:raise ValueError('nonempty zero-fill')
            bits=np.zeros(want,dtype=np.uint8)
        elif mode==FILL1:
            if payload:raise ValueError('nonempty one-fill')
            bits=np.ones(want,dtype=np.uint8)
        elif mode in (RAW,ZSTD3,BROTLI5,XOR_ZSTD,XOR_BROTLI):
            if mode==RAW:part=payload
            elif mode==ZSTD3:part=zstd.decompress(payload,block_bytes)
            elif mode==BROTLI5:part=brotli.decompress(payload)
            elif mode==XOR_ZSTD:
                if cfg.variant!='xor':raise ValueError('incorrect XOR mode')
                part=_xor_inverse(zstd.decompress(payload,block_bytes))
            else:
                if mode==XOR_BROTLI:
                    if cfg.variant!='xor':raise ValueError('incorrect XOR mode')
                    part=_xor_inverse(brotli.decompress(payload))
                else:part=brotli.decompress(payload)
            if len(part)!=block_bytes:raise ValueError('wrong decompressed block length')
            bits=np.unpackbits(np.frombuffer(part,dtype=np.uint8))[:want]
        elif mode==SPECTRAL:
            if cfg.variant not in ('hypercomplex','walsh'):
                raise ValueError('invalid spectral variant')
            bits=_spectral_decode(payload,D,cfg.variant)[:want]
        elif mode==BASIS_INDEX:
            if cfg.variant not in ('hypercomplex','walsh','basis_index'):
                raise ValueError('invalid basis mode')
            bits=_spectral_decode(payload,D,'walsh' if cfg.variant=='walsh' else 'hypercomplex',basis=True)[:want]
        else:raise ValueError('invalid block mode')
        bits_out[offset:offset+want]=bits
    if pos!=len(blob):raise ValueError('trailing bytes in framed payload')
    result=np.packbits(bits_out).tobytes()
    if zlib.crc32(result)!=crc:raise ValueError('CRC mismatch')
    return result


_DICT_HEADER=struct.Struct('<4sIII') # magic, original len, dictionary size, crc

def zstd_dictionary_encode(data:bytes,zstd:Zstd,dictionary:bytes,level:int=9)->bytes:
    payload=zstd.compress(data,level,dictionary)
    return _DICT_HEADER.pack(b'ZDHC',len(data),len(dictionary),zlib.crc32(data))+dictionary+payload


def zstd_dictionary_decode(blob:bytes,zstd:Zstd)->bytes:
    if len(blob)<_DICT_HEADER.size:raise ValueError('truncated dictionary frame')
    magic,total,dlen,crc=_DICT_HEADER.unpack_from(blob)
    if magic!=b'ZDHC' or dlen>len(blob)-_DICT_HEADER.size:raise ValueError('invalid dictionary frame')
    d=blob[_DICT_HEADER.size:_DICT_HEADER.size+dlen]
    output=zstd.decompress(blob[_DICT_HEADER.size+dlen:],total,d)
    if zlib.crc32(output)!=crc:raise ValueError('CRC mismatch')
    return output


def compress_named(data:bytes,name:str,zstd:Zstd,dictionary:bytes|None=None)->bytes:
    if name=='raw':return data
    if name.startswith('zstd') and name[4:].isdigit():return zstd.compress(data,int(name[4:]))
    if name.startswith('brotli') and name[6:].isdigit():return brotli.compress(data,quality=int(name[6:]))
    if name.startswith('lzma') and name[4:].isdigit():return lzma.compress(data,preset=int(name[4:]))
    if name.startswith('zlib') and name[4:].isdigit():return zlib.compress(data,int(name[4:]))
    if name.startswith('zstd_dict') and name[9:].isdigit():
        if not dictionary:raise ValueError('need calibration-trained dictionary')
        return zstd_dictionary_encode(data,zstd,dictionary,int(name[9:]))
    raise ValueError('unknown compressor: '+name)


def decompress_named(blob:bytes,name:str,zstd:Zstd,expected:int)->bytes:
    if name=='raw':return blob
    if name.startswith('zstd_dict'):return zstd_dictionary_decode(blob,zstd)
    if name.startswith('zstd'):return zstd.decompress(blob,expected)
    if name.startswith('brotli'):return brotli.decompress(blob)
    if name.startswith('lzma'):return lzma.decompress(blob)
    if name.startswith('zlib'):return zlib.decompress(blob)
    raise ValueError('unknown compressor: '+name)

# Independent family-specific *format* comparator: index the CD basis family.
# Recognition still needs the same transform, but only a known codebook index
# is transmitted, without general spectral (index,value) pairs or block selectors.
_BASIS_HDR=struct.Struct('<4sIIQI')  # magic, D, full D-1-bit blocks, raw byte len, crc


def basis_codebook_encode(data:bytes,D:int)->bytes|None:
    if D not in (64,256,1024,4096,8192):
        raise ValueError('basis codebook dimension outside experimental range')
    bits=np.unpackbits(np.frombuffer(data,dtype=np.uint8))
    per=D-1
    full=len(bits)//per
    out=bytearray(_BASIS_HDR.pack(b'BIDX',D,full,len(data),zlib.crc32(data)))
    for k in range(full):
        part=bits[k*per:(k+1)*per]
        spectral,basis,_=_spectral_encode(part,D,'hypercomplex')
        if basis is None:return None
        out.extend(basis)
    tail=bits[full*per:]
    # The tail is literally part of the standalone format and charged in size.
    out.extend(np.packbits(tail).tobytes())
    return bytes(out)


def basis_codebook_decode(blob:bytes)->bytes:
    if len(blob)<_BASIS_HDR.size:raise ValueError('incomplete basis-codebook header')
    magic,D,full,total,crc=_BASIS_HDR.unpack_from(blob)
    if magic!=b'BIDX' or D not in (64,256,1024,4096,8192):
        raise ValueError('invalid basis-codebook signature')
    if total>1<<28 or full>total*8//(D-1):
        raise ValueError('invalid basis-codebook length')
    per=D-1
    output=np.zeros(total*8,dtype=np.uint8)
    pos=_BASIS_HDR.size
    for k in range(full):
        j,pos=unvarint(blob,pos)
        if pos>=len(blob) or j>=D or blob[pos] not in (0,1):
            raise ValueError('bad basis-codeword index/sign')
        sign=blob[pos];pos+=1
        value=np.zeros(D,np.int64);value[j]= -D if sign else D
        decoded=cd_cascade_inverse_checked(value)
        if decoded[0]!=1 or np.any(np.abs(decoded[1:])!=1):
            raise ValueError('basis index lacks fixed +1 real coefficient')
        output[k*per:(k+1)*per]=((decoded[1:]+1)//2).astype(np.uint8)
    tail_bits=total*8-full*per
    raw_bytes=(tail_bits+7)//8
    if pos+raw_bytes!=len(blob):raise ValueError('basis tail length mismatch')
    if raw_bytes:
        output[full*per:]=np.unpackbits(np.frombuffer(blob[pos:],np.uint8))[:tail_bits]
    decoded=np.packbits(output).tobytes()
    if zlib.crc32(decoded)!=crc:raise ValueError('basis-codebook checksum failed')
    return decoded
