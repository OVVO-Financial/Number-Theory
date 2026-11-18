#!/usr/bin/env python3
"""Encode/decode user files with the experimental actual hypercomplex algebra codecs."""
from __future__ import annotations
import argparse
from pathlib import Path
import zlib
import brotli
from hypercomplex_highdim_test import encode_unit,decode_unit
from hypercomplex_cascade_test import encode_sparse_cascade,decode_sparse_cascade


def encode_file(data,mode,dimension,bipolar,compress):
    if mode=='cascade':
        body,stats=encode_sparse_cascade(data,dimension,bipolar)
    elif mode=='unit':
        body=encode_unit(data,dimension,min(7,dimension-1))
        stats={'unit':1}
    else:raise ValueError('invalid mode')
    if compress=='zlib':return b'HCZ1'+zlib.compress(body,9),stats
    if compress=='brotli':return b'HCB1'+brotli.compress(body,quality=9),stats
    return body,stats


def decode_file(blob):
    if blob.startswith(b'HCZ1'):blob=zlib.decompress(blob[4:])
    elif blob.startswith(b'HCB1'):blob=brotli.decompress(blob[4:])
    if blob.startswith((b'HCSP',b'HCSB')):return decode_sparse_cascade(blob)
    if blob.startswith(b'HCD1'):return decode_unit(blob)
    raise ValueError('unknown hypercomplex codec signature')


def main():
    p=argparse.ArgumentParser(description='Lossless high-dimensional Cayley-Dickson coefficient codec')
    sub=p.add_subparsers(dest='command',required=True)
    en=sub.add_parser('encode')
    en.add_argument('input');en.add_argument('output')
    en.add_argument('--mode',choices=('cascade','unit'),default='cascade')
    en.add_argument('--dimension',type=int,default=4096)
    en.add_argument('--coefficient',choices=('bipolar','binary'),default='bipolar')
    en.add_argument('--compress',choices=('brotli','zlib','none'),default='brotli')
    de=sub.add_parser('decode')
    de.add_argument('input');de.add_argument('output')
    args=p.parse_args()
    input_bytes=Path(args.input).read_bytes()
    if args.command=='encode':
        if args.dimension<2 or args.dimension&(args.dimension-1):p.error('dimension must be a power of two >=2')
        encoded,stats=encode_file(input_bytes,args.mode,args.dimension,
                                  args.coefficient=='bipolar',args.compress)
        if decode_file(encoded)!=input_bytes:raise RuntimeError('codec did not roundtrip')
        Path(args.output).write_bytes(encoded)
        print(f'input={len(input_bytes)} encoded={len(encoded)} bytes; exact decode verified; mode_counts={stats}')
    else:
        recovered=decode_file(input_bytes)
        Path(args.output).write_bytes(recovered)
        print(f'encoded={len(input_bytes)} decoded={len(recovered)} bytes; CRC passed')

if __name__=='__main__':main()