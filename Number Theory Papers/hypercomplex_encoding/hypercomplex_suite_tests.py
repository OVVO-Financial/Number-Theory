#!/usr/bin/env python3
"""Independent regression tests for two actual hypercomplex codecs and algebra."""
from __future__ import annotations
import os, random, unittest, zlib
import numpy as np
from hypercomplex_highdim_test import cd_sign, sign_cd_ref, cd_product_dense, encode_unit, decode_unit
from hypercomplex_cascade_test import (cd_cascade, cd_cascade_inverse_checked,
                                       encode_sparse_cascade, decode_sparse_cascade)

class HighDimensionCodecTests(unittest.TestCase):
    def test_quaternion_paper_phi(self):
        for a in range(-3,4):
            for b in range(-3,4):
                x=np.array([a,b],dtype=np.int64)
                e1x=np.array([cd_sign(1,i^1,2)*x[i^1] for i in range(2)],dtype=np.int64)
                self.assertEqual((x+e1x).tolist(),[a-b,a+b])
        for a,b,c,d in ((1,0,1,1),(-1,2,3,-4),(2,3,-5,7)):
            x=np.array([a,b,c,d],dtype=np.int64)
            e1x=np.array([cd_sign(1,i^1,4)*x[i^1] for i in range(4)],dtype=np.int64)
            self.assertEqual((x+e1x).tolist(),[a-b,a+b,c-d,c+d])

    def test_cd_inverse_normalized_gram(self):
        rng=np.random.default_rng(177)
        for n in (2,4,8,16,32,64,256,4096):
            for _ in range(10):
                vec=rng.integers(-7,8,size=n,dtype=np.int16)
                transformed=cd_cascade(vec)
                decoded=cd_cascade_inverse_checked(transformed)
                np.testing.assert_array_equal(vec,decoded)
                self.assertEqual(int(np.dot(transformed,transformed)),
                                 n*int(np.dot(vec.astype(np.int64),vec.astype(np.int64))))

    def test_unit_codec_roundtrips(self):
        rng=random.Random(45)
        for l in (0,1,2,17,127,512,513,4095,4096,16384):
            payload=rng.randbytes(l)
            for n in (256,4096,65536):
                encoded=encode_unit(payload,n,min(7,n-1))
                self.assertEqual(decode_unit(encoded),payload)

    def test_sparse_codec_roundtrips(self):
        rng=random.Random(86)
        messages=[b'',b'\0',b'\xff',b'\xff'*128,b'\0'*1024,
                  bytes([0xB5])*809, rng.randbytes(2048),rng.randbytes(10000)]
        for n in (256,1024,4096,65536):
            for bipolar in (False,True):
                for data in messages:
                    encoded,stats=encode_sparse_cascade(data,n,bipolar)
                    self.assertEqual(decode_sparse_cascade(encoded),data)
                    self.assertEqual(stats['blocks'],stats['raw']+stats['spectral']+stats['fill'])

    def test_sparse_codeword_mode(self):
        n=1024;row=17
        spectrum=np.zeros(n,dtype=np.int64);spectrum[row]=n
        v=cd_cascade_inverse_checked(spectrum)
        if v[0]==-1:v=-v
        self.assertEqual(v[0],1)
        assert np.all(np.abs(v)==1)
        # 1023 bits padded with one trailing zero in one byte.
        payload=np.packbits(((v[1:]+1)//2).astype(np.uint8)).tobytes()
        encoded,stats=encode_sparse_cascade(payload,n,True)
        self.assertEqual(stats['spectral'],1)
        self.assertEqual(decode_sparse_cascade(encoded),payload)

    def test_corruption_detection(self):
        data=bytes(range(256))*2
        for n in (256,1024):
            raw=bytearray(encode_unit(data,n,7))
            # e_1 input maps to e_(1 XOR 7)=e_6; change a known data-bearing bit.
            from hypercomplex_highdim_test import H_SIZE
            raw[H_SIZE]^=0x02
            with self.assertRaises(ValueError): decode_unit(bytes(raw))
            for bipolar in (False,True):
                code,_=encode_sparse_cascade(data,n,bipolar)
                raw=bytearray(code)
                from hypercomplex_cascade_test import HSZ
                self.assertEqual(raw[HSZ],0,'expected raw-mode block in corruption test')
                raw[HSZ+1]^=0x80
                with self.assertRaises(ValueError):decode_sparse_cascade(bytes(raw))

    def test_sedenion_zero_divisor(self):
        a=np.zeros(16,dtype=np.int8);a[1]=a[10]=1
        b=np.zeros(16,dtype=np.int8);b[4]=1;b[15]=-1
        self.assertTrue(np.all(cd_product_dense(a,b)==0))

if __name__=='__main__':
    unittest.main(verbosity=2)