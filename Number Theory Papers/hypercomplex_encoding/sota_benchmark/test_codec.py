"""Independent codec and fairness regression tests (safe for CI)."""
from __future__ import annotations
import os
import random
import struct
import unittest
import numpy as np
from sota_benchmark.codec import (
    Zstd, HybridConfig, VARIANT, hybrid_encode, hybrid_decode,
    fwht, inverse_fwht, zstd_dictionary_encode, zstd_dictionary_decode,
    _spectral_encode, _spectral_decode, compress_named, decompress_named,
    basis_codebook_encode, basis_codebook_decode
)
from hypercomplex_cascade_test import cd_cascade_inverse_checked


class TestCodec(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z=Zstd()
        # Warm the JIT once before benchmarking; tests verify truth, not speed.
        fwht(np.array([1,-1,1,-1],dtype=np.int8))

    def test_highdim_exact_roundtrip_and_overhead(self):
        rng=random.Random(0xC0DE)
        for n in (0,1,2,3,11,31,32,100,201):
            for data in (bytes(n),b'\xff'*n,bytes(rng.randrange(256) for _ in range(n))):
                for v in VARIANT:
                    cfg=HybridConfig(block_bytes=32,variant=v)
                    stats={}
                    frame=hybrid_encode(data,cfg,self.z,stats)
                    self.assertEqual(hybrid_decode(frame,self.z),data)
                    self.assertGreaterEqual(len(frame),23)
                    self.assertEqual(stats['blocks'],(len(data)*8+254)//255)

    def test_hybrid_truncated_and_corruption(self):
        data=b'abcd'*32
        frame=hybrid_encode(data,HybridConfig(32,'hypercomplex'),self.z)
        with self.assertRaises(ValueError):hybrid_decode(frame[:-1],self.z)
        bad=bytearray(frame);bad[20]^=0xFF  # CRC area of header
        with self.assertRaises(ValueError):hybrid_decode(bytes(bad),self.z)
        with self.assertRaises(ValueError):hybrid_decode(frame+b'\0',self.z)

    def test_spectral_variant_distinction(self):
        rng=np.random.default_rng(23)
        for n in (64,256,1024):
            vec=rng.integers(0,2,size=n-1,dtype=np.uint8)
            data=np.packbits(vec[:((n-1)//8)*8]).tobytes()
            for variant in ('hypercomplex','walsh'):
                s,b,nnz=_spectral_encode(np.unpackbits(np.frombuffer(data,np.uint8)),n,variant)
                # A random string is normally dense and therefore spectral may be rejected.
                self.assertGreater(nnz,0)
                coeff=np.full(n,-1,dtype=np.int8);coeff[0]=1
                bits=np.unpackbits(np.frombuffer(data,np.uint8))
                coeff[1:len(bits)+1]=2*bits.astype(np.int8)-1
                spectrum=fwht(coeff) if variant=='walsh' else __import__('hypercomplex_cascade_test').cd_cascade(coeff)
                recovered=inverse_fwht(spectrum) if variant=='walsh' else cd_cascade_inverse_checked(spectrum)
                np.testing.assert_array_equal(recovered,coeff)

    def test_basis_index_special_family(self):
        for n in (256,1024):
            for idx in (0,1,53,n-1):
                target=np.zeros(n,dtype=np.int64);target[idx]=n
                coeff=cd_cascade_inverse_checked(target)
                if coeff[0]<0:coeff=-coeff
                data=np.packbits(((coeff[1:]+1)//2).astype(np.uint8)).tobytes()
                cfg=HybridConfig(n//8,'basis_index')
                d=hybrid_encode(data,cfg,self.z)
                self.assertEqual(hybrid_decode(d,self.z),data)
                self.assertLess(len(d),len(data)+25)

    def test_independent_basis_index_codebook_format(self):
        from sota_benchmark.benchmark import synthetic_basis_signal
        coded=synthetic_basis_signal(D=1024,blocks=8)
        frame=basis_codebook_encode(coded,1024)
        self.assertIsNotNone(frame)
        self.assertEqual(basis_codebook_decode(frame),coded)
        self.assertLess(len(frame),len(coded))
        with self.assertRaises(ValueError):basis_codebook_decode(frame[:-1])
        # The codebook applies only to the stated basis family.
        self.assertIsNone(basis_codebook_encode(os.urandom(1023),1024))

    def test_zstandard_levels(self):
        for level in (1,3,9,19):
            for data in (b'',b'hello world'*40,os.urandom(128)):
                c=self.z.compress(data,level)
                self.assertEqual(self.z.decompress(c,len(data)),data)

    def test_calibration_only_dictionary_embedded(self):
        samples=[(b'field:amount:'+bytes([65+(j%20)])+b' primary-codec')*100 for j in range(30)]
        dictionary=self.z.train(samples,2048)
        payload=b''.join(samples[:3])+b'field:amount:X'
        frame=zstd_dictionary_encode(payload,self.z,dictionary,9)
        self.assertEqual(zstd_dictionary_decode(frame,self.z),payload)
        self.assertIn(dictionary,frame)
        self.assertGreater(len(frame),len(dictionary)+16)

    def test_conventional_roundtrip(self):
        data=b'abaabb'*90
        for name in ('zstd1','zstd3','zstd9','zstd19','brotli5','brotli9','brotli11','lzma6','lzma9','zlib6','zlib9'):
            encoded=compress_named(data,name,self.z)
            self.assertEqual(decompress_named(encoded,name,self.z,len(data)),data)

    def test_no_hidden_unserialised_dictionary(self):
        data=b'confidential'*12
        cfg=HybridConfig(32,'hypercomplex')
        raw=hybrid_encode(data,cfg,self.z)
        self.assertEqual(hybrid_decode(raw,Zstd()),data)


if __name__=='__main__':unittest.main()

# Allow independent execution via -m unittest discover as well as direct module invocation.
