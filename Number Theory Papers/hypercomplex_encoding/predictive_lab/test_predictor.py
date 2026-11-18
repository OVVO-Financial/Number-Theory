"""Verify losslessness, arithmetic origin, causality, and model metadata."""
import os
import random
import unittest
import numpy as np
from hypercomplex_highdim_test import cd_product_dense
from predictive_lab import predictor as p

class PredictiveCodecTests(unittest.TestCase):
    def test_all_models_random_lengths(self):
        rng=random.Random(90210)
        for D in (2,4,8,16,32):
            for n in (0,1,D-1,D,D+3,5*D+7,32*D+1):
                payload=rng.randbytes(n)
                for mode in p.MODELS:
                    frame,meta=p.encode(payload,mode,D,backend='zstd3')
                    self.assertEqual(p.decode(frame),payload,(D,n,mode,meta))

    def test_brotli_and_adaptive(self):
        rng=random.Random(4998)
        data=rng.randbytes(2048)+b'abc123xyz'*200
        for backend in ('brotli11','zstd9'):
            a,meta=p.encode_adaptive(data,16,backend)
            self.assertEqual(p.decode(a),data)
            self.assertIn(meta['model'],p.MODELS)

    def test_cayley_dickson_product_is_actual_bilinear_product(self):
        for D in (2,4,8,16):
            rng=np.random.default_rng(D)
            for _ in range(15):
                x=rng.integers(-40,41,size=D,dtype=np.int64)
                y=rng.integers(-40,41,size=D,dtype=np.int64)
                raw=cd_product_dense(x,y)
                scaled=p.pair_feature(x,y,1)
                denom=int(np.sqrt(D))*128
                truth=np.where(raw>=0,(raw+denom//2)//denom,-((-raw+denom//2)//denom))
                self.assertTrue(np.array_equal(scaled,truth))

    def test_residual_is_all_pairs_bijection(self):
        a=np.tile(np.arange(256,dtype=np.uint8),256)
        b=np.repeat(np.arange(256,dtype=np.uint8),256)
        r=p.encode_residual(a,b)
        self.assertTrue(np.array_equal(a,p.decode_residual(r,b)))

    def test_fitted_coefficients_cannot_see_the_future(self):
        rng=random.Random(47)
        original=rng.randbytes(16*100)
        split=p.default_prefix_frames(len(original),16)*16
        tampered=original[:split]+rng.randbytes(len(original)-split)
        for mode in p.TRAINED:
            a=p.fit_coefficients(original,16,split//16,mode)
            b=p.fit_coefficients(tampered,16,split//16,mode)
            self.assertTrue(np.array_equal(a,b),(mode,a,b))
        self.assertEqual(p.select_prefix_model(original,16),p.select_prefix_model(tampered,16))

    def test_frame_tampering_is_detected(self):
        data=os.urandom(1024)
        frame,_=p.encode(data,'ar2_cd',16,'zstd9')
        frame=bytearray(frame);frame[p.HEADER.size-8]^=8  # recorded CRC32 byte
        with self.assertRaisesRegex(ValueError,'CRC mismatch'):
            p.decode(bytes(frame))

    def test_predictor_is_causal(self):
        previous=np.arange(16,dtype=np.uint8)
        older=np.arange(16,dtype=np.uint8)[::-1]
        weights=np.array([0,300,100,90],dtype=np.int32)
        x=p.predict(previous,older,'ar2_cd',weights)
        y=p.predict(previous.copy(),older.copy(),'ar2_cd',weights)
        self.assertTrue(np.array_equal(x,y))

if __name__=='__main__':unittest.main()
