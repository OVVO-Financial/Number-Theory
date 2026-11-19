import unittest
import numpy as np
from algebraic_codec import (MODES,encode,decode,as_fixed,residuals,recover,mul,conjugate)
from benchmark_real import selected_signals,per_signal_scores

class LosslessTests(unittest.TestCase):
    def test_all_modes_mixed(self):
        rng=np.random.default_rng(9843)
        base=rng.normal(0,0.3,size=(48,4))
        for mode in MODES:
            x=as_fixed(base)
            for packing in ('i32_row','i32_col','varint_col'):
                for compressor in ('none','zlib','brotli'):
                    if mode=='pointwise': # may overflow extreme ratio, fallback has clipping
                        pass
                    blob=encode(x,mode,packing,compressor)
                    y=decode(blob)
                    self.assertTrue(np.array_equal(x,y),(mode,packing,compressor))

    def test_norm_exact_recovery_perturbed(self):
        t=np.arange(83)*0.05
        s=0.5*np.cos(t)
        xyz=np.stack([s,0.3*np.sin(t),0.1*np.cos(t/2)],axis=1)
        w=np.sqrt(1-np.sum(xyz*xyz,axis=1))
        x=as_fixed(np.column_stack([w,xyz]))
        x[12,0]+=7
        x[38,0]-=13
        for mode in ('norm_delta','norm_ar2'):
            result=recover(residuals(x,mode),mode)
            self.assertTrue(np.array_equal(x,result))

    def test_mul_hamilton(self):
        q=np.array([1,2,3,4],dtype=float)
        r=np.array([5,6,7,8],dtype=float)
        self.assertTrue(np.array_equal(mul(q,r),np.array([-60,12,30,24],dtype=float)))
        self.assertTrue(np.allclose(mul(q,conjugate(q)),[30,0,0,0]))

    def test_corruption(self):
        a=as_fixed(np.tile([1.,0.,0.,0.],(18,1)))
        blob=bytearray(encode(a,'norm_delta','varint_col','none'))
        blob[-1]^=0x44
        with self.assertRaises(Exception):decode(bytes(blob))

    def test_source_not_imaginary_advantage_by_construction(self):
        # Changing only w destroys norm coupling, while multiplication of a
        # real sample has correct imaginary signs.
        t=np.arange(32)*0.01
        xyz=np.stack([.1*np.cos(t),.1*np.sin(t),t*.003],axis=1)
        w=np.sqrt(1-np.sum(xyz*xyz,axis=1))
        exact=as_fixed(np.c_[w,xyz])
        self.assertLess(np.median(abs(residuals(exact,'norm_ar2')[:,0])),3)

if __name__=='__main__':unittest.main()
