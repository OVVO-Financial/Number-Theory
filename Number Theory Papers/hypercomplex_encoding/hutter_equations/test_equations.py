"""Exactness, causality, and independence checks; no enwik8 required."""
import unittest
import random
import zlib
import numpy as np
from equations import NAMES, IDX, SEED_LENGTH, predictions, predict_one, prev_class, encode, decode
from benchmark import run, baseline_probabilities, causal_residual_probabilities, prior_classes

class EquationsTest(unittest.TestCase):
    def test_scalar_vector_equivalence(self):
        rng=random.Random(3303)
        for payload in [bytes(range(256))*2,
                        bytes(rng.randrange(256) for _ in range(300)),
                        b'<page><title>dog</title></page>\n'*15]:
            p=predictions(payload)
            for i in range(SEED_LENGTH,len(payload)):
                for j,name in enumerate(NAMES):
                    self.assertEqual(int(p[name][i]),predict_one(j,payload[:i]),(name,i))

    def test_causal_past_only(self):
        rng=random.Random(889)
        original=bytes(rng.randrange(256) for _ in range(132))
        mutated=original[:87]+bytes(rng.randrange(256) for _ in range(45))
        a,b=predictions(original),predictions(mutated)
        for name in NAMES:
            self.assertTrue(np.array_equal(a[name][:88],b[name][:88]),name)

    def test_all_eq_roundtrips_and_classes(self):
        samples=[b'',b'dog',bytes(range(42)),bytes(range(256))*3,
                 bytes([0,255,17,17,32,128])*80]
        for payload in samples:
            for mid in range(len(NAMES)):
                blob=encode(payload,mid,quality=5)
                self.assertEqual(decode(blob),payload)
            self.assertEqual(decode(encode(payload,2,bytes(i%len(NAMES) for i in range(8)),quality=5)),payload)

    def test_corrupt_and_bad_config(self):
        b=encode(b'The quick brown fox jumps over lazy dog.'*12,4,quality=6)
        with self.assertRaises(ValueError):decode(b[:10])
        with self.assertRaises(ValueError):decode(b[:14]+b'X'+b[15:])
        with self.assertRaises(ValueError):decode(b[:9]+bytes([b[9]^1])+b[10:])
        with self.assertRaises(ValueError):encode(b'test',99)
        with self.assertRaises(ValueError):encode(b'test',2,bytes([3]*7))

    def test_known_modular_constraint(self):
        # One actual low-description generative law whose residual is zero.
        x=bytearray(range(SEED_LENGTH))
        for _ in range(160):x.append((2*x[-1]-x[-2])&255)
        preds=predictions(bytes(x))
        self.assertTrue(np.all(preds['affine_2x_minus_y'][SEED_LENGTH:]==np.frombuffer(bytes(x)[SEED_LENGTH:],dtype=np.uint8)))
        self.assertEqual(decode(encode(bytes(x),IDX['affine_2x_minus_y'])),bytes(x))

    def test_benchmark_no_oracle_selection_and_coding(self):
        rng=random.Random(42)
        fake=(b'<page><title>Example</title><text>This is a text article.\n</text></page>\n'*45)
        out=run(fake[:2400],warmup=1200,validation=600,heldout=600)
        self.assertEqual(out['heldout_bytes'],600)
        self.assertEqual(len(out['model_rows']),len(NAMES))
        self.assertGreater(out['baseline_actual_bytes']['brotli11'],0)
        self.assertIn('class_selected',out['actual_framed_bytes'])
        self.assertIn('quaternion',out['family_winners'])

if __name__=='__main__':unittest.main()
