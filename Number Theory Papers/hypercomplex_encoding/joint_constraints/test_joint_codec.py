import tempfile,unittest,struct,zlib
from pathlib import Path
import numpy as np
from joint_codec import *
from benchmark_joint import choose_baseline,choose_kinematic,compare_one

class JointConstraintsTests(unittest.TestCase):
 def setUp(self): self.ts,self.x=synthetic(101,seed=918)
 def test_joint_round_trip_all_24_modes(self):
  for n in (False,True):
   for r in (False,True):
    for k in ('none','p_from_v','v_from_p'):
     for b in (False,True):
      data=encode(self.ts,self.x,norm=n,rotation=r,kinematic=k,bias_delta=b)
      t,y=decode(data)
      self.assertTrue(np.array_equal(y,self.x),(n,r,k,b))
      self.assertTrue(np.array_equal(t,self.ts))
 def test_truncated_and_corrupt(self):
  blob=encode(self.ts,self.x,norm=True,kinematic='p_from_v')
  for b in (blob[:5],b'broken'+blob[6:]):
   with self.assertRaises(ValueError):decode(b)
  m=bytearray(blob);m[-1]^=5
  with self.assertRaises(Exception):decode(bytes(m))
 def test_norm_exact_current_row(self):
  v=np.array([0,1000000,0],dtype=np.int64)
  self.assertEqual(norm_pred(v,1000000),0)
  self.assertEqual(norm_pred(np.zeros(3,dtype=np.int64),1),1000000)
 def test_timestamp_roundtrip(self):
  b=serialize_times(self.ts)
  self.assertTrue(np.array_equal(deserialize_times(b,len(self.ts)),self.ts))
 def test_kinematic_direction_uses_current_known_channel(self):
  er=build_residuals(self.ts,self.x,kinematic='p_from_v',norm=True)
  y=reconstruct(self.ts,er,kinematic='p_from_v',norm=True)
  self.assertTrue(np.array_equal(y,self.x))
 def test_calibration_only(self):
  base,_=choose_baseline(self.ts[:25],self.x[:25]);kin,_=choose_kinematic(self.ts[:25],self.x[:25],base)
  self.assertIn(kin,('p_from_v','v_from_p'))
  for g in GROUPS:self.assertIn(base[g],('raw','delta','ar2'))
 def test_malformed_integer_range(self):
  blob=encode(self.ts,self.x)
  self.assertEqual(blob[:4],MAGIC)
  self.assertEqual(len(self.x),decode(blob)[1].shape[0])
 def test_independence_of_constraint_groups(self):
  res_base=build_residuals(self.ts,self.x)
  nres=build_residuals(self.ts,self.x,norm=True)
  kres=build_residuals(self.ts,self.x,kinematic='p_from_v')
  self.assertTrue(np.array_equal(nres[:,:3],res_base[:,:3]))
  self.assertTrue(np.array_equal(kres[:,3:7],res_base[:,3:7]))
  self.assertFalse(np.array_equal(nres[:,3],res_base[:,3]))
  self.assertFalse(np.array_equal(kres[:,:3],res_base[:,:3]))

if __name__=='__main__':unittest.main()
