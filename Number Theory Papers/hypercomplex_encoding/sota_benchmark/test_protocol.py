"""Split-integrity, frozen-config and end-to-end benchmark smoke tests."""
from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path
from sota_benchmark.benchmark import calibrate, evaluate, get_manifest, DEFAULT_MANIFEST, digest


class TestProtocol(unittest.TestCase):
    def test_manifest_disjoint(self):
        m=get_manifest(DEFAULT_MANIFEST)
        train={x['path'] for x in m['calibration']}
        hold={x['path'] for x in m['heldout']}
        self.assertFalse(train & hold)
        self.assertGreaterEqual(len(m['heldout']),3)

    def test_calibration_does_not_read_heldout_and_freezes(self):
        # Deliberately missing heldout file: calibration must still succeed.
        example=json.loads(DEFAULT_MANIFEST.read_text())
        example['heldout']=[{'id':'unavailable_future_file','path':'sota_benchmark/does-not-exist.test','family':'real_text'}]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            manifest=root/'manifest.json';manifest.write_text(json.dumps(example))
            config=root/'calibration.json'
            info=calibrate(manifest,config,smoke=True)
            self.assertTrue(info['heldout_files_not_accessed'])
            self.assertEqual(len(info['training_files']),2)
            self.assertEqual(info['max_bytes'],3072)
            self.assertEqual(digest(__import__('base64').b64decode(info['dictionary_base64'])),info['dictionary_sha256'])
            # Evaluation needs the heldout file and cannot silently substitute the train set.
            with self.assertRaises(FileNotFoundError):
                evaluate(manifest,config,root/'results',smoke=True)

    def test_frozen_config_rejects_calibration_mutation(self):
        m=json.loads(DEFAULT_MANIFEST.read_text())
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            mf=root/'manifest.json';mf.write_text(json.dumps(m))
            config=root/'calibration.json'
            info=calibrate(mf,config,smoke=True)
            info['training_files'][0]['sha256']='0'*64
            config.write_text(json.dumps(info))
            with self.assertRaises(ValueError):
                evaluate(mf,config,root/'results',smoke=True,allow_missing=True)


if __name__=='__main__':unittest.main()
