# Matched lossless codec results on EuRoC estimated quaternion trajectories

**Date:** November 18, 2026

**CI verified:** https://github.com/OVVO-Financial/Number-Theory/actions/runs/44432265843

All four real trajectory tests, exact round trips, corruption tests, and matching compressor comparisons passed on GitHub Actions. The result tables are in `EUROC_FULL_CODEC_RESULTS.csv`; the source file names and SHA256 are in `EUROC_FULL_CODEC_MANIFEST.json`.

## Bit accounting

The source represents pose orientation quaternions in a public visual-inertial **trajectory estimate**, not raw Vicon or IMU observations. The quaternion coefficients `(w,x,y,z)` were parsed as floats and rounded to signed 32-bit integers at 1e-6 precision. Each tested codec is lossless on these fixed-point integers, including all corrections, framing, CRC32, and compressor bytes. The original CSV or float64 byte patterns are **not** reconstructed. The first chronological quarter of each trajectory was used only to freeze packing and entropy backend; performance is reported on the remaining three quarters.

## Complete measured bytes: same model-specific codec frame and calibrated compressor

| Source | Held-out samples | Direct raw | AR2 | Hamilton-product | Norm + AR2 | Saving vs AR2 |
|---|---:|---:|---:|---:|---:|---:|
| MH_01_mono | 2,745 | 23,663 | 16,841 | 16,838 | **13,315** | **20.94%** |
| MH_03_mono | 1,924 | 16,827 | 12,187 | 12,179 | **9,724** | **20.21%** |
| MH_05_mono | 1,684 | 15,703 | 10,331 | 10,329 | **8,166** | **20.96%** |
| MH_05_stereo | 1,705 | 14,657 | 10,506 | 10,504 | **8,311** | **20.89%** |
| **Total** | **8,058** | **70,850** | **49,865** | **49,850** | **39,516** | **20.75%** |

Every recorded quaternion result selected **Brotli-11 with signed-varint column packing** on its calibration prefix. The same possible compressor/packing choices were available to every model, so this is a matched backend comparison. The measured improvement is **not** a claim to beat all general-purpose SOTA codecs or standardized quaternion representations.

## Critical negative controls

| Signal | AR2 combined bytes | Norm + AR2 combined bytes | Norm gain |
|---|---:|---:|---:|
| Shuffled scalar w | 58,515 | 57,465 | +1.79% |
| Position XYZ + quaternion scalar w (not a unit quaternion) | 55,099 | 65,002 | -17.97% |
| Time-shuffled complete quaternions | 84,369 | 66,733 | +20.90% |

Interpretation: the norm advantage almost disappears when the scalar component is disconnected from the vector components, and becomes a penalty on unrelated four-channel data. It persists after time shuffling because it is **a same-row geometric constraint**, not a forecasting effect. The multiplication-only predictor was essentially tied with AR2.

## Limitations and next comparisons

- The source is already normalized *by design*: quaternion representation lives approximately on the 3-sphere. This is a familiar constraint, **not a discovery that imaginary numbers carry more physical information per bit**.
- Rounding to 1e-6 is quantization. Exact reproduction of the original source text, float64 bits, and other modalities was **not** tested in this run.
- A minimal three-coordinate quaternion codec, established orientation formats (e.g. smallest-three), axis-angle, and rotation-vector parametrizations must be benchmarked before declaring competitive gains.
- Higher Cayley-Dickson dimensions did not contribute to this result: the useful structure is the standard 4D quaternion norm.
- Negative controls cannot establish mathematical novelty; they only show that the observed coding advantage depends on the source relation.

## Reproduction
From `constraint_discovery/` with Python3.12 and `numpy brotli zstandard` installed:
```sh
python -m unittest -v test_algebraic_codec.py
python local_controls.py
python benchmark_euroc.py --download --out results/euroc
```
The script records original upstream SHA256s and performs exact decoding before declaring success. Original third-party dataset recordings are not bundled.

## Source
https://github.com/uzh-rpg/rpg_trajectory_evaluation
The project cites Zhang and Scaramuzza, IROS 2018.