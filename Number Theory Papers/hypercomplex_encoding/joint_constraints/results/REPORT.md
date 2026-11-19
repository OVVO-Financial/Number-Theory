# Joint algebraic-constraint discovery on EuRoC estimated motion states

All reported sizes are standalone **quantized-integer** lossless frames: 16 state channels at scale 1e6, timestamp varints, model flags, CRC32, and Brotli9 residual payload. Source CSV decimal text/float64 byte patterns are **not** preserved.

Each sequence is split chronologically: first 25% trains temporal modes and the preferred p/v kinematic direction; last 75% is evaluated without tuning.

| Sequence | Samples | Temporal baseline | Quaternion norm N | Kinematic K | N+K | All N+K+R+B | Combined vs temporal |
|---|---:|---:|---:|---:|---:|---:|---:|
| MH_01_easy | 2,730 | 65,444 | 62,059 | 62,273 | 58,805 | 58,782 | +10.18% |
| MH_03_medium | 1,974 | 47,855 | 45,907 | 45,183 | 43,151 | 43,147 | +9.84% |
| MH_05_difficult | 1,667 | 39,717 | 38,924 | 37,841 | 36,193 | 36,188 | +8.89% |

## Independent constraint family and incremental contribution

* **N**: quaternion unit norm, reconstruct scalar coordinate from current vector part, then transmit exact correction.
* **K**: translation–velocity kinematics, using measured contemporaneous velocity to predict position, or measured position to predict velocity; choose direction exclusively on calibration prefix.
* **R**: temporal quaternion Hamilton-product recurrence, a separate motion assumption from the same-row norm.
* **B**: near-static gyro/accelerometer biases, predicting bias from its prior decoded value. This is a temporal prior, not an exact physical equality or a hypercomplex-specific result.

Rotation-matrix orthogonality and unit-quaternion conjugation are intentionally **not** counted as independent equations: they follow from the same orientation norm.

## Per-sequence incremental net compressed bytes

| Dataset | N gain | K gain after N | R gain after N+K | B gain after N+K+R |
|---|---:|---:|---:|---:|
| MH_01_easy | +3,385 | +3,254 | +23 | +0 |
| MH_03_medium | +1,948 | +2,756 | +4 | +0 |
| MH_05_difficult | +793 | +2,731 | +5 | +0 |

A positive gain means fewer bytes; a negative gain means the proposed constraint harms compression. **Independence of physical assumptions does not guarantee additive information savings**, especially when baseline temporal prediction already exploits the same regularity.

## Negative controls

| Control | Temporal | N | K | B |
|---|---:|---:|---:|---:|
| q_w_disconnected | 68386 | 67971 | 65115 | 68386 |
| velocity_disconnected | 70283 | 66851 | 75267 | 70283 |
| bias_disconnected | 85708 | 80762 | 81011 | 85708 |

## Critical scope

- This is EuRoC *estimated ground-truth states*, not raw inertial observation channels. Estimates may be normalized and smoothed by the estimator.
- Kinematics are approximate in sampled/estimated data; corrections are always transmitted, so no quantized coefficient is literally omitted unless its correction is zero.
- The 16D state vector is ordinary position/quaternion/velocity/bias data; higher Cayley–Dickson dimension contributes no new saving.
- The comparison is against a strong temporal predictor with the same bit framing/backend, not an exhaustive competition with PCL or specialized state-estimate codecs.
- Dataset originals are downloaded by the reproducible script, with immutable SHA256s recorded in source_manifests.json.

## Reproduction

```sh
python -m pip install -r requirements.txt
python -m unittest -v test_joint_codec
python benchmark_joint.py --download --stride 10 --datasets MH_01_easy MH_03_medium MH_05_difficult
```
