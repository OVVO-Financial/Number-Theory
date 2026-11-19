# Independent algebraic constraints: real EuRoC flight-state results

**Run:** [GitHub Actions 44435544867](https://github.com/OVVO-Financial/Number-Theory/actions/runs/44435544867) — all steps green, including independent decoder checks.  
**Data:** OpenVINS-corrected EuRoC `MH_01_easy`, `MH_03_medium`, `MH_05_difficult` state-estimate CSV files; input URLs and raw-file SHA256s are recorded in `results/source_manifests.json`.  
**Transform:** Fixed-point integer representation at **1e-6 units** of 16 state coefficients (position XYZ, quaternion wxyz, velocity XYZ, six inertial biases), plus nanosecond timestamps. This is lossless for those integers, **not** for original CSV text / float64 files. 

## Principal result: two distinct physical constraints save bits cumulatively

Each real sequence was downsampled by a fixed `--stride 10`, split chronologically, and the first 25% was used only to choose ordinary predictor modes for each disjoint coefficient group and the kinematic prediction direction. All results use the final 75%; each packet was decoded independently and checked bit-for-bit against the quantized input.

| Held-out flight | Samples | Frozen temporal baseline | Norm only (N) | Kinematics only (K) | N+K | N+K+R+B |
|---|---:|---:|---:|---:|---:|---:|
| MH_01_easy | 2,730 | 65,444 | 62,059 | 62,273 | **58,805** | 58,782 |
| MH_03_medium | 1,974 | 47,855 | 45,907 | 45,183 | **43,151** | 43,147 |
| MH_05_difficult | 1,667 | 39,717 | 38,924 | 37,841 | **36,193** | 36,188 |
| **Total** | **6,371** | **153,016** | **146,890** | **145,297** | **138,149** | **138,117** |

The two-constraint method saves **14,867 bytes (9.72%)** against the frozen temporal baseline. The combined saving is greater than the sum of the single-constraint improvements by 1,022 bytes due to nonlinear effects in the common entropy coder and joint residual statistics; this must not be interpreted as 1,022 bits of newly created information.

## Genuine independence vs. redundant identities

- **N, quaternion unit norm:** predicts the quaternion scalar w from same-row xyz and prior sign, then transmits the exact integer correction. It operates on orientation.
- **K, position/velocity kinematics:** using the actual measured velocity, predicts position by trapezoidal integration; alternatively velocity from position difference. Training selects `p_from_v` for all three flights. It operates on translational state, not on orientation. These are **separate source equations**.
- **R, constant relative quaternion rotation:** Hamilton-product predictor for successive orientations. It adds only **32 total bytes** of gain after N+K (23/4/5 B), not a compelling independent compression result. It is a dynamic assumption about the same orientation channel, not a second algebraic norm identity.
- **B, near-constant estimated IMU biases:** the ordinary temporal baseline already selected delta coding for bias channels. Consequently B adds **0 bytes** in the three real sequences. It is an ordinary time-series assumption, not a new equation.

Quaternion conjugation, q·q̇=0, SO(3) matrix orthogonality and multiplicative unit norm were not counted as independent constraints: they follow from the same normalization manifold.

## Negative controls (MH_01 holdout)

| Control | Temporal baseline | Norm N | Kinematic K | Bias B |
|---|---:|---:|---:|---:|
| Quaternion scalar w independently permuted | 68,386 | 67,971 | 65,115 | 68,386 |
| Position–velocity alignment disrupted | 70,283 | 66,851 | **75,267** | 70,283 |
| Six bias channels independently time-permuted | 85,708 | 80,762 | 81,011 | 85,708 |

The kinematic predictor becomes **worse than ordinary temporal coding when velocities are independently shuffled**. This supports that its gain on the recorded flights relies on a real position–velocity relation. The norm gain is much smaller when w is disconnected from xyz. These are diagnostic controls, not causal proofs.

## What this does and does not prove

1. **Proved computationally on these measured files:** two independent source-specific physical relations can jointly lower complete, lossless *fixed-point* storage with exact correction residuals, beyond tuned single-channel raw/delta/AR2 controls.
2. **Not proved:** that an entire coefficient can be omitted without a correction, that all four candidate constraints add value, that arbitrary Cayley–Dickson multiplication yields a new saving, or that this beats specialized state-vector codecs or SOTA compressors.
3. These EuRoC files represent already-estimated and normalized states. The quaternion norm and position/velocity consistency are known physics used in estimation. This is **applied constraint-aware coding**, not a new fundamental theorem.
4. **A next genuinely independent group** would require synchronized raw gyroscope and accelerometer measurements: test quaternion angular kinematics `q_dot = (1/2) q * omega` and translational acceleration `v_dot = R(q) * (a - b_a) + g` with exact residuals. Those sensor channels were **not available in this experiment** and no gain is claimed for them.

## Reproduce

```sh
python -m pip install -r requirements.txt
python -m unittest -v test_joint_codec
python benchmark_joint.py --download --stride 10 --datasets MH_01_easy MH_03_medium MH_05_difficult --out results/euroc
```

The full 17-ablation matrix, exactness flags, encode/decode timings, and provenance are tracked in `results/results.csv`, `results/negative_controls.json` and `results/source_manifests.json`. The CI artifact for run 44435544867 also retains the synthetic controlled comparison.
