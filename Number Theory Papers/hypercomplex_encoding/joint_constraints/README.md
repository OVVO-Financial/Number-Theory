# Multiple independent algebraic constraints on real motion data

This experiment tests whether **two physically independent source constraints**
provide cumulative, exact fixed-point compression on EuRoC estimated drone states:

1. Unit quaternion norm predicts the current scalar component from XYZ (N).
2. Translational kinematics predicts the current position from measured velocity or the reverse (K).

It separately investigates quaternion angular-motion recurrence (R) and
inertial bias near-constancy (B). R is an additional temporal hypothesis on the
same orientation coefficients; B is an ordinary random-walk prior, not a new
hypercomplex identity. The last two are deliberately not claimed as independent
same-sample algebraic equations.

Predictions are exact **only after transmitting integer correction residuals**.
Input is quantized to signed integer micro-units (1e-6). The method is not
lossless on original CSV text/float64 representations.

## Held-out, immutable evaluation

Each source's first chronological quarter calibrates the ordinary predictor for
four disjoint coefficient groups and the p/v kinematic encoding direction.
The last 75% is evaluated with those modes frozen. The timestamp vector itself
is included in every standalone encoded file, as are flags, checksum,
Brotli9 residual bytes, and residual headers.

All 17 configurations use the **same lossless container, residual packing, and
entropy coder**, with exact independent decoder checks. Shuffled/wrong-channel
controls establish whether a measured saving truly depends on the claimed
relation rather than on implicit representation choices.

## Data

OpenVINS corrected EuRoC groundtruth-state CSV mirrors:
https://github.com/rpng/open_vins/tree/master/ov_data/euroc_mav

Primary authority: https://projects.asl.ethz.ch/datasets/euroc-mav/

Sample datasets: MH_01_easy, MH_03_medium, MH_05_difficult.
The CSV columns must be the standard EuRoC 17-field order: nanosecond time,
position XYZ, quaternion wxyz, velocity XYZ, gyro bias XYZ, accel bias XYZ.
Raw dataset files are downloaded by GitHub Actions and NOT committed.
A source SHA256 for every observed file is recorded in source_manifests.json.

## Reproduce

```sh
python -m pip install -r requirements.txt
python -m unittest -v test_joint_codec
python benchmark_joint.py --synthetic --out /tmp/joint-synthetic
python benchmark_joint.py --download --stride 10 --datasets MH_01_easy MH_03_medium MH_05_difficult
```

Synthetic results verify implementation; only GitHub Actions' downloaded real
recordings count as real-world evidence. A positive result for N+K on held-out
sources would demonstrate cumulative compression, not mathematical novelty or
SOTA performance. If K is negative or zero, the result will report that gap.
