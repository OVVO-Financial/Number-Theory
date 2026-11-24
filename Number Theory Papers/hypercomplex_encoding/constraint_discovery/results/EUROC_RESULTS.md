# Real-data arithmetic-constraint experiment: first verified results

Date: 2026-11-18

## Source and target of the exact reconstruction
Four publicly available estimated trajectories from `uzh-rpg/rpg_trajectory_evaluation`, produced from the EuRoC visual/inertial data. These are **trajectory estimates** (not raw IMU samples and not the ground-truth files).
Each sample was parsed as `(w,x,y,z)` and rounded to nearest integer after multiplying by 1,000,000.
The result is **lossless on these quantized integer samples**, not on the original CSV or IEEE float64 bytes.
The last 75% of each trajectory is the evaluation suffix, with no parameters fitted from it.

## Varint coding results

| File | Samples | AR2 | Quaternion multiplication | Norm + AR2 | Norm saving |
|---|---:|---:|---:|---:|---:|
| MH_01_mono | 2,745 | 21,386 | 21,406 | 18,788 | 12.15% |
| MH_03_mono | 1,924 | 15,071 | 15,067 | 13,217 | 12.30% |
| MH_05_mono | 1,684 | 13,070 | 13,086 | 11,488 | 12.10% |
| MH_05_stereo | 1,705 | 13,223 | 13,219 | 11,619 | 12.13% |
| **Combined** | **8,058** | **62,750** | **62,778** | **55,112** | **12.17%** |

These values are signed zigzag-varint byte counts for every coefficient plus the same **13-byte hypothetical framing allowance**. In the connector-based experiment, the arithmetic prediction and inverse were tested exactly, and byte lengths were computed from the varint rule. A physically serialized 18-byte-header codec and matched zlib/Brotli/Zstd comparisons are provided as a separate reproducible follow-up script (`benchmark_euroc.py`) and should not be confused with the table above.

## Mathematically informative negative controls

| File | AR2, shuffled w | Norm + AR2, shuffled w | AR2, nonquaternion 4D | Norm + AR2, nonquaternion 4D |
|---|---:|---:|---:|---:|
| MH_01_mono | 24,146 | 23,802 | 21,718 | 26,704 |
| MH_03_mono | 17,053 | 16,937 | 15,604 | 18,736 |
| MH_05_mono | 14,730 | 14,643 | 13,726 | 16,905 |
| MH_05_stereo | 14,924 | 14,811 | 14,197 | 17,279 |

The matching source quaternion norm was within approximately 4.1-4.3e-7 of unity at median. The estimated scalar-coordinate residual entropy fell from roughly 10.3-10.8 bits under AR2 to roughly 1.4-2.0 bits after norm reconstruction (before the corrector is entropy-coded).

## What this means

1. **A real source-specific constraint is present:** unit-quaternion orientation lies close to the three-dimensional hypersphere in R4. Estimating the fourth coordinate from the other three makes its residual cheap.
2. **This is not a new imaginary-number capacity theorem:** the same saving follows from standard unit-quaternion geometry without higher Cayley-Dickson dimensions.
3. **Quaternion multiplication was not superior to AR2:** the tested constant-relative-motion recurrence was almost exactly tied.
4. **No SOTA claim:** this result compares lossless signed-varint encodings on quantized samples. It does not demonstrate a win over Brotli, Zstandard, or a specialized 3-parameter rotation codec.

## Provenance (Git blob SHA is immutable)

- `MH_01_mono` SHA `ae9c92a274b7369b42b53a9626712f8d6630a3b3`: https://github.com/uzh-rpg/rpg_trajectory_evaluation/blob/master/results/euroc_mono_stereo/laptop/vio_mono/laptop_vio_mono_MH_01/stamped_traj_estimate.txt
- `MH_03_mono` SHA `19e0b5d45ba4a91a762ab562fb2469b25b46941b`: https://github.com/uzh-rpg/rpg_trajectory_evaluation/blob/master/results/euroc_mono_stereo/laptop/vio_mono/laptop_vio_mono_MH_03/stamped_traj_estimate.txt
- `MH_05_mono` SHA `350f88bb53868efb57623e2bfd7af31782ff83c8`: https://github.com/uzh-rpg/rpg_trajectory_evaluation/blob/master/results/euroc_mono_stereo/laptop/vio_mono/laptop_vio_mono_MH_05/stamped_traj_estimate.txt
- `MH_05_stereo` SHA `cb8854883ff70227c2c4afb07ac68e952bc176bc`: https://github.com/uzh-rpg/rpg_trajectory_evaluation/blob/master/results/euroc_mono_stereo/laptop/vio_stereo/laptop_vio_stereo_MH_05/stamped_traj_estimate.txt

## Reproduction
`python -m pip install numpy brotli zstandard`
`python -m unittest -v test_algebraic_codec.py`
`python benchmark_euroc.py --download --out results/euroc` (with network access)
`python local_controls.py` (constructed positive/negative controls, clearly separate).
No externally downloaded source data is redistributed in the package.
