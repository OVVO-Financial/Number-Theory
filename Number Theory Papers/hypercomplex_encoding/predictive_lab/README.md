# Causal hypercomplex multiplication prediction lab

This extension tests a **different encoding hypothesis** from the earlier signed spectral cascade: whether the *actual Cayley-Dickson product of two earlier data blocks* improves the prediction of the next block, leaving a more compressible residual. Every predictor is causal and every candidate has a working lossless decoder.

## Mathematical experiment

Divide a file into contiguous blocks of `D` bytes, D a power of two. Center byte coordinates as `u_t = x_t - 128`. For the trained model, use

```
feature_t = round( CD_product(u_(t-1), u_(t-2)) / (128 * floor(sqrt(D))) )
pred_t = clip(128 + round_Q(theta0 + theta1*u_(t-1) + theta2*u_(t-2)
                                         + theta3*feature_t), 0, 255)
```

`CD_product` is the **actual Cayley-Dickson multiplication**, imported from `hypercomplex_highdim_test.py`; the signs depend on the basis multiplication rule `e_i e_j = s(i,j)e_(i XOR j)`. The two inputs to multiplication are always *prior* blocks. It is nonlinear in these data coordinates, unlike the previous invertible spectral cascade.

`theta` is fitted by ridge regression on an **initial prefix only**, quantized to signed fixed-point Q10 coefficients, and serialized in the file header. The prefix contains at most 25% of the complete blocks and is capped at 256 blocks; it is compressed and charged to the total frame size. A held-out suffix is predicted sequentially from the preceding *decoded* blocks. Residual bytes use an invertible modulo-256 zigzag mapping followed by the selected compressor. Format version, dimension, predictor ID, coefficients, both compressed streams, original length and CRC32 are all stored in the finished file. There is **no external training state or free side information**.

## Controls

- `ar2`: global scalar intercept and two previous coordinate values, no bilinear feature.
- `ar2_cd`: same plus signed Cayley-Dickson bilinear feature (one extra coefficient).
- `ar2_xorconv`: matched unsigned XOR-index convolution, one extra coefficient.
- `ar2_pointwise`: matched coordinatewise bilinear product, one extra coefficient.
- `ar2_circular`: matched cyclic convolution, one extra coefficient.
- `ar2_shifted`: matched single shifted bilinear product, one extra coefficient.
- `lag1`, `delta2`, `xor_lag`: fixed, low-cost causal predictors.
- `auto`: predictor chosen **only from within-prefix validation bytes** (first 75% prefix train; last 25% validation); then refitted on the full prefix. No test-suffix information enters model selection.

Direct zstd 3/9, Brotli 11 and LZMA 6 are independently measured whole-file baselines. There are no oracle choices in the primary `auto` result. Additional fixed-mode and dimension-sweep results are diagnostics and should never be reported as a deployable min-over-test-set codec.

## Run

From `Number Theory Papers/hypercomplex_encoding`:

```sh
python -m pip install -r requirements.txt
python -m unittest -v predictive_lab.test_predictor
python -m predictive_lab.benchmark --smoke --out /tmp/hc-predictive-smoke --repeats 1
python -m predictive_lab.benchmark --backend brotli11 --out /tmp/hc-predictive-brotli --repeats 2
python -m predictive_lab.sensitivity --out /tmp/hc-predictive-dimensions
```

Both high-level commands verify byte-for-byte reconstruction for *every model*. `benchmark.py` records source SHA-256, encoded size, encode/decode time, MB/s, raw suffix entropy, residual entropy, raw suffix zstd9 size, and residual zstd9 size. The full PDF and complete module sources are read directly from the repository. The benchmark's `--smoke` mode caps each input at 4,096 bytes and is **not** a full-file benchmark.

## Interpretation and reference results

See [RESULTS_2026-11-18.md](RESULTS_2026-11-18.md) and frozen `results/` tables. On six natural/derived test files, the new CD predictive model did **not** beat direct Brotli 11, and on these inputs often produced less compressible residuals than simple linear predictors. Increasing block dimensionality from 4 through 64 did not change that conclusion.

On a **deliberately CD-generated** nonlinear control at D=16, the Cayley-Dickson model *did* find a causal saving relative to linear AR2, unsigned XOR convolution, and direct Brotli. This control is generated from the exact algebra with coefficients selected for an illustrative nondegenerate stable process; it is **not independent natural data**, and it should not be conflated with held-out evidence. On an independent linear-AR control, the simpler linear model wins by a few bytes.

This only tests a single global-scalar-coefficient quadratic feature and the stated normalizations; it does **not** establish that all quaternion, Clifford, or higher-dimensional predictive models fail. More promising target domains would have actual quaternion orientation states or known algebraic recurrence, with externally sourced held-out data and strong domain-specific predictors as baselines. No GPU measurements were performed.
