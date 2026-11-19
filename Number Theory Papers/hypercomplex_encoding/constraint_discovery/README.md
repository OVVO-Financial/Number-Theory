# Arithmetic constraints of imaginary coefficients on real telemetry

This experiment **does not** ask whether imaginary axes can hold bits. It asks whether the *source distribution* obeys algebraic relations making some values easier to represent.

We evaluate two quaternion-specific constraints on **held-out real Vicon orientation data** (NanoBench, https://github.com/syediu/nanobench):

1. **Unit norm:** with fixed-point coordinates `(w,x,y,z)`, predict `w` from `sqrt(max(0,1-x²-y²-z²))` and its previous sign, then encode the *exact correction*. `x,y,z` are predicted using delta or AR2. For real measurements the norm is only approximate and the correction is essential.
2. **Hamilton multiplication:** predict `q[t] = (q[t-1] * conjugate(q[t-2]) / ||q[t-2]||²) * q[t-1]`. This uses only previously decoded samples, and stores exact residuals. Norm or fixed relative rotation is **not** assumed to be exact.

Baselines: raw 4-channel int32, channel delta, AR2, and pointwise multiplicative extrapolation. All methods share the same 18-byte file framing, CRC, deterministic source precision, and the same calibrated packing/compression selection (zlib9 / Brotli11 / Zstandard9 when installed, row/column int32 or column zigzag-varints). Coding/decoding is checked bit-for-bit at the **quantized int32** level. This **does not** losslessly preserve original CSV text or IEEE float64 values.

Controls: independently shuffled scalar component `w`, time-shuffled quaternion rows, and non-quaternion 4D position/gyroscope measurements. The source CSV is never checked in, but exact filenames, download URLs, hashes, and row counts are recorded in `manifest.json`.

## Run locally (with network)

```bash
python -m pip install numpy brotli zstandard
python -m unittest -v test_algebraic_codec.py
python local_controls.py  # constructed positive/negative controls; NOT real data
python benchmark_real.py --download --max-rows 1600 --max-flights 4 --out results/real
```

To use existing CSV files instead, `python benchmark_real.py --csv /path/to/flight.csv --out results/real`.

**Interpretation:** Any improvement over AR2 on held-out quaternion sequences is evidence of useful source-specific arithmetic structure, not proof the Cayley-Dickson/complex representation has more intrinsic storage capacity. Better quaternion parametrizations and specialized compressors remain important controls in a mature competitive benchmark.

The signal source is independent of this repository. NanoBench is BSD-3-Clause licensed; this PR downloads only public files for CI and does not redistribute them.
