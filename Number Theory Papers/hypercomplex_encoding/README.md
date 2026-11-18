# Full high-dimensional imaginary-coefficient encoding experiment

**Source:** *Real-Pair Representations of Complex Numbers and Quaternions*, November 4, 2026.

This package actually executes Cayley-Dickson imaginary-basis multiplication and extends the paper's `phi` and `psi` mappings by composing successive `(1+e_j)` multiplications through 1,048,576 real dimensions. It also tests associative Clifford vector encoding, exact multi-bit imaginary coefficients, true lossless file codecs, spectral sparsity and compression. **This is not just a tensor Walsh-Hadamard experiment.**

## Quick results

- Exact bit encoding and unit multiplication up to **1,048,575 imaginary axes**.
- Exact 1/2/4/8/16-bit per-axis encodings at the same maximum dimension.
- True signed spectral cascade and inverse, full random vector at 1,048,576 dimensions.
- Dense hypercomplex multiplication through dimension 4,096; sparse hypercomplex multiplication through dimension 1,048,576.
- 16-dimensional zero divisor: `(e1+e10)(e4-e15)=0` with nonzero inputs.
- Synthetic basis-codeword of **131,072 bytes** encoded to **29 bytes**, losslessly, versus **274 bytes** with direct Brotli.
- No improvement on the tested natural document/text/code/bitmap or random data.
- Seven independent regression tests pass, including CRC corruption detection.

See [`hypercomplex_highdim_full_report.md`](hypercomplex_highdim_full_report.md) for detailed methods, limitations, tables, and interpretations.

## Install dependencies

```bash
python -m pip install -r requirements.txt
```

All benchmarks run locally without Internet access after dependencies are available. Numba compilation is excluded from per-operation timing measurements.

## Reproduce the experiment

```bash
sh run_full_suite.sh
```

Or run individual phases:

```bash
python -m unittest -v hypercomplex_suite_tests
python hypercomplex_highdim_test.py
python hypercomplex_cascade_test.py --phase codec
python hypercomplex_cascade_test.py --phase spectrum
python hypercomplex_widths_and_errors.py
python hypercomplex_million_dim_demo.py
python hypercomplex_combine_results.py
```

## Encode and decode your own file

```bash
python hypercomplex_codec_cli.py encode realpair_text.txt out.hcb --mode cascade --dimension 4096 --coefficient bipolar --compress brotli
python hypercomplex_codec_cli.py decode out.hcb reconstructed.txt
cmp realpair_text.txt reconstructed.txt
```

Alternative reversible signed imaginary-unit multiplication mode:

```bash
python hypercomplex_codec_cli.py encode realpair_text.txt out.hcz --mode unit --dimension 4096 --compress zlib
python hypercomplex_codec_cli.py decode out.hcz reconstructed.txt
```

Both modes use actual Cayley-Dickson multiplication, with the cascade mode mixing coordinates. The compressed output includes checksum, dimensional metadata, block tags and original length. **Do not assume arbitrary high-dimensional algebra elements are invertible:** zero divisors exist from dimension 16 onward.

## Input data provenance

The original paper is already tracked one directory above as `../Real_Pair_Representation_of_Complex_Numbers_and_Quaternions.pdf` (not duplicated here). The checked-in extracted text and source-code fixtures reproduce the text/code benchmarks. The optional raster `realpair_renders/page-3.png` was used in the historical bitmap benchmark; if absent, the bitmap dataset is skipped rather than synthesized. Random and basis-aligned data are generated deterministically from the published seeds in the scripts.

## Notes on comparisons

Reference CSVs in `results/` are frozen snapshots from the November 18, 2026 experiment. The benchmark chooses the smallest result across the tested dimensions and binary/bipolar representations per whole file. The configuration-search time is excluded from selected-codec time measurements. Compression baselines are direct zlib9, Brotli9 and LZMA6; no claim is made of beating all possible codecs or existing specialized codebooks. The 29-byte result comes from a deliberately chosen spectral codeword; its low-description-complexity origin is disclosed. No GPU benchmark is included.
## Scope of the included CI

The GitHub Actions job runs the 7 algebra/codec regression tests, checks a small CLI round trip, and reproduces the 1,048,576-dimensional synthetic codeword. The larger cross-corpus benchmarks are intentionally opt-in via `sh run_full_suite.sh` because they are resource-intensive and can depend on optional rendered paper-page fixtures. The code is experimental: no GPU run, no universal compression gain, and no factorization speedup are claimed.

Research results are in [the full report](hypercomplex_highdim_full_report.md) and frozen [benchmark CSVs](results/).

## Held-out SOTA compression benchmark

The [compression benchmark](sota_benchmark/README.md) adds a **file-disjoint calibration/held-out protocol** to the original high-dimensional research: Zstandard levels 1/3/9/19 and calibration-trained dictionaries, Brotli up to quality 11, LZMA, a per-block hybrid with exact framing and CRC, encode/decode throughput and sampled memory, and ablations against ordinary Walsh-Hadamard, XOR prediction, dictionary coding, and a specialized basis-index codebook.

Reproduce the complete suite with `sh sota_benchmark/run.sh full`; CI runs `sh sota_benchmark/run.sh smoke`-equivalent checks and uploads its frozen configuration and results as an artifact. The full six-file run is available as an explicit GitHub Actions workflow-dispatch option, not part of every PR CI run. This benchmark does **not** claim broad SOTA competitiveness: previous synthetic basis-aligned results are evaluated separately from truly file-disjoint source-code and research-text heldout inputs.

## Causal hypercomplex prediction experiment

The new [predictive multiplication laboratory](predictive_lab/README.md) tests a different arithmetic hypothesis: whether the **actual Cayley-Dickson product of two previously decoded data blocks** predicts the next block well enough to compress its residual. Unlike the signed spectral cascade, this adds a genuinely nonlinear past-state interaction. The model is trained only from an initial prefix *transmitted in the file*, with quantized parameters and framing fully charged to the final size. Its controls include linear AR2, XOR, unsigned XOR convolution, coordinatewise and cyclic bilinear products, and direct Brotli11/Zstd9/LZMA6.

**Outcome:** none of the tested natural code/text/PDF data beat direct Brotli11; the causal nonlinear CD predictor typically produced *less compressible* residuals than ordinary AR2. A separate **deliberately CD-generated** synthetic control did benefit (12,335 B CD vs 12,609 B AR2 vs 12,754 B direct Brotli11), demonstrating representational feasibility rather than general compression superiority. The [full analysis](predictive_lab/RESULTS_2026-11-18.md), [raw benchmarks](predictive_lab/results/) and [dimension sweep](predictive_lab/sensitivity.py) are reproducible.

Run `python -m unittest -v predictive_lab.test_predictor` and `python -m predictive_lab.benchmark --smoke --out /tmp/hc-predictive-smoke` from this directory. CI runs these automatically; full benchmarks are available through workflow dispatch. No GPU acceleration or natural-source SOTA advantage is claimed.
