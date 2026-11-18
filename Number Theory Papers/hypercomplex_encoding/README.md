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
