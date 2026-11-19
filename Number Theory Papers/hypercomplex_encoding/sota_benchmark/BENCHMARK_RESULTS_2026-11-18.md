# First full held-out SOTA compression results — November 18, 2026

**Protocol:** Calibrate on two development files, freeze one block size per codec and a development-only Zstandard dictionary, then evaluate six disjoint full-file held-out inputs plus two clearly labeled synthetic diagnostics. Exact byte reconstruction, full encoded size and CRC framing, timing, sampled RSS, dictionary byte cost and metadata are enforced.

**Reproduction:** [GitHub Actions run 44390105121](https://github.com/OVVO-Financial/Number-Theory/actions/runs/44390105121) (all three jobs passed), artifact [hypercomplex-sota-full (14468565424)](https://github.com/OVVO-Financial/Number-Theory/actions/runs/44390105121/artifacts/14468565424). This is a snapshot of the branch at the indicated run, not an independent, multi-machine SOTA leaderboard.

## Full-file holdout comparison

| Dataset | Input bytes | Best standard compressed bytes | Frozen hypercomplex bytes | Hypercomplex size change |
| --- | ---: | ---: | ---: | ---: |
| Hypercomplex algebra source | 19,987 | 5,788 (Brotli-11) | 12,174 | +110.33% |
| Hypercomplex cascade source | 15,689 | 4,705 (Brotli-11) | 9,824 | +108.80% |
| Independent C++ factorizer | 125,209 | 32,775 (Brotli-11) | 75,864 | +131.47% |
| Hybrid research notes | 9,693 | 3,598 (Brotli-11) | 7,052 | +96.00% |
| Complex-space method notes | 26,155 | 9,177 (Brotli-11) | 18,417 | +100.69% |
| Separate Fermat-sieve PDF binary | 163,135 | 156,120 (Brotli-11) | 163,982 | +5.04% |

All twelve independent SOTA codec/split tests and seven original algebra/codec regression tests passed CI; the six complete natural files decoded bit-for-bit.

**Key ablation:** on *every one* of the six held-out natural files, the frozen hypercomplex and fixed-block Walsh, XOR-predictive, no-spectral, and basis-index alternatives produced **identical final byte sizes**. The hypercomplex chooser selected only `brotli5`, `zstd3`, or `raw` modes, and **zero spectral blocks**. For example, the 125,209-byte C++ file used 204 Brotli-5 blocks and 41 Zstandard-3 blocks. The spectral transform therefore contributed no observed size saving in these six sources.

## Deliberately constructed family — not held-out natural data

| Constructed input | Original | Best of conventional tested | Frozen hypercomplex | Standalone basis-index codebook |
| --- | ---: | ---: | ---: | ---: |
| Basis family at D=1,024 (cross-dimension) | 3,069 B | 816 B (Zstd-19) | 1,135 B | 91 B |
| Basis family aligned to frozen D=4,096 | 6,143 B | 648 B (Brotli-11) | 88 B | **61 B** |

The constructed aligned-family result demonstrates exact high-dimensional sparsity, but the specialized codebook outperforms the more general hypercomplex format. This does **not** establish universal compression superiority or a novel entropy bound.

## Throughput and measurement accounting

The Python/Numba hypercomplex encoder achieved approximately **0.42–0.43 MB/s** on these natural sources (GitHub Actions runner); decoding was approximately **34–44 MB/s**. This is far slower at encoding than the optimized native Zstandard implementations. Brotli-11 is slower than low-level Zstd but compressed each of these natural inputs more tightly. These numbers include the online per-block mode search but exclude JIT warmup; the configuration grid search and Zstandard dictionary-training time are reported separately.

Frozen configuration: **D=4,096** (`block_bytes=512`), selected using only calibration sources. The trained dictionary was 4,096 bytes and its entire byte content was included once per compressed file in dictionary-backed baselines. Search time was **1,181.449 ms** and dictionary-training time **10.699 ms** in this CI run. RSS results are sampled process-wide and not independently attributable to individual codecs; `measurements.csv` also records peak Python allocations.

## Research interpretation

This benchmark **does not show SOTA competitiveness on ordinary inputs**. It demonstrates a valid reversible hypercomplex encoding and a useful special-family sparsity phenomenon, but no size win on the present natural holdout corpus and no speed advantage. Only six research-oriented files (mostly code/text, plus one PDF) are tested; true external telemetry, bitmask streams, media, independent structured-scientific datasets, and GPU implementations remain untested.

The reproducible code, corpus manifest, calibration and evaluation commands are in [`sota_benchmark/README.md`](README.md). The complete per-file/per-codec JSON/CSV and frozen calibration configuration are attached to the cited Actions run. The full 178-row run includes four fixed-block causal ablations as well as separately calibrated versions.
