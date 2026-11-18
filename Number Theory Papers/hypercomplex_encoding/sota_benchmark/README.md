# Held-out SOTA compression benchmark (PR #5 extension)

This is an **experimental, lossless** benchmark for the Cayley-Dickson signed cascade in the parent directory. The earlier impressive basis-codeword result was **in-distribution**: input data were constructed from this transform. This extension evaluates whether the same mechanism carries over to independently selected real inputs, with a reproducible train/held-out split and conservative accounting.

## Five required comparisons

1. **File-disjoint held-out data.** `corpus_manifest.json` lists calibration data (`realpair_text.txt` and a preexisting factorization Python program) and five separately identified held-out source files (two algebra/codec modules, the independent C++ hybrid factorizer, and two independent research notes). In normal CI, all five are required and the loader records the SHA-256 of **both the complete source and the evaluated bytes**. Sources are tracked in the Number-Theory repository. With `--smoke`, each file is capped at 3,072 bytes, an intentionally limited CI check; full mode uses all bytes. Neither dictionary training nor hyperparameter calibration opens the held-out files. The test is a limited real-source corpus, **not a representative cross-domain compression benchmark**. Bring new, independently collected bitmask/telemetry/media data by editing the manifest **before** a new evaluation; no post-hoc cherry-picking.
2. **Strong conventional baselines.** Whole-file `zstd1`, `zstd3`, `zstd9`, `zstd19`, `brotli5`, `brotli9`, `brotli11`, `lzma6`, `lzma9`, and `zlib6`/`zlib9` (zlib excluded from CI smoke to bound runtime). **Zstandard dictionary** `zstd_dict3`/`zstd_dict9` trains on calibration data only using libzstd COVER; a standalone dictionary frame embeds its full dictionary **in each file**, including framing, rather than assuming a free shared dictionary. Dictionary-training time and dictionary size are reported separately. The largest allowable zstd level may vary by libzstd version: our chosen levels are supported by libzstd 1.5.x. No internet access is required after `libzstd1`, `numpy`, `numba`, `brotli`, `psutil`, `Pillow` are installed.
3. **Framed per-block hybrid.** Every file contains format/version, block size and algebra dimension, full input length, CRC32, and per-block mode tags and varint payload lengths. Each block chooses the smallest complete serialized representation among raw bits, zero/one fills, `zstd3`, `brotli5`, and (when selected) genuinely hypercomplex sparse coefficients. Encoding searches all allowed candidates per block and its time **is included**. The bitstream is split at *D-1-bit boundaries*, not prematurely rounded to bytes, to preserve the actual imaginary-coordinate encoding. Both CRC checks and exact byte-level equality are required.
4. **Throughput and memory.** Each method measures median **full-file** encoding and decoding time (warmup excluded), MB/s, and a separate memory replay reporting sampled process RSS and peak Python allocations. RSS sampling every millisecond misses brief native peaks and shared Python/Numba runtime is not attributed per codec; figures are diagnostics, not a machine-independent memory ceiling. Record hardware, versions, source SHA-256, exact bytes, and all measurements in JSON and CSV. Calibration configuration search and dictionary training times are **explicit** and not hidden in the measured frozen per-file invocation.
5. **Ablation and specialized competitors.** Hold the same per-block raw/zstd/Brotli fallbacks fixed and replace the spectral step with (a) standard Walsh-Hadamard, (b) bytewise XOR prediction, (c) no spectral transform, and (d) a restricted basis-index spectral mode. Additionally compare a standalone **basis-family codebook-index** format (only applicable when an entire D-1-bit block matches a single constructed Cayley-Dickson basis pattern). The latter has no hidden shared dictionary, serializes indices/signs and final tail; it is deliberately family-specific and demonstrates that hypercomplex's 29-byte constructed signal does **not** establish an information-theoretic or specialized-compressor breakthrough. Zstandard+trained dictionary is a separate conventional dictionary ablation.

## Anti-leakage protocol

There are **two separate commands**. `calibrate` examines only the development files, trains the Zstandard dictionary, measures the search time, and freezes each ablation's one global block size using the minimum *aggregate calibration byte count* (time breaks ties). Its JSON includes the exact train-file SHA-256s and base64 serialized dictionary. `evaluate` verifies those hashes, opens only the held-out sources, and applies the fixed configurations without retuning by file. The per-file *minimum among many conventional compressors* is reported as an **oracle diagnostic**, clearly labeled; it is not an implementable zero-cost codec selector. Results for constructed basis-family inputs have `label=synthetic_diagnostic` and are never included in natural-data conclusions.

All methods must decode bit-for-bit. This harness uses small research-grade files and is not an official SOTA benchmark, an objective cross-platform throughput ranking, or a GPU test.

## Run

From `Number Theory Papers/hypercomplex_encoding`:

```sh
python -m pip install -r requirements.txt
python -m unittest -v sota_benchmark.test_codec sota_benchmark.test_protocol
python -m sota_benchmark.benchmark calibrate --manifest sota_benchmark/corpus_manifest.json --config /tmp/hc-frozen.json
python -m sota_benchmark.benchmark evaluate --manifest sota_benchmark/corpus_manifest.json --config /tmp/hc-frozen.json --output /tmp/hc-sota-results --synthetic --repeats 3
```

CI's bounded version:

```sh
python -m sota_benchmark.benchmark calibrate --smoke --config /tmp/hc-frozen.json
python -m sota_benchmark.benchmark evaluate --smoke --config /tmp/hc-frozen.json --output /tmp/hc-sota-smoke --synthetic --repeats 2
```

`--allow-missing` exists **only** for standalone experiments with an incomplete local repository snapshot. The CI suite does **not** set this option. Full benchmark outputs (JSON/CSV/Markdown) are stored as GitHub Actions artifacts. Pass `--no-memory` for a quicker size-and-throughput-only exploratory run, marking the memory fields absent. `--max-bytes` is optional, but any truncated held-out run must be reported as such and is not a full-corpus result.

### Format and benchmark limitations

The hybrid codec's spectral payload stores signed coefficient values and delta-encoded coordinate indices; its cost is charged in full. The test uses the Cayley-Dickson `cd_cascade` implementation from the main research suite and a separately implemented Walsh-Hadamard transform. Both are mathematically invertible, but sparse compressibility is distribution-dependent. Nonassociativity/zero divisors in dimensions `>=16` are not circumvented: the exact inverse is the tested sequence of basis-left-multiplication **linear operators**, not division by a general hypercomplex element. CRC detects most accidental corruption but does not offer cryptographic integrity or error correction. Dictionary sizes may dominate small files, and native low-level compression libraries are expected to be faster than interpreted whole-file block selection.

See `benchmark.py`, `codec.py`, `test_codec.py`, `test_protocol.py` and machine-readable `results.json` / `measurements.csv` from the workflow. No claim of SOTA competitiveness is made without a broad, truly independent benchmark.
