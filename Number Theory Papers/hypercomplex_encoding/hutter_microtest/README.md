# Hutter-style quaternion context prediction microbenchmark

**Scope:** Small, falsifiable test of whether quaternion-derived byte features improve conditional probability estimates of a baseline order-4 byte model. **Not** a Hutter Prize entry or a deployable lossless compressor.

## Corpus

Use the unmodified 100MB enwik8 corpus from Matt Mahoney's Large Text Compression Benchmark, extracted as enwik8. The workflow downloads enwik8.zip, checks the entire 100,000,000-byte corpus against SHA-1 57b8363b814821dc9d47aa4d41f58733519076b2, and evaluates the first 380,000 bytes. This is a byte-identical subset of the Wikipedia XML snapshot underlying enwik9.

## Comparison

The byte stream is not modified. All predictors use **only previously decoded bytes**. The baseline is adaptive order-4 character n-gram interpolation/backoff. Hamilton feature: Hamilton product of the last eight bytes interpreted as two signed 4-component quaternions, quantized to 512 categories. Comparators: paired-XOR, paired-addition, deterministic 8-byte hash; exactly 512 categories apiece.

First 200k bytes warm up adaptive models, next 60k select each feature's mixture weight from the preregistered grid 0/.05/.1/.2/.4/.8/1, last 120k are held out. Model statistics continue causal updates within heldout. Hyperparameters other than selected mixture weights are fixed before seeing the validation/test results.

**Reported values are ideal minus-log2 conditional code lengths, not actual compressed streams.** No entropy decoder or compressor binary is included and model/code costs are omitted. The standalone Brotli-11 and zlib-9 output sizes on the same held-out segment are reference controls, not directly comparable to externally warmed predictive codelengths.

Run: python -m unittest -v test_benchmark ; python benchmark.py --input enwik8 --output results

A quaternion-specific improvement must survive comparisons with ordinary byte contexts **and** equally sized XOR/additive/hash features, on held-out data. The Hutter Prize evaluates complete lossless enwik9 archive plus the decoder, not ideal log loss.