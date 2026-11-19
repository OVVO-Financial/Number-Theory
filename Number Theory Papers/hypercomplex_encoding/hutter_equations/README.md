# Direct algebraic-equation discovery for the Hutter text corpus

**Question:** Do algebraic equations among previous bytes constrain the next byte of naturally occurring Wikipedia text enough to reduce lossless coding cost, especially equations from complex/quaternion multiplication?

**This is not the old context-label experiment.** Each of 30 prespecified rules explicitly predicts a byte `pred_t = f(byte_(t-1), ..., byte_(t-16))`, with operations modulo 256. A residual `r_t = (actual_t - pred_t) % 256` permits exact reconstruction. Rules cover constant-character controls, copy/equality, linear difference, XOR, quadratic, actual complex multiplication and Hamilton quaternion multiplication. The scalar and vector equations are tested separately; their arithmetic is independent of context hashing.

## Measurements

1. **Ideal adaptive log-loss:** A causal order-4 interpolated byte context predictor gives `p0(next_byte | past)`. Each algebraic equation gives a normalized prediction distribution by tracking the conditional counts of its observed residuals in 8 fixed preceding-byte categories. A validation-selected convex mixture of these two distributions is scored on the heldout suffix. Every symbol, including heldout symbols, enters the model *after* that symbol's score, which an entropy decoder can reproduce. Each formula has identical residual-context table capacity.
2. **Real lossless sizes:** The heldout suffix is encoded as exact mod-256 residual bytes using the selected formula, compressed with Brotli quality 11, and placed into a standalone frame with model ID, raw seed, CRC32, length, and (when used) an 8-byte learned class map. Decode is byte-for-byte checked. This is a different measurement from ideal probabilistic bits; the ideal mixture does **not** have a range coder here.
3. **Multi-equation tests:** Validation chooses the strongest ordinary, complex and quaternion rules; a *second* convex-mixture comparison checks whether quaternion equations add benefit beyond the ordinary source-derived constraint. This helps rule out merely relabeling an ordinary equality/difference pattern.
4. **Validation-chosen class map:** Among the 30 models, choose the residual-bias-favored equation separately for preceding-byte classes (lowercase, uppercase, digit, whitespace, XML punctuation, control, non-ASCII, other). Store all eight choices explicitly in the independent lossless frame. This is not a free oracle on heldout data.

## Source and split

The CI workflow downloads the **official enwik8 100,000,000-byte file** from `https://mattmahoney.net/dc/enwik8.zip`; verifies full-file SHA-1 `57b8363b814821dc9d47aa4d41f58733519076b2`; and runs **two independent offset windows** at 0 and 2,000,000 bytes, each with 200k warmup, 60k validation and 120k heldout bytes. The holdout is chronological and was not used to select equation variants, λ values or class mappings. The formula pool and class rules are preregistered in source code. Input/output hashes go into `result.json`.

The result is small-scale evidence about enwik8, which is related to but not the full 1GB enwik9 Hutter Prize task. Two windows are not sufficient for leaderboard claims. Decoder executable size, full trainer, and a combined probabilistic coder are **not** included in the ideal-number comparison, and the framed codec must be compared on compressed bytes instead.

## Reproduce

```bash
python -m pip install numpy brotli
python -m unittest -v test_equations
curl -fsSL https://mattmahoney.net/dc/enwik8.zip -o /tmp/enwik8.zip
unzip -p /tmp/enwik8.zip enwik8 > /tmp/enwik8
printf '57b8363b814821dc9d47aa4d41f58733519076b2  /tmp/enwik8\n' | sha1sum -c -
python benchmark.py --input /tmp/enwik8 --offsets 0 2000000 --source-label official-verified-enwik8 --output results
```

`results/offset_0/` and `results/offset_2000000/` contain JSON, all equation-level CSV rows and a Markdown summary. `results/SUMMARY.json` provides the concise cross-offset outcomes. GitHub Actions publishes the files as an artifact.

## Interpretations and exclusions

- An equation's exact-hit rate matters. If its direct predicted byte is usually wrong, an improvement in conditional probability may come from the residual histogram, not an exact arithmetic relation in language.
- Constant-space and constant-`e` are explicit negative controls. More expressive algebra must improve **beyond** such controls, matched copy/add/XOR rules, and Brotli, after all metadata.
- A short regular input generated to satisfy one of the equations confirms the codec's capability, not source-specific natural-language structure.
- XML and linguistic grammars may have genuine constraints, but this particular experiment does not implement an XML parser or an external language model; it tests this finite family of arithmetic equations only.
- No GPU, enwik9 score, or Hutter Prize competitiveness is claimed.
