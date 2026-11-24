# Full high-dimensional imaginary-coefficient encoding experiment

**Run date:** November 18, 2026  
**Source:** Fred Viole, *Real-Pair Representations of Complex Numbers and Quaternions* (November 4, 2026)  
**Status:** Fully implemented, exact arithmetic through specified tested dimensions; reproducible tests and lossless file codecs included.

## Executive findings

- **Yes, full-dimensional encoding is mathematically and computationally possible.** The Cayley-Dickson basis contains one real axis and D-1 independent imaginary axes. Exact bit recovery and reversible imaginary-unit multiplication passed at all 14 dimension levels through **D=1,048,576**. Multi-bit coefficients of widths **1, 2, 4, 8, and 16** passed independently through the same maximum dimension.
- **Your real-pair transformation is literally multiplication by 1+i.** `(1+i)(a+bi)=(a-b)+i(a+b)`. For quaternion coefficients, the same left multiplication produces the paper's blockwise `psi(q)=(a-b,a+b,c-d,c+d)`.
- **A true hypercomplex mixing transform also works.** Compose successive left multiplications by `(1+e_1), (1+e_2), (1+e_4), ...`. Each stage has an exact inverse on the integer image; the resulting signed Hadamard-type matrix satisfies `M_D^T M_D=D I` and was tested to a million dimensions. Beyond the octonions this is a succession of linear operators, not generally a single associative hypercomplex product.
- **Strong structure-dependent compression exists.** A 1,048,575-bit synthetic codeword aligned with one hypercomplex basis pattern encoded to **29 bytes** with exact decoding. Direct Brotli quality 9 needed **274 bytes**; raw input was **131,072 bytes** including one byte-padding bit.
- **No general-purpose compression win was found.** None of the tested ordinary text, source code, PDF, document bitmap, or random data beat the best direct zlib/Brotli/LZMA baseline. The special wins use data constructed from the transform basis.
- **Algebraic arithmetic deteriorates with dimension.** The sedenion example `(e1+e10)(e4-e15)=0` establishes the existence of zero divisors. Arbitrary nonzero hypercomplex multiplication is unsafe as a reversible key; the tested butterfly sequence avoids that problem.

## 1. Mathematical definitions and the connection to the paper

Let `D=2^k` and `e_0=1`, with imaginary axes `e_1,...,e_(D-1)`. Encode a binary message of D-1 bits by `X=sum_{i=1}^{D-1} b_i e_i`, `b_i in {0,1}`, or alternatively by `X=e_0+sum_{i=1}^{D-1} (2b_i-1)e_i` for the bipolar model. The latter uses the real coefficient as a fixed +1 reference, **not** as a free information bit.
Cayley-Dickson multiplication uses `(a,b)(c,d)=(ac - conjugate(d)b, da + b conjugate(c))` at each doubling. Basis multiplication has `e_i e_j = s(i,j) e_(i XOR j)` with `s(i,j)` in `{+1,-1}`; the code calculates the sign recursively, independently verifies 21,844 basis-product pairs, and uses Numba for actual multiplication.
The original two-coordinate real-pair map is `T(a,b)=(a-b,a+b)`. Multiplication by `1+e_1` realizes this map exactly at dimension 2 and the paper's quaternion blockwise map at dimension 4. Let `L_j(x)=e_j*x`. For tested Cayley-Dickson basis units, `L_j^2=-I`, so each `B_j=I+L_j` has `B_j^T B_j=2I` and inverse `(I-L_j)/2`. The product `M_D=B_(D/2)...B_4 B_2 B_1` has `M_D^T M_D=D I`, and can be inverted by reversing the stages, dividing by two at each exact stage. The code verifies this integrality throughout the round trip.
For D<=64 the entire matrix was independently constructed from basis columns and checked to contain only +/-1 entries and satisfy the Gram identity exactly. Larger dimensions were tested with exact integer roundtrips and norm identities. The construction is a signed Hadamard-type transformation *derived from actual higher-dimensional imaginary-basis multiplication*; it is not a claim of new information or superior asymptotic computation.

## 2. Exact encoding and dimension scaling

| Total algebra dimension D | Imaginary axes | Tight one-bit storage | Dense int8 coefficients | Unit multiply time (ms)* |
|---:|---:|---:|---:|---:|
| 4 | 3 | 1 B | 4 B | 0.001 |
| 8 | 7 | 1 B | 8 B | 0.001 |
| 16 | 15 | 2 B | 16 B | 0.001 |
| 256 | 255 | 32 B | 256 B | 0.003 |
| 4,096 | 4,095 | 512 B | 4,096 B | 0.096 |
| 65,536 | 65,535 | 8,192 B | 65,536 B | 1.651 |
| 1,048,576 | 1,048,575 | 131,072 B | 1,048,576 B | 28.968 |

*Median of repeated Numba-compiled CPU runs; timings are prototype diagnostics, not a controlled cross-machine benchmark or a GPU result. The fixed unit is `e_7` where available. All values passed exact bit recovery.

The million-dimensional 1-bit case represents **1,048,575 binary coefficients**, requiring **131,072 packed bytes**. Materializing one signed 8-bit coefficient per axis costs **1,048,576 bytes**; 64-bit floating coefficients would cost **8,388,608 bytes**. These are storage representations of the *same information*.

## 3. Increasing the bits stored in each imaginary coefficient

| Bits per imaginary coefficient | Information capacity in 1,048,575 axes | Tight byte count | Unit multiplication/inversion |
|---:|---:|---:|---|
| 1 | 1,048,575 bits | 131,072 B | Exact / exact |
| 2 | 2,097,150 bits | 262,144 B | Exact / exact |
| 4 | 4,194,300 bits | 524,288 B | Exact / exact |
| 8 | 8,388,600 bits | 1,048,575 B | Exact / exact |
| 16 | 16,777,200 bits | 2,097,150 B | Exact / exact |

Using 16-bit coefficients and D=1,048,576, the full successive imaginary mixing transform also inverted **exactly**. Its largest output coefficient in the seeded test had magnitude **129,324,076**, which exceeds the original 16-bit coefficient range: transported arithmetic requires adequate intermediate precision. The exact norm identity was checked with Python arbitrary-size integer summation to avoid overflow of a 64-bit squared-energy accumulator.

### Direct multiplication versus signed-unit scrambling

- Dense binary-coefficient Cayley-Dickson multiplication was executed through **D=4,096**; the D=4,096 case required 212.088 ms in this CPU prototype, with 4,388,944 nonzero term products, and yielded 4,046 nonzero output coefficients.
- True sparse-sparse multiplication was executed up to **D=1,048,576** with 64 nonzero input coordinates on each side: 4,096 term products, 4,085 nonzero outputs, 2.703 ms.
- Fully dense multiplication of two million-dimensional vectors **was not attempted**; its straightforward algorithm would require on the order of a trillion basis term pairs. The million-dimensional *unit multiplication* and *multi-stage mixing* were actually executed and verified.

## 4. Algebra verification and failure modes

| Algebra dimension | Property | Result |
|---:|---|---|
| 2 (complex) | Associative, multiplicative norm | Confirmed |
| 4 (quaternions) | Associative, multiplicative norm | Confirmed |
| 8 (octonions) | Nonassociative, multiplicative norm | Nonassociativity witnessed by `(e1 e2)e4` versus `e1(e2 e4)` |
| 16 (sedenions) | Zero divisors and failed composition law | `(e1+e10)(e4-e15)=0`, both operands nonzero |
| 32, 64, 128 | Multiplicative norm generally fails | Counterexamples found |

Randomized integer coefficient checks of norm multiplicativity:

| D | Trials | Failures of N(xy)=N(x)N(y) |
|---:|---:|---:|
| 2 | 100 | 0 |
| 4 | 100 | 0 |
| 8 | 100 | 0 |
| 16 | 100 | 92 |
| 32 | 100 | 99 |
| 64 | 100 | 100 |
| 128 | 40 | 40 |

These random tests demonstrate counterexamples, not a statistical property of all inputs. The algebraic zero-divisor example is an explicit exact witness, not a probabilistic finding.

### Associative Clifford alternative

For `Cl_(0,n)` with n anticommuting imaginary generators, storing an n-bit **vector** is straightforward, but arbitrary multiplication leaves the vector subspace. The full algebra has **2^n** basis coordinates; only the scalar and grade-2 components of a vector product were materialized. At n=1,024 generator bits, there are **523,776 possible bivectors**, of which **205,368** were nonzero in one seeded random product. This avoids the Cayley-Dickson nonassociativity, but is not an ordinary fixed-size n-vector arithmetic format under multiplication.

## 5. Actual file compression: fair end-to-end comparisons

The codec applies the algebraically defined exact multi-stage transform to chunks with D-1 imaginary message coefficients, supports **binary (0/1)** and **bipolar (-1/+1)** coefficient models, and chooses raw/fill/sparse-spectral block encodings. Each output includes block tags, dimensional metadata, original length, and CRC32. The completed encoded file is compressed with the smaller of zlib 9 and Brotli 9. Each dataset was exhaustively tested over block dimensions 256, 1,024, 4,096, and 65,536, selecting the smallest result. For conventional comparison, the original bytes were compressed directly with zlib 9, Brotli 9, and LZMA 6. Search over all candidate configurations is not included in individual encoding-time figures.

| Data | Original bytes | Best direct compressor | Direct size | Best hypercomplex cascade | Result |
|---|---:|---|---:|---:|---|
| Paper extracted text | 17,582 | brotli9 | 5,636 B | 6,930 B | 23.0% larger |
| Python source | 17,634 | brotli9 | 5,395 B | 6,359 B | 17.9% larger |
| Original paper PDF | 355,787 | brotli9 | 311,547 B | 316,023 B | 1.4% larger |
| Document 1-bit bitmap | 23,808 | lzma6 | 2,284 B | 2,318 B | 1.5% larger |
| Random bytes | 65,536 | brotli9 | 65,540 B | 65,584 B | 0.1% larger |
| Synthetic 1,024D basis codewords | 65,472 | lzma6 | 10,016 B | 1,133 B | 88.7% smaller |
| Synthetic 4,096D basis codewords | 65,520 | lzma6 | 5,272 B | 373 B | 92.9% smaller |

**Interpretation:** For the five listed real/random inputs, no sparse spectral blocks were selected by the best codec. The special inputs are algorithmically constructed from this *exact same hypercomplex transformation*: one spectral coefficient per block, so the codec dramatically outperforms the general compressors tested. These are in-sample signal-model demonstrations, not independent natural data or evidence of universal compression.

## 6. One-million-dimensional binary codeword (independent demonstration)

- Generated a signed hypercomplex codeword with D=1,048,576, selected basis index e_123,456, and fixed real coefficient +1.
- The source file holds D-1 meaningful imaginary bit coefficients, packed into 131,072 bytes (one trailing padding bit).
- The transformed binary coefficients yield a single nonzero spectral coefficient; the decoder also handles the single padded bit through a fill block.
- Encoding, inverse recovery, and CRC32 verification passed; SHA-256 is `c383021767df4f09ce4dc64d709f9d99868bc8d97cfae7ed1b5186c5585c6c79`.

| Method | Total bytes |
|---|---:|
| Original packed source | 131,072 |
| zlib level 9 | 898 |
| Brotli quality 9 | 274 |
| LZMA preset 6 | 844 |
| Raw hypercomplex spectral file, including metadata | 29 |
| Hypercomplex + Brotli | 33 |

A tighter specialized **basis-index** encoder can represent a known codeword from a D-sized codebook using only `log2(D)` bits for its index (plus framing and protocol metadata). In this sense the above result is not an information-theoretic breakthrough: its input was explicitly chosen from a low-description-complexity family. The algebra gives a principled way to *recognize and transform* that family, not a way to encode arbitrary random strings in fewer bits than their entropy.

## 7. Error behavior

1,000 independent bit-flip trials in D=4,096 signed-unit bit-packed representations showed exactly one recovered data-bit error per one transmitted bit flip, when the flipped coordinate corresponded to a message bit. No native error correction was observed. CRC32 makes corrupted files detectable, but CRC alone cannot repair them. Multi-stage spectral transforms are exactly invertible in noiseless integer arithmetic; neither they nor the algebra change the underlying channel-capacity constraint.

## 8. What has and has not been established

**Established:** actual high-dimensional imaginary coefficients can hold and recover independent bits or multi-bit digits; products by basis units and sequential (1+e_j) transformations are implemented with correct signs; the tested mixing map has an exact inverse, a scaled orthogonal Gram identity, and exceptional sparsity on particular algebraically constructed sign patterns; runnable file codecs and corruption tests exist.

**Not established:** general compression superiority; extra information capacity without paying coefficient precision/storage; GPU acceleration; arithmetic speedup on random large multiplications; safe use of every nonzero high-dimensional key; universal error correction; or a novel entropy lower bound. The dense multiplication test stops at D=4,096, while the million-dimensional tests include real unit products, sparse products, and full sequential mixing—not fully dense binary-by-binary multiplication.

## 9. Reproduction

Dependencies: Python 3.10+ (tested with current Python 3.13), `numpy`, `numba`, `brotli`, and `Pillow`. Tests use fixed seeds and do not require network access. From the extracted package folder:

```bash
python -m unittest -v hypercomplex_suite_tests
python hypercomplex_highdim_test.py
python hypercomplex_cascade_test.py --phase codec
python hypercomplex_cascade_test.py --phase spectrum
python hypercomplex_widths_and_errors.py
python hypercomplex_million_dim_demo.py
python hypercomplex_combine_results.py
```

Try a user file (lossless encode/decode):

```bash
python hypercomplex_codec_cli.py encode INPUT.bin OUTPUT.hcb --mode cascade --dimension 4096 --coefficient bipolar --compress brotli
python hypercomplex_codec_cli.py decode OUTPUT.hcb RESTORED.bin
cmp INPUT.bin RESTORED.bin
```

The results are single-machine CPU/Numba prototype measurements, with short repeated timing samples. The compressor tests involve whole-file oracle selection across known configurations, so comparisons indicate attainable compressed sizes rather than training-free online throughput. No GPU kernel was run.

## 10. Next experiment suggested by the full-dimensional result

Use the **bipolar** high-dimensional transform on a measured data source with a hypothesized low-support generator (e.g., actual structured bitmask streams or sparse spectral code families), and evaluate on held-out messages with a full encoder-format overhead and speed comparison against Brotli, LZMA, and a specialized codebook-index compressor. Only after confirming a real data distribution that is sparse in this basis is SIMD/GPU implementation likely to be worthwhile.
