#!/usr/bin/env python3
"""Aggregate verified phases and author the full-dimensional research report."""
from pathlib import Path
import json, platform, csv
R=Path(__file__).resolve().parent
load=lambda f:json.loads((R/f).read_text())
unit=load('hypercomplex_highdim_results.json')
spectral=load('hypercomplex_cascade_spectrum_only.json')
codec=load('hypercomplex_cascade_bipolar_results.json')
width=load('hypercomplex_widths_results.json')
million=load('hypercomplex_million_dim_codeword_results.json')
full={'source_paper':'Fred Viole, Real-Pair Representations of Complex Numbers and Quaternions (2026-11-04)',
      'date':'2026-11-18',
      'definition':'Actual Cayley-Dickson algebra products; left multiplication by basis imaginary units; sequential (1+e_2^r) multiplications and their exact inverses; Clifford grade-1 alternative.',
      'unit_test':unit,'signed_cascade_spectrum':spectral,
      'signed_cascade_codec':codec,'coefficient_widths':width,
      'million_dim_codeword':million}
(R/'hypercomplex_full_results.json').write_text(json.dumps(full,indent=2)+'\n')

def md_num(x): return f'{x:,}'
def fmt_ms(x): return f'{x:.3f}'
dim_set=[x for x in unit['dimensions'] if x['dimension_D'] in (4,8,16,256,4096,65536,1048576)]
norm_set=unit['algebra']['norm_trials']
comp_set=codec['compression']
consts={x['dataset']:x for x in comp_set}
wide=[x for x in width['multi_bit_widths'] if x['dimensions']==1048576]
sp_million=[x for x in spectral['spectrum_scaling'] if x['dimension']==1048576]
col=[]
col.append('# Full high-dimensional imaginary-coefficient encoding experiment\n')
col.append('**Run date:** November 18, 2026  \n**Source:** Fred Viole, *Real-Pair Representations of Complex Numbers and Quaternions* (November 4, 2026)  \n**Status:** Fully implemented, exact arithmetic through specified tested dimensions; reproducible tests and lossless file codecs included.\n')
col.append('## Executive findings\n')
col.append('- **Yes, full-dimensional encoding is mathematically and computationally possible.** The Cayley-Dickson basis contains one real axis and D-1 independent imaginary axes. Exact bit recovery and reversible imaginary-unit multiplication passed at all 14 dimension levels through **D=1,048,576**. Multi-bit coefficients of widths **1, 2, 4, 8, and 16** passed independently through the same maximum dimension.')
col.append('- **Your real-pair transformation is literally multiplication by 1+i.** `(1+i)(a+bi)=(a-b)+i(a+b)`. For quaternion coefficients, the same left multiplication produces the paper\'s blockwise `psi(q)=(a-b,a+b,c-d,c+d)`.')
col.append('- **A true hypercomplex mixing transform also works.** Compose successive left multiplications by `(1+e_1), (1+e_2), (1+e_4), ...`. Each stage has an exact inverse on the integer image; the resulting signed Hadamard-type matrix satisfies `M_D^T M_D=D I` and was tested to a million dimensions. Beyond the octonions this is a succession of linear operators, not generally a single associative hypercomplex product.')
col.append('- **Strong structure-dependent compression exists.** A 1,048,575-bit synthetic codeword aligned with one hypercomplex basis pattern encoded to **29 bytes** with exact decoding. Direct Brotli quality 9 needed **274 bytes**; raw input was **131,072 bytes** including one byte-padding bit.')
col.append('- **No general-purpose compression win was found.** None of the tested ordinary text, source code, PDF, document bitmap, or random data beat the best direct zlib/Brotli/LZMA baseline. The special wins use data constructed from the transform basis.')
col.append('- **Algebraic arithmetic deteriorates with dimension.** The sedenion example `(e1+e10)(e4-e15)=0` establishes the existence of zero divisors. Arbitrary nonzero hypercomplex multiplication is unsafe as a reversible key; the tested butterfly sequence avoids that problem.')
col.append('\n## 1. Mathematical definitions and the connection to the paper\n')
col.append('Let `D=2^k` and `e_0=1`, with imaginary axes `e_1,...,e_(D-1)`. Encode a binary message of D-1 bits by `X=sum_{i=1}^{D-1} b_i e_i`, `b_i in {0,1}`, or alternatively by `X=e_0+sum_{i=1}^{D-1} (2b_i-1)e_i` for the bipolar model. The latter uses the real coefficient as a fixed +1 reference, **not** as a free information bit.')
col.append('Cayley-Dickson multiplication uses `(a,b)(c,d)=(ac - conjugate(d)b, da + b conjugate(c))` at each doubling. Basis multiplication has `e_i e_j = s(i,j) e_(i XOR j)` with `s(i,j)` in `{+1,-1}`; the code calculates the sign recursively, independently verifies 21,844 basis-product pairs, and uses Numba for actual multiplication.')
col.append('The original two-coordinate real-pair map is `T(a,b)=(a-b,a+b)`. Multiplication by `1+e_1` realizes this map exactly at dimension 2 and the paper\'s quaternion blockwise map at dimension 4. Let `L_j(x)=e_j*x`. For tested Cayley-Dickson basis units, `L_j^2=-I`, so each `B_j=I+L_j` has `B_j^T B_j=2I` and inverse `(I-L_j)/2`. The product `M_D=B_(D/2)...B_4 B_2 B_1` has `M_D^T M_D=D I`, and can be inverted by reversing the stages, dividing by two at each exact stage. The code verifies this integrality throughout the round trip.')
col.append('For D<=64 the entire matrix was independently constructed from basis columns and checked to contain only +/-1 entries and satisfy the Gram identity exactly. Larger dimensions were tested with exact integer roundtrips and norm identities. The construction is a signed Hadamard-type transformation *derived from actual higher-dimensional imaginary-basis multiplication*; it is not a claim of new information or superior asymptotic computation.')
col.append('\n## 2. Exact encoding and dimension scaling\n')
col.append('| Total algebra dimension D | Imaginary axes | Tight one-bit storage | Dense int8 coefficients | Unit multiply time (ms)* |')
col.append('|---:|---:|---:|---:|---:|')
for r in dim_set:
 col.append(f"| {md_num(r['dimension_D'])} | {md_num(r['imaginary_axes'])} | {md_num(r['raw_packed_bytes'])} B | {md_num(r['dense_int8_coefficient_bytes'])} B | {fmt_ms(r['unit_product_ms_median'])} |")
col.append('\n*Median of repeated Numba-compiled CPU runs; timings are prototype diagnostics, not a controlled cross-machine benchmark or a GPU result. The fixed unit is `e_7` where available. All values passed exact bit recovery.\n')
col.append('The million-dimensional 1-bit case represents **1,048,575 binary coefficients**, requiring **131,072 packed bytes**. Materializing one signed 8-bit coefficient per axis costs **1,048,576 bytes**; 64-bit floating coefficients would cost **8,388,608 bytes**. These are storage representations of the *same information*.')
col.append('\n## 3. Increasing the bits stored in each imaginary coefficient\n')
col.append('| Bits per imaginary coefficient | Information capacity in 1,048,575 axes | Tight byte count | Unit multiplication/inversion |')
col.append('|---:|---:|---:|---|')
for r in wide:
 col.append(f"| {r['bits_per_coefficient']} | {md_num(r['input_bits'])} bits | {md_num(r['tight_storage_bytes'])} B | Exact / exact |")
col.append('\nUsing 16-bit coefficients and D=1,048,576, the full successive imaginary mixing transform also inverted **exactly**. Its largest output coefficient in the seeded test had magnitude **129,324,076**, which exceeds the original 16-bit coefficient range: transported arithmetic requires adequate intermediate precision. The exact norm identity was checked with Python arbitrary-size integer summation to avoid overflow of a 64-bit squared-energy accumulator.\n')
col.append('### Direct multiplication versus signed-unit scrambling\n')
col.append(f"- Dense binary-coefficient Cayley-Dickson multiplication was executed through **D=4,096**; the D=4,096 case required {fmt_ms(unit['dense_products'][-1]['cd_product_ms_median'])} ms in this CPU prototype, with {md_num(unit['dense_products'][-1]['pair_multiply_contributions'])} nonzero term products, and yielded {md_num(unit['dense_products'][-1]['nnz_output'])} nonzero output coefficients.\n- True sparse-sparse multiplication was executed up to **D=1,048,576** with 64 nonzero input coordinates on each side: 4,096 term products, {md_num(unit['sparse_products'][-1]['output_nonzero'])} nonzero outputs, {fmt_ms(unit['sparse_products'][-1]['sparse_product_ms_median'])} ms.\n- Fully dense multiplication of two million-dimensional vectors **was not attempted**; its straightforward algorithm would require on the order of a trillion basis term pairs. The million-dimensional *unit multiplication* and *multi-stage mixing* were actually executed and verified.")
col.append('\n## 4. Algebra verification and failure modes\n')
col.append('| Algebra dimension | Property | Result |')
col.append('|---:|---|---|')
col.append('| 2 (complex) | Associative, multiplicative norm | Confirmed |')
col.append('| 4 (quaternions) | Associative, multiplicative norm | Confirmed |')
col.append('| 8 (octonions) | Nonassociative, multiplicative norm | Nonassociativity witnessed by `(e1 e2)e4` versus `e1(e2 e4)` |')
col.append('| 16 (sedenions) | Zero divisors and failed composition law | `(e1+e10)(e4-e15)=0`, both operands nonzero |')
col.append('| 32, 64, 128 | Multiplicative norm generally fails | Counterexamples found |')
col.append('\nRandomized integer coefficient checks of norm multiplicativity:\n')
col.append('| D | Trials | Failures of N(xy)=N(x)N(y) |')
col.append('|---:|---:|---:|')
for r in norm_set:
 col.append(f"| {r['dimension']} | {r['trials']} | {r['norm_multiplicativity_failures']} |")
col.append('\nThese random tests demonstrate counterexamples, not a statistical property of all inputs. The algebraic zero-divisor example is an explicit exact witness, not a probabilistic finding.\n')
col.append('### Associative Clifford alternative\n')
col.append('For `Cl_(0,n)` with n anticommuting imaginary generators, storing an n-bit **vector** is straightforward, but arbitrary multiplication leaves the vector subspace. The full algebra has **2^n** basis coordinates; only the scalar and grade-2 components of a vector product were materialized. At n=1,024 generator bits, there are **523,776 possible bivectors**, of which **205,368** were nonzero in one seeded random product. This avoids the Cayley-Dickson nonassociativity, but is not an ordinary fixed-size n-vector arithmetic format under multiplication.\n')
col.append('## 5. Actual file compression: fair end-to-end comparisons\n')
col.append('The codec applies the algebraically defined exact multi-stage transform to chunks with D-1 imaginary message coefficients, supports **binary (0/1)** and **bipolar (-1/+1)** coefficient models, and chooses raw/fill/sparse-spectral block encodings. Each output includes block tags, dimensional metadata, original length, and CRC32. The completed encoded file is compressed with the smaller of zlib 9 and Brotli 9. Each dataset was exhaustively tested over block dimensions 256, 1,024, 4,096, and 65,536, selecting the smallest result. For conventional comparison, the original bytes were compressed directly with zlib 9, Brotli 9, and LZMA 6. Search over all candidate configurations is not included in individual encoding-time figures.\n')
col.append('| Data | Original bytes | Best direct compressor | Direct size | Best hypercomplex cascade | Result |')
col.append('|---|---:|---|---:|---:|---|')
for ds in ('paper_text','python_source','paper_PDF','document_bitmap','random_control','CD_basis_binary_1024','CD_basis_binary_4096'):
 r=consts[ds]
 label={'paper_text':'Paper extracted text','python_source':'Python source','paper_PDF':'Original paper PDF','document_bitmap':'Document 1-bit bitmap','random_control':'Random bytes','CD_basis_binary_1024':'Synthetic 1,024D basis codewords','CD_basis_binary_4096':'Synthetic 4,096D basis codewords'}[ds]
 sign='smaller' if r['pct_larger_than_baseline']<0 else 'larger'
 col.append(f"| {label} | {md_num(r['input_bytes'])} | {r['best_conventional_name']} | {md_num(r['best_conventional_bytes'])} B | {md_num(r['best_cd_cascade_bytes'])} B | {abs(r['pct_larger_than_baseline']):.1f}% {sign} |")
col.append('\n**Interpretation:** For the five listed real/random inputs, no sparse spectral blocks were selected by the best codec. The special inputs are algorithmically constructed from this *exact same hypercomplex transformation*: one spectral coefficient per block, so the codec dramatically outperforms the general compressors tested. These are in-sample signal-model demonstrations, not independent natural data or evidence of universal compression.\n')
col.append('## 6. One-million-dimensional binary codeword (independent demonstration)\n')
col.append(f"- Generated a signed hypercomplex codeword with D={md_num(million['dimension'])}, selected basis index e_{md_num(million['axis'])}, and fixed real coefficient +1.\n- The source file holds D-1 meaningful imaginary bit coefficients, packed into {md_num(million['sizes_bytes']['original'])} bytes (one trailing padding bit).\n- The transformed binary coefficients yield a single nonzero spectral coefficient; the decoder also handles the single padded bit through a fill block.\n- Encoding, inverse recovery, and CRC32 verification passed; SHA-256 is `{million['sha256']}`.\n")
col.append('| Method | Total bytes |')
col.append('|---|---:|')
for desc,key in [('Original packed source','original'),('zlib level 9','zlib9'),('Brotli quality 9','brotli9'),('LZMA preset 6','lzma6'),('Raw hypercomplex spectral file, including metadata','hypercomplex_raw'),('Hypercomplex + Brotli','hypercomplex_brotli9')]:
 col.append(f"| {desc} | {md_num(million['sizes_bytes'][key])} |")
col.append('\nA tighter specialized **basis-index** encoder can represent a known codeword from a D-sized codebook using only `log2(D)` bits for its index (plus framing and protocol metadata). In this sense the above result is not an information-theoretic breakthrough: its input was explicitly chosen from a low-description-complexity family. The algebra gives a principled way to *recognize and transform* that family, not a way to encode arbitrary random strings in fewer bits than their entropy.\n')
col.append('## 7. Error behavior\n')
col.append('1,000 independent bit-flip trials in D=4,096 signed-unit bit-packed representations showed exactly one recovered data-bit error per one transmitted bit flip, when the flipped coordinate corresponded to a message bit. No native error correction was observed. CRC32 makes corrupted files detectable, but CRC alone cannot repair them. Multi-stage spectral transforms are exactly invertible in noiseless integer arithmetic; neither they nor the algebra change the underlying channel-capacity constraint.\n')
col.append('## 8. What has and has not been established\n')
col.append('**Established:** actual high-dimensional imaginary coefficients can hold and recover independent bits or multi-bit digits; products by basis units and sequential (1+e_j) transformations are implemented with correct signs; the tested mixing map has an exact inverse, a scaled orthogonal Gram identity, and exceptional sparsity on particular algebraically constructed sign patterns; runnable file codecs and corruption tests exist.\n')
col.append('**Not established:** general compression superiority; extra information capacity without paying coefficient precision/storage; GPU acceleration; arithmetic speedup on random large multiplications; safe use of every nonzero high-dimensional key; universal error correction; or a novel entropy lower bound. The dense multiplication test stops at D=4,096, while the million-dimensional tests include real unit products, sparse products, and full sequential mixing—not fully dense binary-by-binary multiplication.\n')
col.append('## 9. Reproduction\n')
col.append('Dependencies: Python 3.10+ (tested with current Python 3.13), `numpy`, `numba`, `brotli`, and `Pillow`. Tests use fixed seeds and do not require network access. From the extracted package folder:\n')
col.append('```bash\npython -m unittest -v hypercomplex_suite_tests\npython hypercomplex_highdim_test.py\npython hypercomplex_cascade_test.py --phase codec\npython hypercomplex_cascade_test.py --phase spectrum\npython hypercomplex_widths_and_errors.py\npython hypercomplex_million_dim_demo.py\npython hypercomplex_combine_results.py\n```\n')
col.append('Try a user file (lossless encode/decode):\n')
col.append('```bash\npython hypercomplex_codec_cli.py encode INPUT.bin OUTPUT.hcb --mode cascade --dimension 4096 --coefficient bipolar --compress brotli\npython hypercomplex_codec_cli.py decode OUTPUT.hcb RESTORED.bin\ncmp INPUT.bin RESTORED.bin\n```\n')
col.append('The results are single-machine CPU/Numba prototype measurements, with short repeated timing samples. The compressor tests involve whole-file oracle selection across known configurations, so comparisons indicate attainable compressed sizes rather than training-free online throughput. No GPU kernel was run.\n')
col.append('## 10. Next experiment suggested by the full-dimensional result\n')
col.append('Use the **bipolar** high-dimensional transform on a measured data source with a hypothesized low-support generator (e.g., actual structured bitmask streams or sparse spectral code families), and evaluate on held-out messages with a full encoder-format overhead and speed comparison against Brotli, LZMA, and a specialized codebook-index compressor. Only after confirming a real data distribution that is sparse in this basis is SIMD/GPU implementation likely to be worthwhile.\n')
(R/'hypercomplex_highdim_full_report.md').write_text('\n'.join(col)+'\n')
print('Wrote full report',len(col),'paragraphs, consolidated results', (R/'hypercomplex_full_results.json').stat().st_size,'bytes')