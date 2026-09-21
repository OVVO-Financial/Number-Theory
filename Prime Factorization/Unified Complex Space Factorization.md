# Unified Complex Space Factorization

This note describes [`julia/complex_space_factorization.cpp`](../julia/complex_space_factorization.cpp),
a single engine that applies the **corrected Fermat terminal-digit sieve** and combines
every method developed across the complex-space papers into one deterministic routine.

Build and run:

```sh
g++ -std=c++17 -O3 -march=native julia/complex_space_factorization.cpp \
    -lgmpxx -lgmp -pthread -o csf

./csf 3168987219233877513774136225800517      # factor
./csf --selftest                            # verify every claim below
```

---

## 1. The sieve correction, and what it actually is

### 1.1 Setup

For odd `N = p*q` with `p <= q`, write

```
R = (p+q)/2 ,   i = (q-p)/2 ,   N = R^2 - i^2 = (R-i)(R+i)
```

`R` and `i` always have opposite parity, and which one is odd is fixed by `N mod 4`:

| lane | condition | `R` | `i` |
| ---- | --------- | --- | --- |
| A | `N = 1 (mod 4)` | odd | even |
| B | `N = 3 (mod 4)` | even | odd |

### 1.2 What was corrected

[`Complete_Fermat_Sieve_Verified.pdf`](../Number%20Theory%20Papers/Complete_Fermat_Sieve_Verified.pdf)
corrects a table indexed by `N mod 10` alone. Because `N mod 10` does **not** determine
`N mod 4` (e.g. `11` and `21` are both `1 mod 10`, but `11 = 3 mod 4` and `21 = 1 mod 4`),
such a table must draw from **both** lanes: **8 admissible `(R,i) mod 10` classes per
digit, not 4.** A 4-class version has silently dropped one lane, and any bound derived
from it is invalid.

The `--selftest` reproduces the published tables by direct enumeration of
`N = (R-i)(R+i)` and confirms them exactly:

```
N = 1,3,7,9 (mod 10)  ->  8 classes   (4 per lane)
N = 5 (mod 10)        -> 18 classes   (9 per lane)
```

In an implementation we always know `N mod 4`, so we select the lane directly and use
its 4 classes (9 when `5 | N`). The engine's `fermat_digit_pairs()` is verified
class-for-class against brute force for **all ten odd residues mod 20** — including
the `5 | N` case that the 8-class table excludes.

### 1.3 The structural result: the digit sieve *is* a quadratic-residue sieve

The terminal-digit sieve asks which `(R mod 10, i mod 10)` classes can occur. Since
`i^2 = R^2 - N`, that is nothing other than the statement

> `R^2 - N` must be a square modulo 10,

and the lane condition is the same statement modulo 4. Putting them together:

```
        the corrected terminal-digit sieve
                     ==
      the quadratic-residue sieve at modulus 20
```

`--selftest` section 2 verifies this equality of admissible sets for all ten odd
residues `mod 20`. It is not an analogy — the two sieves admit exactly the same `R`.

### 1.4 Why that matters

Once the sieve is recognised as "`a^2 - N` is a square mod `m`", the modulus is a free
parameter. The engine therefore replaces the single modulus 20 with a CRT wheel:

```
primary   2^6 * 3^3 * 5^2 * 7 * 11 * 13 * 17   = 735,134,400
secondary 19, 23, 29, 31, 37, 41, 43, 47       (tested per surviving candidate)
```

Measured fraction of `a` surviving (from `--selftest` section 3):

| sieve | modulus | density | reduction |
| ----- | ------- | ------- | --------- |
| corrected terminal-digit sieve | 20 | 0.2 | **5x** |
| full wheel + secondary | 735,134,400 (+8) | 2.38e-06 | **419,644x** |

The admissible-residue list is built by CRT composition, so construction costs
`O(number of admissible residues)` rather than `O(modulus)`.

---

## 2. The combined method

Section 4 of *Complex Space Factorization* observes that trial division, vertical
Fermat and horizontal Fermat are each most efficient on a different part of the factor
strip, and proposes running them simultaneously. In this engine that becomes literal:
**all of the Fermat-type streams are one routine** pointed at different `M`.

`fermat_scan_*` finds `x` with `x^2 - M` a perfect square, then offers
`gcd(|x-y|, N)` and `gcd(x+y, N)`:

| stream | `M` | `x` | covers |
| ------ | --- | --- | ------ |
| vertical Fermat (`vf`) | `N` | `a` | bottom of the strip, from `ceil(sqrt(N))` up |
| horizontal Fermat (`hf`) | `-N` | `b` | middle of the strip, from `b = 1` up |
| Lehman multiplier (`lm`) | `4kN` | `a` | everything else; the completeness guarantee |
| trial division (`td`) | — | `p` | top of the strip |
| iterated average (`ia`) | — | — | Sec. 5 GCD accelerator |

They race across the thread pool and the first split wins.

### 2.1 The Lehman stream

Multiplying `N` by `k` rescales the complex space so that a different part of the
factor strip is brought down into the near-vertical region where Fermat's method is
efficient. If `N` has no prime factor below `N^(1/3)`, then for some `k <= N^(1/3)`

```
ceil(sqrt(4kN)) <= a <= ceil(sqrt(4kN)) + N^(1/6) / (4 sqrt(k))
```

contains a solution of `a^2 - b^2 = 4kN`. Combined with the trial-division stream
covering every prime below `N^(1/3)`, the engine is **deterministic and `O(N^(1/3))`**.

This matters: racing trial division against vertical Fermat *alone* gives only
`O(N^(3/8))`. Writing `p = sqrt(N)/c`, `q = c*sqrt(N)`, the trial-division cost is
`sqrt(N)/c` and the vertical-Fermat depth is `N^(1/4)(sqrt(c) - 1/sqrt(c))^2 / 2`;
these balance at `c ~ N^(1/8)`, giving `N^(3/8)`. The multiplier sweep removes that
worst case.

### 2.2 The iterated-average stream (Sec. 5 / 5.2)

With `a0 = ceil(sqrt(N))`, the iterated averages are

```
IA_k = (IA_{k-1} + N)/2 = (N*(2^k - 1) + a0) / 2^k  =  Num_k / 2^k
```

The paper's test is `GCD(IA_k - b/2^k, N)`. Writing `T = (Num_k - b)/2^k`, the target
`T` is an integer exactly when `b = Num_k (mod 2^k)`. Intersecting that with the
terminal-digit condition `b = db (mod 10)` gives, by CRT,

```
b = r (mod 5 * 2^k)
```

and since `b` advances by `5*2^k` while `T` advances by `5`, the scan visits **every
5th integer** below `IA_k` — which is exactly the Sec. 5.2 claim, and the reason the
count of admissible `b` per series is `2^(k-1)`. The engine implements both the
descending and ascending branches over all `k` and all admissible terminal digits.

The sieve is correct. The underlying identity, however, turns out to hold only when
`ceil(sqrt(N))` already equals `a = (p+q)/2` — see **3.3**, which derives the condition
and measures it. That is why this stream is off by default.

---

## 3. Measurements

All timings: 4 threads, `--b1=10000`, on the same machine. `j = a - ceil(sqrt(N))` is
the vertical-Fermat scan depth — the honest difficulty measure for a Fermat-type method.

### 3.1 Sieve ablation

Identical engine, identical streams — only the sieve changes. This isolates exactly
what the generalisation of the terminal-digit sieve is worth.

| N | `j` | no sieve | digit sieve only (`mod 20`) | full wheel |
| - | --- | -------- | --------------------------- | ---------- |
| 96-bit  | 1.0e6 | 63 ms | 8 ms | 15 ms |
| 96-bit  | 1.0e7 | 603 ms | 111 ms | **9 ms** |
| 112-bit | 1.0e7 | 620 ms | 75 ms | **38 ms** |
| 112-bit | 1.0e8 | 5,934 ms | 1,080 ms | **26 ms** |
| 128-bit | 1.0e8 | 6,254 ms | 764 ms | **149 ms** |
| 127-bit | 3.0e8 | 19,024 ms | 3,452 ms | **48 ms** |

At the deepest case the corrected digit sieve alone is worth **5.5x**, and generalising
it to the wheel is worth a further **72x** — **396x** end to end. Below about `j = 1e6`
the wheel's construction cost dominates and the plain digit sieve can be ahead; the
engine sizes the wheel by the bit-length of `N` for that reason.

Run-to-run variation at a fixed `j` (e.g. 149 ms vs 48 ms above) comes from where the
true `a` falls relative to chunk boundaries across the four threads, not from the sieve.

### 3.2 The five streams, in isolation

`--only=<stream>` runs a single stream. `j = a - ceil(sqrt(N))` is the vertical scan depth.

| N | `j` | vertical Fermat | horizontal Fermat | iterated average | trial division | Lehman |
| - | --- | --------------- | ----------------- | ---------------- | -------------- | ------ |
| 8051 (paper) | 0 | 0 ms | 0 ms | 0 ms | 0 ms | 0 ms |
| 40-bit | 523 | 0 ms | 0 ms | *not found* | 3 ms | 42 ms |
| 56-bit | 9,881 | **4 ms** | 6 ms | *not found* | 309 ms | 1,748 ms |
| 96-bit | 1.0e6 | **14 ms** | 329 ms | *not found* | >120 s | >120 s |
| 112-bit | 1.0e7 | **43 ms** | 68,985 ms | *not found* | >120 s | >120 s |

Two things to read off this table.

**Vertical Fermat dominates horizontal Fermat, quadratically.** Writing
`p, q = sqrt(N)(1 -+ e)`, the horizontal scan length is `b ~ e*sqrt(N)` while the
vertical scan depth is

```
j = a - ceil(sqrt(N)) = (sqrt(q) - sqrt(p))^2 / 2  ~  b^2 / (2 sqrt(N))
```

so `j / b ~ b / (2 sqrt(N))`, which is far below 1 in exactly the balanced regime the
complex-space picture targets. At 112 bits, `b = 1.06e12` against `j = 1.0e7` — and the
measured gap is 68,985 ms against 43 ms. Sec. 4 of the paper is right that horizontal
Fermat covers a different part of the strip, but as a *search order* it is strictly worse
than the vertical scan, which is why it is off by default.

**Trial division and Lehman are the other end of the strip.** They are the wrong tool for
balanced `N` and the only tool for unbalanced `N`. That is precisely why the engine races
all of them rather than choosing.

### 3.3 The iterated-average method requires `j = 0`

The Sec. 5 identity does not hold in general. Substituting `a0 = ceil(sqrt(N)) = a - j`
into the descending target:

```
T   = IA_k - b/2^k = N - (N - a0 + b)/2^k
N - a0 + b = N - (a - j) + b = N - (a - b) + j = N - p + j

so  gcd(T, N) = gcd( (p(q-1) + j) / 2^k , N )
```

and since `2^k` is coprime to `p`,

```
p | T   <=>   p | j
q | T   <=>   q | (p - j)
```

With `0 <= j < p` — which holds whenever the semiprime is not wildly unbalanced — both
conditions collapse to **`j = 0`**: the iterated-average GCD returns a factor exactly when
`ceil(sqrt(N))` already equals `a = (p+q)/2`.

Measured over 300 random semiprimes:

| | identity yields a factor |
| - | - |
| `j == 0` | **51 / 51** |
| `j > 0` | **0 / 249** |

The paper's worked example `N = 8051` has `a = 90 = ceil(sqrt(8051))`, i.e. `j = 0`, which
is why it works there and in every re-iteration of the average. This also answers the open
question in Sec. 5 — *"why these partial products ... are residing near or on iterated
averages ... is still a mystery"*. It is not a property of the iterated averages at all:
when `j = 0` the quantity being tested is `N - p`, and `gcd(N - p, N) = p` identically.
Halving it `k` times keeps it a multiple of `p` for as long as `2^k` divides `q - 1`.

The `j = 0` case is precisely the case where plain Fermat succeeds on its very first
trial, so the identity carries no search advantage over the vertical scan.

The Sec. 5.2 remainder sieve is a separate and correct observation: the CRT of
`b = Num_k (mod 2^k)` with `b = db (mod 10)` really does leave every 5th integer below
`IA_k`, and the engine implements it. But what it accelerates is a blind GCD scan, not a
structurally guided one. The `ia` stream is therefore **off by default** and is not
exhaustive; running `--only=ia` on a composite it cannot split reports
`INCOMPLETE` on stderr and exits 2 rather than silently presenting `N` as its own factor.

### 3.4 Lehman window scan

The multiplier windows are short — `N^(1/6) / (4 sqrt(k))` — so materialising a
residue-list wheel per `k` costs far more to build, and to skip past, than the window
itself contains. The stream instead applies the same quadratic-residue condition
incrementally, one modulus at a time, straight across the window:

| N | wheel rebuilt per `k` | incremental small-window sieve |
| - | --------------------- | ------------------------------ |
| 40-bit | 42 ms | **10 ms** |
| 56-bit | 1,748 ms | **19 ms** |

A 92x improvement on the second case, from choosing the sieve representation to match
the length of the interval being scanned.

### 3.5 Against the existing routine

`paper_test_gmp_parallel_scalable_chunks` versus `csf`, both on **default settings**,
4 threads, 60 s cap. Both solve the same problem and both are correct where they finish
(the 40-bit case returns `993913 * 1059467` from either program).

| case | `j` | existing | csf |
| ---- | --- | -------- | --- |
| 40-bit | 523 | 497 ms | **6 ms** |
| 56-bit | 9,881 | >60 s | **9 ms** |
| 64-bit | 37,682 | >60 s | **8 ms** |
| 96-bit | 1.0e6 | >60 s | **21 ms** |
| 112-bit | 1.0e7 | >60 s | **45 ms** |
| 128-bit | 1.0e8 | >60 s | **154 ms** |

This should be read fairly. The existing program is a Sec. 5.2 accelerator with roughly
twenty tunable parameters, and it was not tuned here; a better parameter set would move
these numbers. But its core search is the iterated-average GCD scan, and §3.3 shows why
that scan has no structural edge once `j > 0` — which is every row below the first.

---

## 4. Verification

`./csf --selftest` runs five checks, all of which must pass:

1. **Corrected sieve vs brute force.** `fermat_digit_pairs()` is compared class-for-class
   against direct enumeration of `N = (R-i)(R+i)`, for all ten odd residues `mod 20`.
   This includes `N = 5 (mod 10)`, where the answer is 9 classes per lane, not 0.
2. **Digit sieve = QR sieve at modulus 20.** The two admissible sets are computed
   independently and compared.
3. **Sieve strength.** Reports the measured density of both sieves.
4. **Wheel soundness.** For 400 random semiprimes, the genuine `a = (p+q)/2` must survive
   the wheel. A sieve that ever dropped it would make the engine silently incomplete.
   Measured: 0 drops in 400.
5. **Factorization correctness.** A suite spanning both parity lanes, every `N mod 10`
   including `5 | N`, perfect squares, perfect powers, primes, highly composite numbers,
   strong pseudoprimes, and balanced semiprimes up to 112 bits. Products are re-multiplied
   and every returned factor is primality-tested.

Additionally, every split found by a Fermat-type stream is checked against the squared
complex space identity `sqrt(N^2 + (2ab)^2) = a^2 + b^2` before being accepted.

---

## 5. What this does and does not claim

**What it is.** The fastest deterministic realisation of the complex-space program in
this repository: a strip traversal whose sieve is the corrected terminal-digit sieve
taken to its natural conclusion, with a worst case of `O(N^(1/3))`.

**What it is not.** `O(N^(1/3))` is exponential in the *length* of `N`. For a balanced
semiprime of any cryptographic size this engine — like every Fermat-type method — is not
competitive with the quadratic sieve or the number field sieve, which are subexponential.
The wheel's ~420,000x is a constant factor, not a change of complexity class.

The regime where this approach genuinely wins is the one the complex-space picture was
built for: `N` whose factors are close to `sqrt(N)`, where the vertical scan depth
`j = a - ceil(sqrt(N))` is small. There it is fast, deterministic, needs no randomness,
and produces an exactly certifiable answer.

**On the deterministic claim.** Sec. 7 of *Complex Space Factorization* lists determinism
as a benefit over Pollard's rho. That is preserved here: no stream uses randomness, and
the `td` + `lm` pair makes termination provable rather than heuristic.
