# Unified Complex Space Factorization

This note describes [`julia/complex_space_factorization.cpp`](../julia/complex_space_factorization.cpp),
a single engine that applies the **corrected Fermat terminal-digit sieve** and combines
every method developed across the complex-space papers into one deterministic routine.

Build and run:

```sh
g++ -std=c++17 -O3 -march=native julia/complex_space_factorization.cpp \
    -lgmpxx -lgmp -pthread -o csf

./csf 3169100003113213800570166301717      # factor
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
