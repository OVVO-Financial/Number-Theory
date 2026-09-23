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

Let `j = a - ceil(sqrt(N))` be the vertical Fermat depth. The Sec. 5 identity does not
hold in general, and the reason is sharper than it first looks.

**The construction collapses.** With `IA_k = N - (N - a0)/2^k`, the target is

```
T_k = IA_k - b/2^k = N - (N - a0 + b) / 2^k
```

`T_k` is only defined when that division is exact, i.e. `N - a0 + b = 2^k * X`. But `N` is
odd, so `gcd(2^k, N) = 1` and dividing by `2^k` cannot change a gcd with `N`:

```
gcd(T_k, N) = gcd(N - a0 + b, N) = gcd(a0 - b, N)      for EVERY k
```

Every level of the iterated average returns the identical gcd. Verified:

| N | `j` | `gcd(a0-b,N)` | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 |
| - | --- | ------------- | --- | --- | --- | --- | --- | --- |
| 8051 | 0 | 83 | 83 | 83 | 83 | 83 | 83 | — |
| 143 | 0 | 11 | 11 | 11 | — | — | — | — |
| 10967535067 | 0 | 104723 | 104723 | 104723 | 104723 | — | — | — |
| 1053018024371 | 523 | 1 | — | — | — | — | — | — |

The `—` entries are where `T_k` stops being an integer, which happens at exactly
`k > v2(q-1)`: for `N = 8051`, `q - 1 = 96 = 2^5 * 3`, so `k <= 5`. The `k` iterations are
not `k` independent chances at a factor — they are one arithmetic fact re-tested `k` times,
until the halving exhausts the powers of 2 in `q - 1`.

**What the test really is.** Substituting `a0 = a - j`:

```
a0 - b = (a - b) - j = p - j          a0 + b = (a + b) - j = q - j
```

so the descending branch asks whether `p - j` divides `N`, and the ascending branch asks
the same of `q - j`. Unrolled through `N - p = p(q-1)`:

```
gcd(T_k, N) = gcd( p(q-1) + j , N ),   and since gcd(2^k, p) = 1:

    p | T_k  <=>  p | j              q | T_k  <=>  q | (p - j)
```

With `0 <= j < p` the first forces `j = 0` and the second forces `j = p`, a contradiction.
The side condition `j < p` holds unless `q/p >= 3 + 2*sqrt(2) = 5.828` — checked on 3,000
semiprimes, with no exceptions — and that is also the regime where trial division finds
`p` at once.

Measured over 300 random semiprimes:

| | identity yields a factor |
| - | - |
| `j == 0` | **51 / 51** |
| `j > 0` | **0 / 249** |

**Why `j = 0` is the wrong thing to be good at.** Since
`a - sqrt(N) = (sqrt(q) - sqrt(p))^2 / 2`, the condition `j = 0` means
`(sqrt(q) - sqrt(p))^2 < 2`, i.e.

```
b  <~  sqrt(2) * N^(1/4)
```

(tested as a predictor on 2,000 semiprimes: zero mismatches). But `j = 0` is the definition
of "Fermat succeeds on its first trial" — `a = ceil(sqrt(N))` is candidate number one. The
identity fires exactly when a single `isqrt` would already have finished, and is silent
everywhere else.

**The Sec. 5 open question.** The paper asks why the partial products "are residing near or
on iterated averages". They are not, and it has nothing to do with the imaginary
coefficient. When `j = 0` the quantity tested is `N - p = p(q-1)`, a multiple of `p` by
construction; `gcd(N - p, N) = p` is an identity, not a discovery. Halving keeps it a
multiple of `p` while `2^k | (q-1)`, and the halving is invisible to the gcd because `N` is
odd. The apparent alignment between `b/2^k` and `IA_k` is the same `2^k` cancelling against
itself on both sides.

**What Sec. 5.2 gets right.** The remainder sieve is correct and is the best part of the
section. `T_k` is an integer iff `b = Num_k (mod 2^k)`; intersecting with `b = d (mod 10)`
gives `b = r (mod 5*2^k)` by CRT, so `b` advances by `5*2^k` while `T` advances by `5` —
exactly the "every 5th integer" claim, with `5*2^3 = 40` and `5*2^4 = 80` matching the
worked example, and `2^(k-1)` admissible `b` per series. The engine implements it.

What it accelerates is the problem. Since the identity is inert for `j > 0`, the Sec. 7.3
routine is a blind search for any multiple of `p` near `IA_k`, costing about `p/2` gcds,
which Sec. 5.2 cuts by 5. Trial division does the same `O(p)` work with a division per
step rather than a gcd.

**Why the seed-sweeping implementations still factor things.** They sweep
`a_s = a0 + 2s`; when a seed lands on the true `a` that seed has `j = 0` and the identity
fires. The seed-sweep version is therefore vertical Fermat in disguise — the same `a`
ladder, with each rung verified by a chain of gcds instead of one perfect-square test.
That is the gap measured in 3.5.

The `ia` stream is consequently **off by default** and is not exhaustive; running
`--only=ia` on a composite it cannot split reports `INCOMPLETE` on stderr and exits 2
rather than silently presenting `N` as its own factor.

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

### 3.6 The yellow path (Sec. 4, Figure 8) as a fused stream

The "key complex number" `(r + (r-3)i)`, `r = ceil(sqrt(N))`, and the 135-degree line
through it give a single index that drives all three methods at once. With `m = r-3`
and `S = r+m`:

```
R_j = r + j      i_j = m - j      P_j = R - i = 2j + 3      R + i = S  (constant)
```

Three facts make the step very cheap (all verified in `--selftest`, 0 violations over
18,000 `(N,j)` pairs):

1. `V_j = R_j^2 - N` and `H_j = i_j^2 + N` advance by **addition alone**:
   `V += 2R + 1`, `H -= 2i - 1`.
2. They are **linearly linked**: `H_j = V_j - P_j*S + 2N`, so one state variable
   carries both Fermat tests.
3. Wheel residues advance by increment/decrement.

So the inner loop is a few additions, two table lookups and one small division — no
multiplication, no square root in the common case.

**The path is exactly `sqrt(N)/4` steps and provably complete.** It stops when
`P_j * S > N`, i.e. `j_max = (N/S - 3)/2`; the measured ratio `j_max / sqrt(N)`
converges to `0.2500`. Coverage: the trial-division leg catches any `p <= 2 j_max + 3
~ sqrt(N)/2`, and for `p > sqrt(N)/2` one has
`j_true = (sqrt(q) - sqrt(p))^2 / 2 <= sqrt(N)/4`, so the vertical leg catches it. The
two legs cover each other exactly — which is what Figure 8 shows geometrically.
Verified on 400 semiprimes: 0 misses. For `N = 309` the whole path is `j = 0,1,2,3`.

### The cost is known before you start

Because `Q = R + i = S` is invariant along the path, the strip's length follows from **one
division**. Worked on `N = 8051`, `r = ceil(sqrt 8051) = 90`, `m = 87`, `S = 177`:

| point | `P = R - i` | `Q = R + i` | `P*Q` | |
| --- | --- | --- | --- | --- |
| `90 + 87i` | 3 | 177 | 531 | below N |
| `100 + 77i` | 23 | 177 | 4,071 | below N |
| `110 + 67i` | 43 | 177 | 7,611 | below N |
| `112 + 65i` | 47 | 177 | 8,319 | **above N — red** |

`Q` never moves. So the path ends at the largest odd `P` with `P*S <= N`, i.e.
`P <= N/S = 8051/177 = 45.48`, giving `P = 45`, and `P = 2j+3` puts that at `j = 21` —
a strip of **22**, `j = 0..21`.

That is the practical point: `N/S` is a single division, and it tells you *exactly* how
many steps the deterministic factorization will take before a single step is walked. Most
factoring methods cannot say what they will cost until they are done.

(Starting from `floor(sqrt N)` instead indexes the same strip from one point lower and is
also complete — checked over 20,000 semiprimes, 0 uncovered either way — but `ceil` is the
convention here, and the two agree on the length to within 1.)

### How the strip grows

`Q_j = R_j + i_j = r + m = S` is **constant** along the 135° line — that is what makes it
a line — while `P_j = 2j + 3` climbs. So the path meets the red boundary, the hyperbola
`P·Q = N`, at exactly

```
j_max = floor( (floor(N/S) - 3) / 2 ),      S = 2r - 3,   r = ceil(sqrt N)
```

Writing `r = sqrt(N) + e` with `e = ceil(sqrt N) - sqrt(N)` in `[0,1)` and expanding:

```
j_max = sqrt(N)/4 - 9/8 - e/4 + O(1/sqrt N)
```

so the strip is **`sqrt(N)/4` long — `2^(b/2 - 2)` for a b-bit N**. For `N = 309` that is
`j = 0,1,2,3`: a strip of 4, which is the whole path.

Floating point cannot check this past ~96 bits: at 128 bits `sqrt(N)/4` is near `2^62` and
a double resolves it only to about 512, while the quantity being bounded is `O(1)`.
Multiplying through by 4 makes it exact integer arithmetic at any size:

```
4 * j_max = isqrt(N) - c
```

with `c` in `[4, 9]` over 3,006 values of N from 12 to 512 bits — a bounded constant, the
width coming from the two floors. `--selftest` check 6 asserts this.

| bits | strip length | 1 thread | 4 threads |
| --- | --- | --- | --- |
| 9 (`N=309`) | 4 | 56 ns | 16 ns |
| 32 | 1.6e4 | 162 µs | 46 µs |
| 48 | 4.2e6 | 41 ms | 12 ms |
| 64 | 1.1e9 | 10.6 s | 3.0 s |
| 80 | 2.7e11 | 45 min | 13 min |
| 96 | 7.0e13 | 8.1 d | 2.3 d |
| 112 | 1.8e16 | 5.7 yr | 1.6 yr |
| 2048 | 1e308 | — | — |

at 9.89 ns/step marginal (setup cancelled between two runs of very different length).

**Why the constant is 1/4, and not something arbitrary.** Write the smaller factor as
`p = sqrt(N)/t`, `t >= 1`. The trial-division leg reaches `2*j_max + 3 = N/S -> sqrt(N)/2`,
so it covers `t >= 2` outright. For `1 <= t < 2` the vertical leg needs depth

```
j_true = (p+q)/2 - ceil(sqrt N)  ~  sqrt(N) * (t-1)^2 / (2t)
```

which increases on `[1,2]` and equals `sqrt(N)/4` exactly at `t = 2`. The two legs meet at
the worst case and neither is longer than it has to be — that is the red boundary. Measured
over 8,000 semiprimes: 0 uncovered, and among those beyond the trial leg's reach the worst
`j_true/sqrt(N)` is **0.2495** against the bound 0.2500, the tightest case landing 126
steps inside a 61,512-step path.

The engine's own step counts confirm the model end to end — `--only=yp --threads=1` walks
exactly `min(j_true, (p-3)/2) + 1` steps, ratio 1.000 on all five cases checked.

**The strip is not one step longer than it needs to be.** There exist N where
`j_true == j_max` exactly, so the factor sits on the very last step:

| N | | strip | steps walked |
| --- | --- | --- | --- |
| 1,007,509 | 503 × 2,003 | 250 | 250 |
| 1,289,923 | 569 × 2,267 | 283 | 283 |
| 17,005,309 | 2,063 × 8,243 | 1,030 | 1,030 |

`--selftest` pins these, so shortening the strip by even one step fails the suite. Finding
them takes some care: the tight corner is *not* at `t = 2`, since there `p = sqrt(N)/2` is
exactly the trial leg's reach and the vertical leg is never asked. The binding cases sit
just *past* the reach, and they are rare enough that a random sweep over a wide range found
none in 400,000 tries — these came from an exhaustive scan of `p` in `[500, 40000)`.

On by default, ablated with `--no-yp`, isolated with `--only=yp`, and given the whole
thread pool with `--parallel`. `8051` resolves in 1 step; `798607` in 50
(`p = 101`, so `2j+3 = 101` at `j = 49`).

Because every quantity along the path is a closed form in `j`, a bracket can begin at
any `j` for O(1) cost — there is no dependency chain to preserve. `--parallel` uses
that: it hands the path the whole pool and stands the other streams down, which is
sound because the two legs cover each other, so the path is complete on its own.
Measured min-of-7, `--b1=1`:

| case | 1 thread | 2 | 3 | 4 |
| --- | --- | --- | --- | --- |
| 62-bit | 11 ms | 9 ms | 8 ms | 8 ms |
| 96-bit | 28 ms | 23 ms | 21 ms | 19 ms |
| 112-bit | 117 ms | 62 ms | 43 ms | **33 ms** (3.5×) |

Shallow cases are dominated by fixed setup, so the deeper the scan the closer it gets
to linear.

**But fusing loses to racing, and the reason is the sieve.** Locking `a = r + j` to
`p = 2j + 3` means the wheel can suppress the `isqrt` but cannot *skip the iteration*:

```
yellow path cost = min( j_true , p/2 )            raw iterations
race cost        = min( p/2 , j_true / ~4000 )    the vertical scan skips
```

so the race is never worse. Measured (4 threads, `--b1=1000`):

| regime | yellow path | vertical only | race |
| ------ | ----------- | ------------- | ---- |
| balanced `p ~ N^0.50` | 11 ms | 11 ms | **10 ms** |
| mid `p ~ N^0.42` | 57 ms | 212 ms | **10 ms** |
| mid `p ~ N^0.38` | 57 ms | 1,286 ms | **8 ms** |
| unbalanced `p ~ N^0.30` | 12 ms | 1,489 ms | **7 ms** |

The fused path does beat the *vertical stream alone* by 20-100x in the unbalanced
regimes — its trial-division leg gets there first. That is exactly the complementarity
Figure 8 describes. It just does not beat running the same three legs at independent
rates.

### 3.7 A bug this exercise found

Before this comparison the trial-division stream was capped at `N^(1/3)` whenever the
Lehman stream was enabled, on the reasoning that Lehman's theorem only requires every
prime below `N^(1/3)` to be covered. That is sound asymptotically and wrong in practice:
Lehman's per-candidate constant is far worse than one division, so a `p` modestly above
`N^(1/3)` was handed to the slow stream. At `p ~ N^0.38` the full race *lost to the
fused path*, 127 ms against 61 ms.

Trial division now always runs to `sqrt(N)`. It costs nothing asymptotically — Lehman
still bounds the worst case at `O(N^(1/3))` and the two race — and the measured effect
is large:

| regime | capped at `N^(1/3)` | uncapped |
| ------ | ------------------- | -------- |
| `p ~ N^0.42` | 20 ms | **10 ms** |
| `p ~ N^0.38` | 127 ms | **8 ms** |

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
