# Complex Trial Multiplication
This method of factorization involves no division, rather tests products `p*q==n` and adjusts `p` or `q` accordingly.

The basis of the method is simply: 

   * if `p*q<n`, raise `q`. 
   * If `p*q>n`, lower `p`.
    
   * All while continuing where the last `p` and the last `q` left off.

The complex Fermat sieves, [described here](../Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf)
define the sequences for `p` and `q` via their complex mapping.

## The two walks, and where they start

The reference implementation runs **two walks from one start point, simultaneously**:

```
ascending :  if p*q > n  ->  TMIS  += 10   else  TMRS  += 10
descending:  if p*q < n  ->  TMdIS -= 10   else  TMdRS -= 10
```

The moves are on `R` and `i` by 10 -- *not* on `p` and `q` independently. Both
branches of the ascending walk raise `q = R + i` by 10, and both branches of the
descending walk lower it by 10. So `q` is monotone in each: **up in one, down in
the other.** That is the ascend and the descend.

**The start is a point on the red/blue boundary, entered from the real axis.**

Fix a real value `R` and climb vertically -- `R` held, `i` rising. The product

```
p*q = (R - i)(R + i) = R^2 - i^2
```

falls monotonically, so for `R >= sqrt(n)` the climb begins red (`i = 0` gives
`R^2 >= n`) and ends blue. The crossing is the seed:

```
i = floor(sqrt(R^2 - n))          the last red point at this R
```

For `n = 309` at `R = 21`: `(21-x)(21+x) = 309` gives `x^2 = 132`, `x = 11.489`,
so `441 - 11^2 = 320` is red and `441 - 12^2 = 297` is blue. The seed is
**`21 + 11i`**. Note what this costs: `R` and `i` *are* the coordinates, so a
seed is one integer square root and **no division**.

### Which real values

```
R in [ ceil(sqrt n), floor(sqrt 2n) ]          1 .. 1.41421 sqrt(n)
```

The interval is forced, not chosen. With `R = (p+q)/2` and `i = (q-p)/2` on the
hyperbola, `q = R + sqrt(R^2 - n)`, so:

| `R / sqrt(n)` | `i / sqrt(n)` | `p / sqrt(n)` | `q / sqrt(n)` |
|---|---|---|---|
| 1.00000 | 0.00000 | 1.00000 | 1.00000 |
| 1.06066 | 0.35355 | 0.70711 | 1.41421 |
| 1.41421 | 0.99999 | 0.41421 | 2.41421 |

That is the whole balanced arc -- both factors, nothing outside it -- and
`sqrt(2n)` is its top **exactly**, with no floating-point constant anywhere.

The real axis and the `q` axis are different axes, and it matters: the RSA
window is stated on `q` as `[1, 1.41421] sqrt(n)`, and the real values carrying
it are `[1, 1.06066] sqrt(n)`.

For `n = 309` the interval is `[18, 24]` -- seven integer real values, seven
seeds:

| R | `R^2 - n` | `i` | `p` | `q` | `p*q` | next `i` |
|---|---|---|---|---|---|---|
| 18 | 15 | 3 | 15 | 21 | 315 red | 308 blue |
| 19 | 52 | 7 | 12 | 26 | 312 red | 297 blue |
| 20 | 91 | 9 | 11 | 29 | 319 red | 300 blue |
| 21 | 132 | 11 | 10 | 32 | 320 red | 297 blue |
| 22 | 175 | 13 | 9 | 35 | 315 red | 288 blue |
| 23 | 220 | 14 | 9 | 37 | 333 red | 304 blue |
| 24 | 267 | 16 | 8 | 40 | 320 red | 287 blue |

### Spacing the seeds

Spreading seeds along that interval is the parallel axis, and each seed runs
**both** walks. But they cannot be spaced evenly in `R`, because the walk out of
a seed is paid in `q` and

```
dq/dR = 1 + R / i          with i = sqrt(R^2 - n)
```

diverges at the balanced corner, where `i -> 0`. At 90 bits, moving the real
part by **one** integer at the bottom of the interval jumps `q` by `5.2e6` --
518,000 walk steps out of that single seed -- while a seed at the top moves `q`
by 2.41 and does nothing. Evenly spaced in `R`, the first seed carries the whole
arc and the rest idle. On a GPU that is the entire launch waiting on lane 0.

So: take every integer real value while they are few (`309` -> 7 seeds, `8051`
-> 37), and past that space them for equal work, which means equal in `q`. Those
are real values on the same interval with the same crossing beneath them --
`R = (q + n/q)/2`, `i = (q - n/q)/2` is that crossing written from the other
side -- sampled densely where the arc turns and sparsely where it runs straight.

### What this replaced

CTM used to seed every walk at the first red point of the 135-degree traversal,
`(r + j_max + 1, m - j_max - 1)` -- `112 + 65i` for `n = 8051` -- where the
strip meets the hyperbola at `q = S ~ 2 sqrt(n)`. Measured, `--only=ctm`,
multiplications to the factor:

| n | bits | from `q = S` | from the real interval | |
|---|---|---|---|---|
| 422491 | 19 | 26,923 | 1 | 13,462x |
| 247401437 | 28 | 37,489 | 4,453 | 8.4x |
| 23550892333 | 35 | 112,619 | 8,195 | 13.7x |
| 14866993139051 | 44 | 3,114,811 | 16,388 | 190x |
| 1183910463999919 | 51 | 27,560,289 | 8,195 | 3,363x |
| 634309085142375407 | 60 | 637,225,959 | 24,582 | 25,921x |
| 120686502968132437433 | 67 | 8,789,279,687 | 1 | 4.4e9x |

Both seedings find the factor in all seven; the gap is distance to it, and it
widens with `n` because a balanced factor sits at the bottom of the interval
while `q = S` starts a fixed fraction of `sqrt(n)` away.

The Fermat bracket and trial division still work the 135-degree traversal, from
`(r, m) = 90 + 87i` toward its crossing. CTM no longer shares that line with
them -- it has its own geometry, and they do not collide.

### Sieving both directions

Every move is `+-10` on `R` or on `i`, so both the increases and the decreases
preserve the terminal digits. A stream synced into an admissible `(R mod 10,
i mod 10)` class at the start stays in it for the whole walk -- which is what
lets the paired Fermat sieve for this `N`'s classification apply to a
two-directional search unchanged.

For `n = 8051`: `n mod 4 = 3` (4k-1), last digit 1, so the paired tables give
`(R, i)` classes `(0,3) (0,7) (4,5) (6,5)`. The factor `83 * 97` is `R = 90`,
`i = 7` -- class `(0,7)`.

The first seed is `R = ceil(sqrt 8051) = 90`, and `90^2 - 8051 = 49`, so
`i = 7` exactly:

```
 90+ 7i  p=83 q= 97   8051 = n     found on the seed itself
```

`(0,7)` is admissible, so the pair needs no rounding at all. The engine tries
`(0,3)` first and spends 4 multiplications in total, against 9 descending from
`112 + 65i`.

Both directions are still needed. A seed is only the entry point: its ascending
walk covers `q` above it and its descending walk covers `q` below. Seed `c`
ascends as far as seed `c+1`'s `q` and descends as far as seed `c-1`'s, with a
pad at each handover to absorb the class rounding.

### Endpoints walk inward only

At the two ends of the interval one direction has nothing to do:

* **`R = ceil(sqrt n)`, the balanced corner: ascend only.** There is nothing
  below it, and that is a theorem rather than a heuristic. Every factor pair has
  `R = (p+q)/2 >= sqrt(pq) = sqrt(n)` by AM-GM, so `ceil(sqrt n)` is the
  smallest real part any pair can have; `q` rises with `R` along the boundary,
  so no factor sits below the first seed's `q` either. Descending off seed 0
  could only re-test an empty strip.
* **The top of the interval: descend only.** Above it is the tail, which the
  `q`-chunked region owns outright -- and under `--rsa` there is no tail at all,
  because the real values stop at `1.06066 sqrt(n)`, the one carrying
  `q = 1.41421 sqrt(n)`, and the window ends there.

The saving is small (1-3% where it shows, since seed 0's descending span is
`floor(sqrt(r^2 - n)) ~ n^0.25` against an arc of `sqrt(n)`), but it is work
that provably cannot find anything.

### What does NOT work: meeting at the midpoint

Seed `c` ascending and seed `c+1` descending cover the same `q` range, so it
looks as though each should stop at the midpoint between them, halving the work
for identical coverage. It measures exactly that way -- **1.98x fewer
multiplications** -- and it is wrong.

The two walks are **not two halves of one sweep**. Ascending from seed `c` and
descending from seed `c+1` trace *different* `(p, q)` trajectories across the
same range, and a trajectory only lands on the true `p` if it is allowed to run
the whole way. Cut at the midpoint, coverage drops from 690/690 to 543/690 --
147 real factorizations lost, every one of them at `q/p` between 3 and 5, well
inside the interval, and every one found by the full-span form.

The apparent redundancy is not redundancy. Each direction runs its full span.

### Measured

`--only=ctm --b1=1 --threads=1`, multiplications to the factor, seeding from
`q = S` against seeding from the real interval:

| n | p * q | q/p | from `q = S` | from the interval | |
| - | ----- | --- | ---- | ---- | --- |
| 143 | 11 * 13 | 1.2 | 9 | **1** | 9x |
| 2923 | 37 * 79 | 2.1 | 245 | **5** | 49x |
| 8051 | 83 * 97 | 1.2 | 690 | **4** | 173x |
| 288419 | 379 * 761 | 2.0 | 23,530 | **690** | 34x |
| 5115191 | 1597 * 3203 | 2.0 | 32,902 | **3,415** | 9.6x |
| 1007509 | 503 * 2003 | 4.0 | 32,769 | **5,684** | 5.8x |
| 17821157 | 3019 * 5903 | 2.0 | 33,023 | **6,060** | 5.4x |
| 10967535067 | 104723 * 104729 | 1.0 | 105,152 | **16,387** | 6.4x |
| 978508015703 | 752867 * 1299709 | 1.7 | 575,818 | **186,741** | 3.1x |
| 10963 | 19 * 577 | 30 | **863** | 1,843 | 0.47x |
| 798607 | 101 * 7907 | 78 | **612** | 9,101 | 0.07x |

The last two rows are the cost of the change. A seed sweep that starts at the
balanced corner and works up reaches `q ~ sqrt(n)` first and `q >> sqrt(n)`
last; `q = S ~ 2 sqrt(n)` started partway along the arc and so was nearer to a
lopsided factor.

Both are **out of scope for `--rsa`**, which assumes `1 < q/p < 2`: `q/p = 78`
and `q/p = 30` cannot occur in that window at all. And in the general engine
small-`p` factorizations belong to trial division and the yellow path, which
reach `101` immediately. CTM's comparative advantage is the balanced region,
where nothing else is cheap, and that is what the interval seeding front-loads.
On balanced semiprimes the advantage grows with `n`:

| n | bits | from `q = S` | from the interval | |
|---|---|---|---|---|
| 422491 | 19 | 26,923 | 1 | 13,462x |
| 14866993139051 | 44 | 3,114,811 | 16,388 | 190x |
| 1183910463999919 | 51 | 27,560,289 | 8,195 | 3,363x |
| 634309085142375407 | 60 | 637,225,959 | 24,582 | 25,921x |
| 120686502968132437433 | 67 | 8,789,279,687 | 1 | 4.4e9x |

End to end, with every stream running, the default path is unchanged -- CTM is
one stream among five and rarely the one that wins on a general `n`. What moved
is what CTM does when it is the stream that matters.

### Two things that must not be added

**No `p < 3` early exit.** `p*q < n` when `p` is tiny, so the rule moves `i` and
raises `p` by 10 again -- the walk recovers. Breaking discards the class first,
which loses `n = 143`, whose seed passes through `p = 1` one step before
`p = 11`.

**Not one walk.** An earlier version of the C++ engine had only the ascending
walk. From the crossing, ascending can never reach a factor whose `q` is below
`S` -- 8051's `q = 97` against `S = 177` -- so it could not factor 8051 at all.
Seeding from the real interval does not retire the descending walk: it is what
covers `q` below each seed, and without it the arc between consecutive seeds is
only half swept.

## Three defects in the routine above

Found by porting the Julia faithfully and testing it; all three are fixed in
[`julia/Complex_Trial_Multiplication.jl`](../julia/Complex_Trial_Multiplication.jl).

1. **`real = r + 1` skips `R = r`.** The factor of a very balanced semiprime
   sits at `R = ceil(sqrt(n))` exactly (the `j = 0` case). Starting one above
   syncs that stream to the next admissible `R` in its digit class, already past
   the answer, so the factor is unreachable. `n = 8051` (83*97, `R = 90`) and
   `n = 143` (11*13) both fail this way.

2. **The loop tested only `TMIS[1]`.** Once stream 1 passed `max_im` every other
   stream was cut off with it. `n = 314187` (3*104729) needs `i = 52363` on
   stream 4, right at `max_im`, and was terminated early.

3. **`n % 10 == 5` throws `UndefVarError`.** No branch assigns the sieve arrays
   for that digit. `5 | n` is structurally different -- 9 admissible `(R,i)`
   classes per parity lane rather than 4 -- see
   [`Complete_Fermat_Sieve_Verified.pdf`](../Number%20Theory%20Papers/Complete_Fermat_Sieve_Verified.pdf).

Over a nine-case suite the published routine failed 2; with the fixes, 0.

The `(R,i)` digit tables themselves are **correct** -- verified class-for-class
against brute force for all eight (lane, digit) combinations. They also keep the
pairing, which the C++ engine discards in favour of the b-digit marginal.

## Examples
```julia
@benchmark CTM(798607)
BenchmarkTools.Trial: 
  memory estimate:  18.16 KiB
  allocs estimate:  27
  --------------
  minimum time:     23.606 μs (0.00% GC)
  median time:      25.146 μs (0.00% GC)
  mean time:        41.549 μs (10.93% GC)
  maximum time:     10.666 ms (97.68% GC)
  --------------
  samples:          10000
  evals/sample:     1
  time tolerance:   5.00%
  memory tolerance: 1.00%
```
```julia
@benchmark CTM(978508015703)
BenchmarkTools.Trial: 
  memory estimate:  18.14 KiB
  allocs estimate:  27
  --------------
  minimum time:     461.867 μs (0.00% GC)
  median time:      539.359 μs (0.00% GC)
  mean time:        696.977 μs (0.69% GC)
  maximum time:     10.191 ms (92.18% GC)
  --------------
  samples:          6739
  evals/sample:     1
  time tolerance:   5.00%
  memory tolerance: 1.00%
```

## Julia code
Below is the `julia` code.  There are 4 paired multiplications per iteration.  One immediate efficiency would be from the parallelization of these multiplications within the `while` loop.

``` julia
function CTM(n)
    r = Newton_sqrt(n)
    max_im= div((n-9),6)
    real = r+1

# FERMAT SIEVES
        last_digit = n % 10

        if((n+1)%4==0)
          if(last_digit==1)
            real_sieve=[0,4,6,0,10]
            imaginary_sieve=[3,5,5,7]
          end
          if(last_digit==3)
            real_sieve=[2,8,2,8]
            imaginary_sieve=[1,9,9,1]
          end
          if(last_digit==7)
            real_sieve=[4,6,4,6]
            imaginary_sieve=[3,7,7,3]
          end
          if(last_digit==9)
            real_sieve=[0,2,8,0,10]
            imaginary_sieve=[1,5,5,9]
          end
        else
          if(last_digit==1)
            real_sieve=[1,5,5,9]
            imaginary_sieve=[0,2,8,0]
          end
          if(last_digit==3)
            real_sieve=[3,7,3,7]
            imaginary_sieve=[4,6,6,4]
          end
          if(last_digit==7)
            real_sieve=[1,9,1,9]
            imaginary_sieve=[2,8,8,2]
          end
          if(last_digit==9)
            real_sieve=[3,5,5,7]
            imaginary_sieve=[0,4,6,0]
          end
    end #FERMAT SIEVE

# SYNC REAL TO REAL SEQUENCE
    last_digit_real = real % 10
    if any(last_digit_real!=real_sieve)
      real_init_diff = real_sieve - last_digit_real
      if any(real_init_diff.>=0)
      real_init_diff = real_init_diff[real_init_diff.>=0][1]
      else
      real_init_diff = real_init_diff[end]
      end
      real = real + real_init_diff
    end

# ELIMINATE EXTRA SYNCHING SIEVE ENTRIES FOR POSITIONAL MATCHING MULTIPLICATION (IF N ENDS IN 1, TEST ONLY PRODUCTS POSSIBLY ENDING IN 1)
real_sieve = real_sieve[real_sieve.<=9]

# ALIGN ALL REALS TO SYNCHED POSITIONS
TMRS = (real+(real_sieve - real%10))
lTMRS = length(TMRS)

# MAKE SURE 0 IMAGINARY STARTS AT 10
imaginary_sieve[imaginary_sieve.==0] = 10
TMIS = imaginary_sieve

# SETS LOWER BOUNDARIES FOR (p) AND UPPER BOUNDARIES FOR (q) FOR EACH SIEVE ENTRY 
ceilings = div(n,3)
ceiling_p = [ceilings,ceilings,ceilings,ceilings]
floor_q = [3,3,3,3]


while (TMIS[1] <= max_im)
  for i in 1:lTMRS
    p = TMRS[i] - TMIS[i]
    if(p>ceiling_p[i])
      TMIS[i] = TMIS[i] + 10
    end
    if(p<=1)
        TMRS[i] = TMRS[i] + 10
        p = TMRS[i] - TMIS[i]
    end
    q = TMRS[i] + TMIS[i]
    if(q<floor_q[i])
      TMRS[i]=TMRS[i] + 10
    end
    N = p*q

    if (N == n) return(p, q) end

    if (N > n)
      ceiling_p[i]=p
      TMIS[i] = TMIS[i] + 10
    else
      floor_q[i]=q
      TMRS[i] = TMRS[i] + 10
    end

  end #for

end # while

end
```
