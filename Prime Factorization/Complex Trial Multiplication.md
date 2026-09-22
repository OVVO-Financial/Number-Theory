# Complex Trial Multiplication
This method of factorization involves no division, rather tests products `p*q==n` and adjusts `p` or `q` accordingly.

The basis of the method is simply: 

   * if `p*q<n`, raise `q`. 
   * If `p*q>n`, lower `p`.
    
   * All while continuing where the last `p` and the last `q` left off.

The complex Fermat sieves, [described here](../Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf)
define the sequences for `p` and `q` via their complex mapping.

## Where CTM starts: at balance, not at the crossing

CTM begins at `R = ceil(sqrt(n))` -- `p` and `q` both as close to `sqrt(n)` as
their digit classes allow -- and expands outward from there.

An earlier version of this section argued it should instead begin where the
Fermat bracket stops, at the crossing `p0 = 2*j_max + 3`, and claimed the
crossing partitions the strip so that only factors "above" it are CTM's. **That
is wrong, and `n = 8051` shows why.**

`p` only ever descends, so the seed is the largest `p` the walk can ever test.
For `n = 8051 = 83 * 97`:

```
r = ceil(sqrt 8051) = 90,  m = 87,  S = 177,  j_max = 21,  p0 = 45
```

The factor is `p = 83`, and `83 > 45`. A walk seeded at the crossing starts
*past* the answer and moves away from it, so it can never find it. Seeded at
balance, in class `(3, 7)`:

```
p = largest  <= 90 with p % 10 == 3  ->  83
q = smallest >= 90 with q % 10 == 7  ->  97
83 * 97 = 8051 = n                       found on the first multiplication
```

This is the same failure as defect #1 below -- `real = r + 1` skipping `R = r`,
whose worked example is also `8051` -- reached by a different route. The crossing
is a real geometric point, but it is where the *135-degree traversal* leaves the
strip, not where CTM's reach begins.

From balance the walk covers `p` in `[3, r]` and `q` in `[r, n/3]`, which is
every factorisation of `n`. **CTM is complete on its own**; it does not own a
half. What the crossing does mark is where CTM stops being the cheaper of the
two, since its cost is `(q-p)/10` -- linear in the gap -- against the vertical
leg's quadratic `j_true`.

Running `--only=ctm --b1=1`, every case is found, balanced and unbalanced alike:

| n | p * q | multiplications |
| - | ----- | --------------- |
| 8051 | 83 * 97 | 117 |
| 143 | 11 * 13 | 1 |
| 2923 | 37 * 79 | 130 |
| 288419 | 379 * 761 | 38 |
| 1007509 | 503 * 2003 | 151 |
| 5115191 | 1597 * 3203 | 161 |
| 17821157 | 3019 * 5903 | 289 |
| 798607 | 101 * 7907 | 782 |
| 10967535067 | 104723 * 104729 | 56,622 |
| 978508015703 | 752867 * 1299709 | 323,663 |

## The walk, and what the published Julia does instead

The rule at the top of this page -- *raise `q`* when `p*q < n`, *lower `p`* when
`p*q > n* -- moves `p` and `q` **independently**, each by 10, with the other held
fixed. Those two moves are the ascend and the descend. Consequently `p` is
monotone down, `q` is monotone up, and `R = (p+q)/2` **oscillates**: there is no
monotone `R`.

The Julia below does something different. It bumps `TMRS` or `TMIS` by 10, which
shifts *both* `p` and `q` together:

```
p*q < n  ->  TMRS += 10  ->  p+10 AND q+10
p*q > n  ->  TMIS += 10  ->  p-10 AND q+10
```

so `q = R+i` climbs unconditionally and `q` becomes monotone. Traced side by side
on `n = 798607` the two visit different points. The C++ engine originally ported
the Julia and inherited its monotone `q`; it now implements the rule above.

The rule matters because its cost is `(q-p)/10` -- **linear** in the gap between
the factors, where the vertical Fermat leg's `j_true` is quadratic in it:

| n | p * q | CTM steps | vertical leg `j_true` |
| - | ----- | --------- | --------------------- |
| 309 | 3 * 103 | **10** | 35 |
| 798607 | 101 * 7907 | **781** | 3,110 |
| 1007509 | 503 * 2003 | **150** | 249 |
| 978508015703 | 752867 * 1299709 | 54,684 | **37,092** |

CTM wins as the factors separate and loses as they close up -- which is exactly
the partition this page specifies.

### Sieving the steps

`p` and `q` each keep their terminal digit under a `+-10` step, so a stream seeded
in an admissible `(p mod 10, q mod 10)` class stays in it for the whole walk. The
classes come from the Fermat sieve for this `N`'s classification: the `(R, i)`
table gives `p = R - i` and `q = R + i`, so

```
(p mod 10, q mod 10) = ((a10 - b10) mod 10, (a10 + b10) mod 10)
```

Verified against brute force for all ten odd residues mod 20:

| `N mod 20` | `(R,i)` classes | induced `(p,q)` classes |
| - | - | - |
| 1, 11 | 4 | (1,1) (3,7) (7,3) (9,9) |
| 3, 13 | 4 | (1,3) (3,1) (7,9) (9,7) |
| 7, 17 | 4 | (1,7) (3,9) (7,1) (9,3) |
| 9, 19 | 4 | (1,9) (3,3) (7,7) (9,1) |
| 5, 15 | 9 | (1,5) (3,5) (5,1) (5,3) (5,5) (5,7) (5,9) (7,5) (9,5) |

Several `(R,i)` classes collapse onto the same `(p,q)` class, so the list must be
de-duplicated -- stepping `p` and `q` independently needs the `(p,q)` classes, not
the `(R,i)` ones.

### Seeding q

`q` is seeded just *below* `n/p`, never above: from below the first test reads
`p*q < n` and the walk raises `q` into place, whereas from above it reads
`p*q > n` and drops `p` past the seed, skipping the top of the range.

Rounding `q` down into its class can put it below `p`, and that is harmless --
`p*q < n` there, so `q` climbs past `p` by itself. Guarding it with a skip
discards whole classes: for `n = 798607` it discarded `(1, 7)`, the one holding
`101 * 7907`, and the factor went unfound.

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
