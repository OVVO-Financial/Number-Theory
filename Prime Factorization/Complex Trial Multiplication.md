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

**The start is the first red point of the 135-degree traversal**,

```
(R, i) = (r + j_max + 1,  m - j_max - 1)
```

which for `n = 8051` is `112 + 65i` -- `p = 47`, `q = 177`, `p*q = 8319 > n`, the
first point past the crossing. Because `R + i = r + m = S` all along the
traversal, both walks begin at `q = S`.

The Fermat bracket and trial division start at the **other end of the same
line**, `(r, m) = 90 + 87i` for `n = 8051`, and walk toward the crossing. The two
methods are indexed to one line from opposite ends.

### Sieving both directions

Every move is `+-10` on `R` or on `i`, so both the increases and the decreases
preserve the terminal digits. A stream synced into an admissible `(R mod 10,
i mod 10)` class at the start stays in it for the whole walk -- which is what
lets the paired Fermat sieve for this `N`'s classification apply to a
two-directional search unchanged.

For `n = 8051`: `n mod 4 = 3` (4k-1), last digit 1, so the paired tables give
`(R, i)` classes `(0,3) (0,7) (4,5) (6,5)`. The factor `83 * 97` is `R = 90`,
`i = 7` -- class `(0,7)`. Syncing `112 + 65i` into it gives `(120, 67)`, and
descending:

```
120+67i  p=53 q=187   9911 > n  -> R -= 10
110+67i  p=43 q=177   7611 < n  -> i -= 10
110+57i  p=53 q=167   8851 > n  -> R -= 10
...
 90+ 7i  p=83 q= 97   8051 = n     found, 9 steps
```

The ascending walk from the same point cannot reach it -- `q` only rises from
177, and the factor is at `q = 97`. That is exactly why both are needed.

### Measured

`--only=ctm --b1=1`, multiplications to the factor:

| n | p * q | mults |
| - | ----- | ----- |
| 143 | 11 * 13 | 9 |
| 2923 | 37 * 79 | 245 |
| 8051 | 83 * 97 | 690 |
| 798607 | 101 * 7907 | 612 |
| 10963 | 19 * 577 | 863 |
| 5115191 | 1597 * 3203 | 16,517 |
| 1007509 | 503 * 2003 | 24,578 |
| 288419 | 379 * 761 | 24,845 |
| 17821157 | 3019 * 5903 | 33,023 |
| 10967535067 | 104723 * 104729 | 96,960 |
| 978508015703 | 752867 * 1299709 | 575,803 |

### Two things that must not be added

**No `p < 3` early exit.** `p*q < n` when `p` is tiny, so the rule moves `i` and
raises `p` by 10 again -- the walk recovers. Breaking discards the class first,
which loses `n = 143`, whose seed passes through `p = 1` one step before
`p = 11`.

**Not one walk.** An earlier version of the C++ engine had only the ascending
walk. From the crossing, ascending can never reach a factor whose `q` is below
`S` -- 8051's `q = 97` against `S = 177` -- so it could not factor 8051 at all.

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
