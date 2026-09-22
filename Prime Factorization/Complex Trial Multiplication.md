# Complex Trial Multiplication
This method of factorization involves no division, rather tests products `p*q==n` and adjusts `p` or `q` accordingly.

The basis of the method is simply: 

   * if `p*q<n`, raise `q`. 
   * If `p*q>n`, lower `p`.
    
   * All while continuing where the last `p` and the last `q` left off.

The complex Fermat sieves, [described here](../Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf)
define the sequences for `p` and `q` via their complex mapping.

## Where CTM should start

CTM is the middle-of-the-strip method, so it should begin where the Fermat
bracket stops rather than at `ceil(sqrt(n)) + 1`.

With `r = ceil(sqrt(n))`, `m = r - 3` and `S = r + m`, the Sec. 4 traversal runs
down the 135-degree line through the key complex number `(r + m*i)` and halts at
`j_max = (n/S - 3)/2`, the last point before the key number turns red. Its
terminal point is

```
(R, i) = (r + j_max,  m - j_max)
```

That point is not arbitrary. The stop condition is `(R-i)(R+i) > n`, so the last
admissible point is by definition where the 135-degree line **crosses the factor
strip**, and the agreement tightens with n:

| n | terminal (R, i) | strip `i` at that R |
| - | --------------- | ------------------- |
| 309 | (21, **12**) | 11.489 |
| 8051 | (111, **66**) | 65.345 |
| 798607 | (1116, **669**) | 668.468 |
| 1000036000099 | (1250021, **750012**) | 750011.000 |

For `n = 309` the path is `j = 0..3` -- four steps, ending at `(21, 12)`.

**The crossing partitions the strip exactly.** Measured over a spread of
semiprimes, a factor whose `(R, i)` lies *above* the terminal point is CTM's;
one *below* it is the vertical Fermat scan's, and CTM started at the crossing
can never reach it (it would be walking away from the answer):

| n | factor | region |
| - | ------ | ------ |
| 309 | 3*103 | above -- CTM |
| 798607 | 101*7907 | above -- CTM |
| 8051 | 83*97 | below -- vertical Fermat |
| 17821157 | 3019*5903 | below -- vertical Fermat |

Starting at the crossing rather than at `r + 1` saves exactly the segment the
Fermat bracket already covered. On CTM's own region that is a **1.04x**
reduction in walk steps overall -- up to 1.16x when the factor is near the
crossing, near 1.00x when it is far up the strip, because CTM's cost is
dominated by the distance still to travel. The crossing's real value is as a
**specification of ownership**: it is the unique point at which the Fermat
bracket's coverage ends and CTM's begins, with no overlap and no gap.

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
