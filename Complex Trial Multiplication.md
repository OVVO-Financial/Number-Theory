# Complex Trial Multiplication
This method of factorization involves no division, rather tests products `p*q==n` and adjusts `p` or `q` accordingly.

The basis of the method is simply: 

   * if `p*q<n`, raise `q`. 
   * If `p*q>n`, lower `p`.
    
   * All while continuing where the last `p` and the last `q` left off.

The complex Fermat sieves, [described here](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf)
define the sequences for `p` and `q` via their complex mapping.

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
Below is the `julia` code.  There are 3 or 4 (depending on the Fermat sieve for that particular last digit) paired multiplications per iteration.  One immediate efficiency would be from the parallelization of these multiplications within the `while` loop.

``` julia
function CTM(n)
    r = sqrt(n)
    max_im= div((n-9),6)
    real = ceil(Int,r)

# FERMAT SIEVES
        last_digit = n % 10
    
        if((n+1)%4==0)
          if(last_digit==1)
            real_sieve=[0,4,0,10]
            imaginary_sieve=[3,5,7]
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
            real_sieve=[0,2,0,10]
            imaginary_sieve=[1,5,9]
          end
        else
          if(last_digit==1)
            real_sieve=[1,5,9]
            imaginary_sieve=[0,2,0]
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
            real_sieve=[3,5,7]
            imaginary_sieve=[0,4,0]
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
realS = (real+(real_sieve - real%10))
lrS = length(realS)

# MAKE SURE 0 IMAGINARY STARTS AT 10
imaginary_sieve[imaginary_sieve.==0] = 10

while (imaginary_sieve[1] <= max_im)
      p1 = realS[1] - imaginary_sieve[1]
      q1 = realS[1] + imaginary_sieve[1]
      N1 = p1*q1

      if (N1 == n) return(p1, q1) end

      if (N1 > n)
        imaginary_sieve[1] = imaginary_sieve[1] + 10
      else
        realS[1] = realS[1] + 10
      end

# 2nd MULTIPLICATION
    p2 = realS[2] - imaginary_sieve[2]
    q2 = realS[2] + imaginary_sieve[2]
    
    N2 = p2*q2
    
    if (N2 == n) return(p2, q2) end
    
    if (N2 > n)
      imaginary_sieve[2] = imaginary_sieve[2] + 10
    else
      realS[2] = realS[2] + 10
    end

# 3rd MULTIPLICATION
   p3 = realS[3] - imaginary_sieve[3]
   q3 = realS[3] + imaginary_sieve[3]
   
   N3 = p3*q3
   
   if (N3 == n) return(p3, q3) end
   
   if (N3 > n)
     imaginary_sieve[3] = imaginary_sieve[3] + 10
   else
     realS[3] = realS[3] + 10
   end

# 4th MULTIPLICATION (IF NECESSARY)
    if lrS > 3

    p4 = realS[4] - imaginary_sieve[4]
    q4 = realS[4] + imaginary_sieve[4]
    
    N4 = p4*q4
    
    if (N4 == n) return(p4, q4) end
    
    if (N4 > n)
     imaginary_sieve[4] = imaginary_sieve[4] + 10
    else
     realS[4] = realS[4] + 10
    end
      end

  end # while

end
```
