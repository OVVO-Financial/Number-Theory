# Simultaneous Complex Factorization
This page presents the annotated `julia` code from the discussion / visualization on factors in the complex space [presented here](https://github.com/OVVO-Financial/Number-Theory/blob/master/Complex%20space.md).

## `julia` code
### Step 1
In the introductory step we find the `Newton_sqrt` of the number `(n)` since no precision is needed. 
* `Newton_sqrt` is an integer based function, [available here](https://github.com/OVVO-Financial/Number-Theory/blob/Prime-Factorization/julia/Newton_Square_Root.jl). 

We also define the `maximum imaginary number`, `maximum real number`, the trial division starting point and the minimum Fermat real number `Fermat_real`.
```julia
function SCF{T<:Union{Int64,UInt64,Int128,UInt128,BigInt}}(n::T,modulo::T)
    r = Newton_sqrt(n)
    max_im= div((n-9),6)
    max_real= max_im + 3
    TD = 3
    Fermat_real = r+1
```

### Step 2
In step 2, we define the Fermat sieves.  Quite simply, we match known facts about `real` and `imaginary` values corresponding to `a` and `b` in Fermat's method.  If `n = 4*k +1` we know the `real` part of the complex number is **even**.

Using the even or odd distinction, we can further determine the combinations of `real` and `imaginary` values to the complex numbers that will result in the last digit of `n`.
* 8051 is of the form `4k + 1`.  The only even numbers with +/- odd counterparts that result in a 1 are `0,4,6`. For example `(20-3)(20+3) = (17)(23) = 391`, last digit equals 1.
* You cannot achieve that with `2 or 8` and +/- the same odd number in this example.
```julia
# FERMAT SIEVES
    last_digit = n % 10

    if((n+1)%4==0)
      if(last_digit==1)
        Fermat_real_sieve=[4,6,0,10]
        real_sieve=[0,4,6,0,10]
        imaginary_sieve=[3,5,5,7]
        end
      if(last_digit==3)
        real_sieve=[2,8,2,8]
        Fermat_real_sieve=real_sieve
        imaginary_sieve=[1,9,9,1]
        end
      if(last_digit==7)
        real_sieve=[4,6,4,6]
        Fermat_real_sieve=real_sieve
        imaginary_sieve=[3,7,7,3]
        end
      if(last_digit==9)
        Fermat_real_sieve=[0,2,8,10]
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
        Fermat_real_sieve=real_sieve
    end #FERMAT SIEVE
```
### Step 3
Using our known `imaginary_sieve` values from step 2, we can square them to find the corresponding required square for `a^2 - n`.
* The `imaginary_sieve` to 8051 is `[3,5,7]`.  Squaring this yields `[9,25,49]`, thus all `a^2 - n` have to end in `9 or 5` and this reduces the number of square checks.
```julia
# SQUARE CHECK FERMAT SIEVE
    square_sieve = unique((imaginary_sieve.*imaginary_sieve) % 10)
```
### Step 4
This step is required to sync our [complex trial mulitiplication](https://github.com/OVVO-Financial/Number-Theory/blob/Prime-Factorization/Complex%20Trial%20Multiplication.md) `real` and `imaginary` components to the middle of the [complex space](https://github.com/OVVO-Financial/Number-Theory/blob/master/Complex%20space.md), while aligning with the respective sieves.
* Example: if our complex space starting point estimate ends in a 7 and we need a `[0,4,6]`, we shift up to the next number ending in `0`.

```julia
# SYNC COMPLEX TRIAL MULTIPLICATION IMAGINARY TO SIEVE STARTING POINT
    z=sqrt(n+(max_im*max_real)im)
    TM_imaginary = floor(Int,imag(z))
    last_digit_im = TM_imaginary % 10
    if any(last_digit_im!=imaginary_sieve)
      im_init_diff1 = imaginary_sieve[imaginary_sieve.<=9] - last_digit_im # need a this copy for later
      im_init_diff = imaginary_sieve - last_digit_im
      if any(im_init_diff.>=0)
      im_init_diff = im_init_diff[im_init_diff.>=0][1]
      else
      im_init_diff = im_init_diff[end]
      end
      TM_imaginary = TM_imaginary + im_init_diff
    end

# FIND & SYNC COMPLEX TRIAL MULTIPLICATION REAL
    TM_real = floor(Int,real(z))
    if(TM_real < TM_imaginary) TM_real = TM_imaginary + 3 end
      last_digit_real = TM_real % 10
      if any(last_digit_real!=real_sieve)
        real_init_diff = real_sieve - last_digit_real
        if any(real_init_diff.>=0)
        real_init_diff = real_init_diff[real_init_diff.>=0][1]
      else
        real_init_diff = real_init_diff[end]
      end
        TM_real = TM_real + real_init_diff
    end
```
### Step 5
This step is required to eliminate the `10`'s that were inserted into the `real_sieve` and `imaginary_sieve` when a `0` was present to find the shift difference.

We also increment each array entry if the difference between `real` and `imaginary` is `< 1` to ensure a proper start.
```julia
# ELIMINATE EXTRA SYNCHING SIEVE ENTRIES FOR NEW POSITIONAL MATCHING MULTIPLICATION & FERMAT REALS
    real_sieve = real_sieve[real_sieve.<=9]
    imaginary_sieve = imaginary_sieve[imaginary_sieve.<=9]

    TMRS = TM_real+(real_sieve - TM_real%10)
    TMIS = TM_imaginary+(imaginary_sieve - TM_imaginary%10)
    TMdRS = TM_real+(real_sieve - TM_real%10)
    TMdIS = TM_imaginary+(imaginary_sieve - TM_imaginary%10)
    lTMRS = length(TMRS)

# AVOID ANY INITIAL MIS-SYNCH
for i in 1:lTMRS
  if((TMRS[i]-TMIS[i])<=1)
    TMRS[i]=TMRS[i]+10
  end
  if((TMdRS[i]-TMdIS[i])<=1)
    TMdIS[i]=TMdIS[i]-10
  end
end
```
### Step 6
This step shifts our `Fermat_real` to the nearest corresponding `real_sieve` value.
```julia
# SYNC FERMAT REAL TO SIEVE STARTING POINT
    last_digit_real = Fermat_real % 10
    if any(last_digit_real!=Fermat_real_sieve)
      real_init_diff = Fermat_real_sieve - last_digit_real
      if any(real_init_diff.>=0)
        real_init_diff = real_init_diff[real_init_diff.>=0][1]
      else
        real_init_diff = real_init_diff[end]
      end
      Fermat_real = Fermat_real + real_init_diff
    end
```
### Step 7
Once shifted, we find the other values in the `real_sieve` and create entries to an array.
* Using 8051, if our `real sieve` shifted value was 90, we include 94 and 96 because the `real sieve = [0,4,6]` for 8051.
```julia
# DEFINE FERMAT REALS - EVERY INITIAL MEMBER OF SIEVE
    F_realS = unique(Fermat_real+(Fermat_real_sieve - Fermat_real%10))
    lFrS = length(F_realS)
```
### Step 8
In an effort to reduce square checks further, we generate the `quadratic_residues` for the given entered `modulo`.  If a `Fermat_real` does not have the corresponding residue, it will not be squared and tested.
```julia
# FIND QUADRATIC RESIDUES FOR MODULO
    residues = [0]
    for i in 1:div(modulo,2)
        push!(residues,(i^2)%modulo)
    end
    residues = unique(residues)
```
### Step 9
During trial multiplication, if `p * q < n` we can set lower boundaries on `q`.  Conversely, if `p * q > n` we can set upper boundaries on `p`.  This helps eliminate combinations of `real` and `imaginary` that do not fall within these boundaries.
```julia
# SETS LOWER BOUNDARIES FOR (p) AND UPPER BOUNDARIES FOR (q) FOR EACH SIEVE ENTRY 
    ceilings = div(n,3)
    ceiling_p = [ceilings,ceilings,ceilings,ceilings];ceiling_p_d=ceiling_p
    floor_q = [3,3,3,3];floor_q_d=floor_q
```
### Step 10
Run the tests!  Everything is contained within the `while` loop.

There are 3 tests that occur for each `while` increment:
1. **Complex Trial Multiplication**: For each value in our `real` and `imaginary` array, we multiply their corresponding `real` numbers.  The arrays are paired such that `real[1]` only has to be multiplied using `imaginary[1]`, for the length of the arrays.
    * This can likely use a parallelization since there are only 4 entries in the arrays.
    * Each `while` increment increases the `real` and `imaginary` arrays by 10, since they are aligned with their respective sieves.
2. **Fermat Difference of Squares**: For each value in our `Fermat_real` array we test its `quadratic_residue`, then see if its difference with `n` is in our `square_sieve`
and if those conditions are met we then test the resulting difference if it is a square with the `Newton_sqrt` function, again, since precision is not necessary.
    * Each `while` increment increases the `Fermat_real` array by 10.  If a `quadratic residue` does not exist, it will keep increasing by 10 until one is found.  This is contained in a `for` loop.
3. **Trial Division**: We eliminate the small number of early primes as our sieve and eliminate any odd numbers ending in 5.
    * Each `while` increment increases the `TD` number by 2.  Obviously a more expansive prime check is desired, but this comes at a significant performance cost.
    * Since precision is not needed, we use a `div(n,TD)*TD==n` check.  Also, we didn't need the exact `modulo` of this calculation so it's a really a truncated `mod()` routine.

```julia
while (TMdIS[1] >= 0)
#COMPLEX TRIAL MULTIPLICATION
for i in 1:lTMRS
  p = TMRS[i] - TMIS[i]; p_d = TMdRS[i] - TMdIS[i]
  q = TMRS[i] + TMIS[i]; q_d = TMdRS[i] + TMdIS[i]
  N = p*q; N_d = p_d*q_d

  if (N == n) return("TM ascending",p, q) end; if (N_d == n) return("TM descending",p_d, q_d) end

  if (N > n)
    TMIS[i] = TMIS[i] + 10
  else
    TMRS[i] = TMRS[i] + 10
  end
  if (N_d < n)
    TMdIS[i] = TMdIS[i] - 10
  else
    TMdRS[i] = TMdRS[i] - 10
  end
end #for

# FERMAT DIFF OF SQUARES
  b2s = (F_realS.*F_realS) - n

  if length(intersect(b2s%10,is_square_sieve))>0
  if length(intersect(b2s%modulo,residues))>0

        z= intersect(b2s%10,is_square_sieve).==b2s%10
        b2s=b2s[z];  Factor_real = F_realS[z]
          if length(b2s)>0
            for i in 1:length(b2s)
              b_test = Newton_sqrt(b2s[i])
              if b_test*b_test == b2s[i]
                return("Fermat",Factor_real[i]-b_test, Factor_real[i]+b_test)
              end
            end #for
          end

      end #residues

    end #if b2s in square sieve

    F_realS = F_realS + 10
    for i in div(modulo,10)
      if length(intersect(F_realS%modulo,Fermat_real_sieve))==0
          F_realS= F_realS + 10
      else
          break
      end
    end
    F_realS = F_realS[F_realS.<max_real]

# TRIAL DIVISION
    # QUICK PRIMALITY CHECK - BASIC, BUT ELIMINATES A LARGE PERCENTAGE (EVEN FOR LARGE N)
    if any(TD.%[3,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109]!=0)
      if div(n,TD)*TD==n
        return("TD",TD,div(n,TD))
      else
        TD = TD +2
        if TD % 10 == 5
          TD = TD +2
        end
      end
    end # PRIMALITY TEST

end # WHILE CONDITION

end # SCF FUNCTION
```

Any and all comments / suggestions are welcomed to increase the efficiency and capability of this function.  

### Thanks for your interest!
