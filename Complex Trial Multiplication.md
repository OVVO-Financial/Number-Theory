# Complex Trial Multiplication
This method of factorization involves no division, rather tests products `p*q==n` and adjusts `p` or `q` accordingly.

The basis of the method is: 

   > * if `p*q<n`, raise `q`. 
   > * If `p*q>n`, lower `p`.
    
   > * All while continuing where the last `p` and the last `q` left off.

The complex Fermat sieves, [described here](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Fermat%20Sieve%20Using%20Complex%20Numbers.pdf)
define the sequences for `p` and `q` via their complex mapping.

## Examples
```julia
@time CTM(798607)
  0.000055 seconds (14 allocations: 4.953 KB)
(101,7907)
```
```julia
@time CTM(978508015703)
  0.000638 seconds (14 allocations: 4.953 KB)
(752867,1299709)
```

## Julia code
Below is the `julia` code.  The steps can be substantially reduced by restricting `p` or `q` within its repective sieve.  However, adding this array check seemed to severely diminish performance.

``` julia
function CTM(n)
    r = sqrt(n)
    max_im= (n-9)/6
    real = ceil(Int,r) 
    
# FERMAT SIEVES
    last_digit = n % 10

    if((n+1)%4==0)
      if(last_digit==1)
        real_sieve=[0,4,6,10]
        imaginary_sieve=[3,5,7,13]
      end
      if(last_digit==3)
        real_sieve=[2,8,12]
        imaginary_sieve=[1,9,11]
      end
      if(last_digit==7)
        real_sieve=[4,6,14]
        imaginary_sieve=[3,7,13]
      end
      if(last_digit==9)
        real_sieve=[0,2,8,10]
        imaginary_sieve=[1,5,9,11]
      end
    else
      if(last_digit==1)
        real_sieve=[1,5,9,11]
        imaginary_sieve=[0,2,8,10]
      end
      if(last_digit==3)
        real_sieve=[3,7,13]
        imaginary_sieve=[4,6,14]
      end
      if(last_digit==7)
        real_sieve=[1,9,11]
        imaginary_sieve=[2,8,12]
      end
      if(last_digit==9)
        real_sieve=[3,5,7,13]
        imaginary_sieve=[0,4,6,10]
      end
    end #FERMAT SIEVE

# Sync imaginary to imaginary_sieve sequence
    imaginary = imaginary_sieve[1]
    if imaginary == 0
      imaginary = imaginary_sieve[2]
    end

# Sync real to real_sieve sequence
    last_digit_real = real % 10
    if any(last_digit_real!=real_sieve)
      real_init_diff = real_sieve - last_digit_real
      real_init_diff = real_init_diff[real_init_diff.>0][1]
      real = real + real_init_diff
    end

    while (imaginary <= max_im)
      p = real - imaginary
      q = real + imaginary
   
      if (p * q == n) return(p, q) end
      
      if (p*q > n)
        imaginary = imaginary + 2
        # possible array check to verify any(imaginary!=imaginary_sieve) else increase imaginary
      else
        real = real + 2
        # possible array check to verify any(real!=real_sieve) else increase real
      end
end

end

```
