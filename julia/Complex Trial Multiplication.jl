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

      N = p*q

      if (N == n) return(p, q) end

      if (N > n)
        imaginary = imaginary + 2
      else
        real = real + 2
      end
end

end
