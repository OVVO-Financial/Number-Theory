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
ceilings = div(n,3)
ceiling_p = [ceilings,ceilings,ceilings,ceilings]
floor_q = [3,3,3,3]


while (TMIS[1] <= max_im)
  print([ceiling_p,floor_q])
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

    if (N == n) return(p, q,i) end

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
