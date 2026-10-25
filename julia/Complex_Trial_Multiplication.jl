function CTM(n)
    r = Newton_sqrt(n)
    max_im= div((n-9),6)

# FIX 1: start at r, not r+1.
#
# The factor of a very balanced semiprime sits at R = ceil(sqrt(n)) exactly
# (the j = 0 case).  Starting at r+1 syncs that stream to the NEXT admissible
# R in its digit class, which is already past the answer, so the factor is
# unreachable for good.  n = 8051 (83*97, R = 90 = ceil(sqrt 8051)) and
# n = 143 (11*13) both fail this way.
    real = r

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

# FIX 3: no branch above assigns the sieve when n % 10 is 5 (or n is even),
# so real_sieve / imaginary_sieve are undefined and Julia throws
# UndefVarError.  5 | n is structurally different -- 9 admissible (R,i)
# classes per parity lane instead of 4, because exactly one of p, q carries
# the factor 5 -- see Number Theory Papers/Complete_Fermat_Sieve_Verified.pdf.
    if !(last_digit in [1,3,7,9])
      error("CTM: n must be coprime to 10; strip factors of 2 and 5 first")
    end

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


# FIX 2: retire streams individually.
#
# The loop condition tested only TMIS[1], so once stream 1 ran past max_im
# every other stream was cut off with it -- even one still short of its
# target.  n = 314187 (3*104729) needs i = 52363 on stream 4, right at
# max_im, and was terminated early.  Each stream now retires on its own.
#
# The per-iteration print([ceiling_p,floor_q]) is also removed: it was
# debug output on the hot path.
live = [true for _ in 1:lTMRS]

while any(live)
  for i in 1:lTMRS
    if !live[i] ; continue ; end
    if TMIS[i] > max_im ; live[i] = false ; continue ; end
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
