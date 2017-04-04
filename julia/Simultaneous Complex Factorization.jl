# INTEGER NEWTON SQUARE ROOT APPROXIMATOR FOR FAST SQUARE NUMBER TEST
function Newton_sqrt(N)
    a = 1
    b = N
    while abs(a-b) > 1
        b = div(N,a)
        a = div((a + b), 2)
    end
    return a
end



# SIMULTANEOUS COMPLEX FACTORIZATION --- NEEDS PARALLELIZATION FOR 3 METHODS WITHIN SAME WHILE LOOP
function SCF(n)
    r = sqrt(n)
    max_im= (n-9)/6
    max_real= max_im + 3

    TD = 3

    Fermat_real = ceil(Int,r)

    TM_imaginary = Fermat_real - 2

# FIND TM REAL
    TM_real=Newton_sqrt(n+(TM_imaginary*TM_imaginary))

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

# SQUARE CHECK FERMAT SIEVE
    is_square_sieve = unique((imaginary_sieve.*imaginary_sieve) % 10)


# SYNC COMPLEX TRIAL MULTIPLICATION REAL TO SIEVE STARTING POINT
    last_digit_real = TM_real % 10
    if any(last_digit_real!=real_sieve)
      real_init_diff = real_sieve - last_digit_real
      real_init_diff = real_init_diff[real_init_diff.>=0][1]
      TM_real = TM_real + real_init_diff
    end

    last_digit_im = TM_imaginary % 10
    if any(last_digit_im!=imaginary_sieve)
      im_init_diff = imaginary_sieve - last_digit_im
      im_init_diff = im_init_diff[im_init_diff.>=0][1]
      TM_imaginary = TM_imaginary + im_init_diff
    end

# SYNC FERMAT REAL TO SIEVE STARTING POINT
    last_digit_real = Fermat_real % 10
    if any(last_digit_real!=real_sieve)
      real_init_diff = real_sieve - last_digit_real
      real_init_diff = real_init_diff[real_init_diff.>=0][1]
      Fermat_real = Fermat_real + real_init_diff
    end

while (TM_imaginary <= max_real-TD)
#COMPLEX TRIAL MULTIPLICATION
      p = TM_real - TM_imaginary
      q = TM_real + TM_imaginary
      N = p*q
      if (N == n) return("TM",p, q) end
      if (N > n)
        TM_imaginary = TM_imaginary + 2
      else
        TM_real = TM_real + 2
      end

# FERMAT DIFF OF SQUARES
      b2 = (Fermat_real*Fermat_real)-n
      if any(b2%10.==is_square_sieve)
      b_test = Newton_sqrt(b2)
        if b_test*b_test == b2
        return("Fermat",Fermat_real-b_test, Fermat_real+b_test)
      end
      end
      Fermat_real=Fermat_real+2

# TRIAL DIVISION
      if div(n,TD)*TD==n
        return("TD",TD,div(n,TD))
      else
        TD = TD +2
        if TD % 10 == 5
          TD = TD +2
        end
      end

end # WHILE CONDITION

end # SCF FUNCTION
