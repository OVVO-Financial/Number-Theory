function Newton_sqrt{T<:Union{Int64,UInt64,Int128,UInt128,BigInt}}(N::T)
    a = 1
    b = N
    while abs(a-b) > 1
        b = div(N,a)
        a = div((a + b), 2)
    end
    return a
end




function SCF{T<:Union{Int64,UInt64,Int128,UInt128,BigInt}}(n::T)
    r = Newton_sqrt(n)
    max_im= div((n-9),6)
    max_real= max_im + 3

    TD = 3

    Fermat_real = r+1

# FERMAT SIEVES
    last_digit = n % 10

    if((n+1)%4==0)
      if(last_digit==1)
        Fermat_real_sieve=[4,6,0,10]
        real_sieve=[0,4,0,10]
        imaginary_sieve=[3,5,7]
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
        Fermat_real_sieve=real_sieve
    end #FERMAT SIEVE

# SQUARE CHECK FERMAT SIEVE
    is_square_sieve = unique((imaginary_sieve.*imaginary_sieve) % 10)

# SYNC COMPLEX TRIAL MULTIPLICATION IMAGINARY TO SIEVE STARTING POINT
    TM_imaginary = div((max_im+Fermat_real),2)
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
      TM_real=Newton_sqrt(n+(TM_imaginary*TM_imaginary))
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

# ELIMINATE EXTRA SYNCHING SIEVE ENTRIES FOR NEW POSITIONAL MATCHING MULTIPLICATION & FERMAT REALS
        real_sieve = real_sieve[real_sieve.<=9]
        imaginary_sieve = imaginary_sieve[imaginary_sieve.<=9]

        TM_realS = TM_real+(real_sieve - TM_real%10)
        TM_imagS = TM_imaginary+im_init_diff1
        lTMrS = length(TM_realS)

# DESCENDING COMPLEX TRIAL MULTIPLICATION STARTING POINTS
        TM_imagS_desc=TM_imagS;TM_realS_desc=TM_realS

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

# DEFINE FERMAT REALS - EVERY INITIAL MEMBER OF SIEVE
      F_realS = unique(Fermat_real+(Fermat_real_sieve - Fermat_real%10))
      lFrS = length(F_realS)

#return([Fermat_real,real_sieve,F_realS])  #- good to here

while (TM_imagS_desc[1] >= 0)
#COMPLEX TRIAL MULTIPLICATION
  p1 = TM_realS[1] - TM_imagS[1]; p1_d = TM_realS_desc[1] - TM_imagS_desc[1]
  q1 = TM_realS[1] + TM_imagS[1]; q1_d = TM_realS_desc[1] + TM_imagS_desc[1]
  N1 = p1*q1; N1d = p1_d*q1_d

  if (N1 == n) return(p1, q1) end; if (N1d == n) return(p1_d, q1_d) end

  if (N1 > n)
    TM_imagS[1] = TM_imagS[1] + 10
  else
    TM_realS[1] = TM_realS[1] + 10
  end
  if (N1d < n)
    TM_imagS_desc[1] = TM_imagS_desc[1] - 10
  else
    TM_realS_desc[1] = TM_realS_desc[1] - 10
  end

# 2nd PAIRED MULTIPLICATION
  p2 = TM_realS[2] - TM_imagS[2]; p2_d = TM_realS_desc[2] - TM_imagS_desc[2]
  q2 = TM_realS[2] + TM_imagS[2]; q2_d = TM_realS_desc[2] + TM_imagS_desc[2]
  N2 = p2*q2; N2d = p2_d*q2_d

  if (N2 == n) return(p2, q2) end; if (N2d == n) return(p2_d, q2_d) end

  if (N2 > n)
    TM_imagS[2] = TM_imagS[2] + 10
  else
    TM_realS[2] = TM_realS[2] + 10
  end
  if (N2d < n)
    TM_imagS_desc[2] = TM_imagS_desc[2] - 10
  else
    TM_realS_desc[2] = TM_realS_desc[2] - 10
  end

# 3rd PAIRED MULTIPLICATION
  p3 = TM_realS[3] - TM_imagS[3]; p3_d = TM_realS_desc[3] - TM_imagS_desc[3]
  q3 = TM_realS[3] + TM_imagS[3]; q3_d = TM_realS_desc[3] + TM_imagS_desc[3]
  N3 = p3*q3; N3d = p3_d*q3_d

  if (N3 == n) return(p3, q3) end; if (N3d == n) return(p3_d, q3_d) end

  if (N3 > n)
    TM_imagS[3] = TM_imagS[3] + 10
  else
    TM_realS[3] = TM_realS[3] + 10
  end
  if (N3d < n)
    TM_imagS_desc[3] = TM_imagS_desc[3] - 10
  else
    TM_realS_desc[3] = TM_realS_desc[3] - 10
  end

# 4th PAIRED MULTIPLICATION (IF NECESSARY)
if lTMrS > 3
  p4 = TM_realS[4] - TM_imagS[4]; p4_d = TM_realS_desc[4] - TM_imagS_desc[4]
  q4 = TM_realS[4] + TM_imagS[4]; q4_d = TM_realS_desc[4] + TM_imagS_desc[4]
  N4 = p4*q4; N4d = p4_d*q4_d

  if (N4 == n) return(p4, q4) end; if (N4d == n) return(p4_d, q4_d) end

  if (N4 > n)
    TM_imagS[4] = TM_imagS[4] + 10
  else
    TM_realS[4] = TM_realS[4] + 10
  end
  if (N4d < n)
    TM_imagS_desc[4] = TM_imagS_desc[4] - 10
  else
    TM_realS_desc[4] = TM_realS_desc[4] - 10
  end
end #4th MULTIPLICATION CONDITION

# FERMAT DIFF OF SQUARES
# 1st FERMAT CHECK
      b2s = (F_realS.*F_realS) - n
      
      if any((b2s%10).!=[is_square_sieve])
      z= ((b2s%10).!=[is_square_sieve])
      b2s=b2s[z]
      b2s=b2s[b2s.>0]
      if length(b2s)>0
        for i in 1:length(b2s)
          b_test = Newton_sqrt(b2s[i])
          if b_test*b_test == b2s[i]
            return("Fermat",Fermat_real-b_test, Fermat_real+b_test)
          end
        end #for
        end
      end #if b2s in square sieve
      F_realS = F_realS + 10
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
