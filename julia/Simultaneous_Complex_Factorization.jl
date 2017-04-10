function SCF{T<:Union{Int64,UInt64,Int128,UInt128,BigInt}}(n::T,modulo::T)
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

# ELIMINATE EXTRA SYNCHING SIEVE ENTRIES FOR NEW POSITIONAL MATCHING MULTIPLICATION & FERMAT REALS
    real_sieve = real_sieve[real_sieve.<=9]
    imaginary_sieve = imaginary_sieve[imaginary_sieve.<=9]

    TMRS = TM_real+(real_sieve - TM_real%10)
    TMIS = TM_imaginary+(imaginary_sieve - TM_imaginary%10)
    TMdRS = TM_real+(real_sieve - TM_real%10)
    TMdIS = TM_imaginary+(imaginary_sieve - TM_imaginary%10)
    lTMRS = length(TMRS)


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

# FIND QUADRATIC RESIDUES FOR MODULO
    residues = [0]
    for i in 1:div(modulo,2)
        push!(residues,(i^2)%modulo)
    end
    residues = unique(residues)

# SETS LOWER BOUNDARIES FOR (p) AND UPPER BOUNDARIES FOR (q) FOR EACH SIEVE ENTRY
    ceilings = div(n,3)
    ceiling_p = [ceilings,ceilings,ceilings,ceilings];ceiling_p_d=ceiling_p
    floor_q = [3,3,3,3];floor_q_d=floor_q


while (TMdIS[1] >= 0)
#COMPLEX TRIAL MULTIPLICATION
for i in 1:lTMRS
  p = TMRS[i] - TMIS[i]; p_d = TMdRS[i] - TMdIS[i]
  if(p>ceiling_p[i])
    TMIS[i] = TMIS[i] + 10
  end
  if(p<=1)
      TMRS[i] = TMRS[i] + 10
      p = TMRS[i] - TMIS[i]
  end
  if(p_d>ceiling_p[i])
    TMdIS[i] = TMdIS[i] - 10
  end
  if(p_d<=1)
      TMdRS[i] = TMdRS[i] - 10
      p = TMdRS[i] - TMdIS[i]
  end
  q = TMRS[i] + TMIS[i]; q_d = TMdRS[i] + TMdIS[i]
  if(q<floor_q[i])
    TMRS[i]=TMRS[i] + 10
  end
  if(q_d<floor_q[i])
    TMdRS[i]=TMdRS[i] - 10
  end
  N = p*q; N_d = p_d*q_d

  if (N == n) return("TM ascending",p, q) end; if (N_d == n) return("TM descending",p_d, q_d) end

  if (N > n)
    ceiling_p[i]=p
    TMIS[i] = TMIS[i] + 10
  else
    floor_q[i]=q
    TMRS[i] = TMRS[i] + 10
  end
  if (N_d < n)
    floor_q_d[i]=q_d
    TMdIS[i] = TMdIS[i] - 10
  else
    ceiling_p_d[i]=p_d
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
