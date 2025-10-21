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
    square_sieve = unique((imaginary_sieve.*imaginary_sieve) % 10)

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
    ceiling_p = [ceilings,ceilings,ceilings,ceilings];ceiling_p2=ceiling_p
    floor_q = [3,3,3,3];floor_q2=floor_q

while (TMdRS[1] >= 0)

#COMPLEX TRIAL MULTIPLICATION
for i in 1:lTMRS
  a=TMRS[i];b=TMIS[i]
  a2=TMdRS[i];b2=TMdIS[i]

  p = a - b; p2 = a2 - b2
  if(p>ceiling_p[i])
    TMIS[i] = TMIS[i] + 10
  end
  if(p2<ceiling_p2[i])
    TMdIS[i] = TMdIS[i] - 10
  end
  if(p<=1)
      TMRS[i] = TMRS[i] + 10
  end
  if(p2<=1)
      TMdIS[i] = TMdIS[i] + 10
  end
  a=TMRS[i];b=TMIS[i]
  a2=TMdRS[i];b2=TMdIS[i]

  p = a - b; p2 = a2 - b2
  q = a + b; q2 = a2 + b2

  if(q<floor_q[i])
    TMRS[i]=TMRS[i] + 10
  end
  if(q2>floor_q2[i])
    TMdRS[i] = TMdRS[i] - 10
  end

  q = TMRS[i] + TMIS[i]; q2 = TMdRS[i] + TMdIS[i]

  a=TMRS[i];b=TMIS[i]
  a2=TMdRS[i];b2=TMdIS[i]

  z = complex(a,b)^2
  z2 = complex(a2,b2)^2

  N = real(z);N2=real(z2)

  if (N == n) return("TM ascending",p, q) end; if (N2 == n) return("TM descending",p2, q2) end

  if (N > n)
    ceiling_p[i]=p
    TMIS[i] = TMIS[i] + 10
  else
    floor_q[i]=q
    TMRS[i] = TMRS[i] + 10
  end

  if (N2 < n)
    floor_q2[i]=q2
    TMdIS[i] = TMdIS[i] - 10
  else
    ceiling_p2[i]=p2
    TMdRS[i] = TMdRS[i] - 10
  end

  if ((TMdRS[i]<TMdIS[i]) | ((TMdRS[i]-TMdIS[i])<=1))
    TMdIS[i]=TMdIS[i]-10
  end

 if (TMdIS[i]<1)
    TMdIS[i]=TMdIS[i]+20
  end


end #for

# FERMAT DIFF OF SQUARES
  b2s = (F_realS.*F_realS) - n

  if length(intersect(b2s%10,square_sieve))>0
  if length(intersect(b2s%modulo,residues))>0

        z= intersect(b2s%10,square_sieve).==b2s%10
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
