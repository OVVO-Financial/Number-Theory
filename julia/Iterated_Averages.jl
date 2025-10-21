function VF(n,step_size)
min_real=ceil(sqrt(n))
max_im= (n-9)/6

last_digit = mod(n,10)

iterated_average = [3, (min_real+n)/2]

#Create Iterated Average for log2(max imaginary value) number of points
for i in 1:log2(max_im)
  push!(iterated_average,(iterated_average[end]+n)/2)
end

deleteat!(iterated_average,1) #Drop the 3 from the initial assignment
cut=ceil(Int,log2(length(iterated_average)))  #Use only the first few, relationship breaks at higher iterated average points
resize!(iterated_average,cut)
likely_point = findmax(rem(iterated_average,1))[2] #Highest decimal remainder iterated average is most likely...
likely_point = max(likely_point,2)  #Can always decrease steps from first iterated average...
ascending_iterations = ceil(Int,iterated_average)
descending_iterations = floor(Int,iterated_average)

j=0
while true

i = likely_point     # use [for i in (cut:-1:1)] to sequentially check over all iterated average points

#for i in (cut:-1:1)

    d=descending_iterations[i]-j
      gcd_d=gcd1(d,n)
      if(gcd_d>1)
        factor_2=div(n,gcd_d)
        true_b=abs(gcd_d - factor_2)/2
        return(Dict("Realized (b)" => j, "True (b)" => convert(Int,true_b), "Iter" => i, "Iterated Average" => d, "FACTOR 1" => gcd_d, "FACTOR 2" => factor_2))
        end
      gcd_d1=gcd1(d-1,n)
      if(gcd_d1>1)
        factor_2=div(n,gcd_d1)
        true_b=abs(gcd_d1 - factor_2)/2
        return(Dict("Realized (b)" => j, "True (b)" => convert(Int,true_b), "Iter" => i, "Iterated Average" => d-1, "FACTOR 1" => gcd_d1, "FACTOR 2" => factor_2))
        end
      if last_digit==1 | last_digit==9
      gcd_d2=gcd1(d-2,n)
        if(gcd_d2>1)
        factor_2=div(n,gcd_d2)
        true_b=abs(gcd_d2 - factor_2)/2
        return(Dict("Realized (b)" => j, "True (b)" => convert(Int,true_b), "Iter" => i, "Iterated Average" => d-2, "FACTOR 1" => gcd_d2, "FACTOR 2" => factor_2))
        end
      end #last digit Fermat sieve for additional +- 2 gcd check

    a=ascending_iterations[i]+j
    if a < n - 2
      gcd_a=gcd1(a,n)
      if(gcd_a>1)
        factor_2=div(n,gcd_a)
        true_b=abs(gcd_a - factor_2)/2
        return(Dict("Realized (b)" => j, "True (b)" => convert(Int,true_b), "Iter" => i, "Iterated Average" => a, "FACTOR 1" => gcd_a, "FACTOR 2" => factor_2))
        end
      gcd_a1=gcd1(a+1,n)
      if(gcd_a1>1)
        factor_2=div(n,gcd_a1)
        true_b=abs(gcd_a1 - factor_2)/2
        return(Dict("Realized (b)" => j, "True (b)" => convert(Int,true_b), "Iter" => i, "Iterated Average" => a+1, "FACTOR 1" => gcd_a1, "FACTOR 2" => factor_2))
        end
      if last_digit==1 | last_digit==9
      gcd_a2=gcd1(a+2,n)
        if(gcd_a2>1)
        factor_2=div(n,gcd_a2)
        true_b=abs(gcd_a2 - factor_2)/2
        return(Dict("Realized (b)" => j, "True (b)" => convert(Int,true_b), "Iter" => i, "Iterated Average" => a+2, "FACTOR 1" => gcd_a2, "FACTOR 2" => factor_2))
        end
      end
    end #last digit Fermat sieve for additional +- 2 gcd check

#end #iff using for i in (cut:-1:1) loop

j += step_size  # Deterministic step_size setting should be j+= 5...Haven't run into many issues with j+= 10...
end #while
end #function
