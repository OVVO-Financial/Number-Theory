# Iterated Averages
The iterated averages routine takes advantage of a Fermat property originally described in [this paper](https://www.scribd.com/doc/298822749/Complex-Space-Factorization).
We are familiar with the Fermat factorization of the form `N = a^2 - b^2` resulting in the factors `(a-b),(a+b)`.
Usually `a` is the variable being tested, however, the iterated averages routine searches for `b`.

#### A brief example will illustrate:

For `n=8051`

Our first Fermat test of `a` begins with `ceil(sqrt(8051)) = 90` and low and behold, we have a factor after our first test (it's never this easy...).    
```
8051 = 90^2 - b^2
8051 = (90-7)(90+7)
```

Let us demonstrate what happens to the value of `b=7` after taking the iterated average and using it in a `gcd` check with `n`...

### The first iterated average is

`(90+8051)/2 = 4070.5`

Solving for `b`
```
gcd(8051,(4070.5-b)) = (a-b)
gcd(8051,(4070.5-3.5)) = 83
```
`b=3.5`, half of the previous value `b=7`.  But wait, we can iterate the average again...

### The **second** iterated average is

`(4070.5+8051)/2 = 6060.75`

Solving for `b`
```
gcd(8051,(6060.75-b)) = (a-b)
gcd(8051,(6060.75-1.75)) = 83
```
`b=1.75`, half of the previous value `b=3.5`.  But wait, we can iterate the average again...

### The **third** iterated average is

`(6060.75+8051)/2 = 7055.875`

Solving for `b`
```
gcd(8051,(7055.875-b)) = (a-b)
gcd(8051,(7055.875-0.875)) = 83
```
`b=0.875`, half of the previous value `b=1.75` and less than 1 step was needed.

## Examples
### Input
`VF(n,step_size)` is the iterated averages function call.  `step_size = 5` is the deterministic setting.  Other `step_sizes` multiples of 5 will also yield results (significantly more quickly), but runs the risk of failing. 

### Output
The `Realized (b)` is the final value of `b` while the `True (b)` is the Fermat value of `b`.

`Iter` is the iterated average position used in the test and `Iterated Average` is the value (after difference) used in the `gcd`.

``` julia
@time VF(798607,5)
0.000029 seconds (29 allocations: 2.188 KB)
Dict{String,Int64} with 6 entries:
  "Realized (b)"     => 35
  "Iter"             => 4
  "Iterated Average" => 748713
  "True (b)"         => 3903
  "FACTOR 1"         => 101
  "FACTOR 2"         => 7907
```
```julia
@time VF(978508015703,5)
0.010549 seconds (30 allocations: 2.719 KB)
Dict{String,Int64} with 6 entries:
  "Realized (b)"     => 50745
  "Iter"             => 6
  "Iterated Average" => 963218792667
  "True (b)"         => 273421
  "FACTOR 1"         => 752867
  "FACTOR 2"         => 1299709
```


## Julia code
Below is the `julia` code.  I'm sure it can be optimized further / parallelized!

One optimization step was to use a binary `gcd` function discussed [here](https://github.com/JuliaLang/julia/blob/master/base/intfuncs.jl).
### binary gcd: `gcd1()`
``` julia
function gcd1{T<:Union{Int64,UInt64,Int128,UInt128}}(a::T, b::T)
    a == 0 && return abs(b)
    b == 0 && return abs(a)
    za = trailing_zeros(a)
    zb = trailing_zeros(b)
    k = min(za, zb)
    u = unsigned(abs(a >> za))
    v = unsigned(abs(b >> zb))
    while u != v
        if u > v
            u, v = v, u
        end
        v -= u
        v >>= trailing_zeros(v)
    end
    r = u << k
    # T(r) would throw InexactError; we want OverflowError instead
    r > typemax(T) && throw(OverflowError())
    r % T
end
```
### iterated averages: `VF()`
``` julia
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

i = likely_point     # use [for i in (cut:-1:1)] to sequentially check over all iterated average points...much slower.

#for i in (cut:-1:1)  # Possible distributed / parallelization from all "iterated_average" points

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
      end #last digit Fermat sieve for additional +- 2 gcd check
    end # n < limit condition
    
#end #iff using for i in (cut:-1:1) loop

j += step_size  # Deterministic step_size setting should be j+= 5...Haven't run into many issues with j+= 10...
end #while
end #function
```
