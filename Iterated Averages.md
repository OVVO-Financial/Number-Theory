# Iterated Averages
The iterated averages routine takes advantage of a Fermat property originally described in [this paper](https://www.scribd.com/doc/298822749/Complex-Space-Factorization).
We are familiar with the Fermat factorization of the form `N = a^2 + b^2` resulting in the factors `(a-b),(a+b)`.
Usually `a` is the variable being tested, however, the iterated averages routine searches for `b`.

A brief example will illustrate:

For `n=8051`

Our first Fermat test of `a` begins with `ceil(sqrt(8051)) = 90` and low and behold, we have a factor after our first test (it's never this easy...).    
```
8051 = 90^2 - b^2
8051 = (90-7)(90+7)
```

Let us demonstrate what happens to the value of `b=7` after taking the iterated average and using it in a `gcd` check with `n`...

The first iterated average is

`(90+8051)/2 = 4070.5`

Solving for `b`
```
gcd(8051,(4070.5-b)) = (a-b)
gcd(8051,(4070.5-3.5)) = 83
```
`b=3.5`, half of the previous value `b=7`.  But wait, we can iterate the average again...

The **second** iterated average is

`(4070.5+8051)/2 = 6060.75`

Solving for `b`
```
gcd(8051,(6060.75-b)) = (a-b)
gcd(8051,(6060.75-1.75)) = 83
```
`b=1.75`, half of the previous value `b=3.5`.  But wait, we can iterate the average again...

The **third** iterated average is

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
The `Actual diff=` is the final value of `b` while the `true diff=` is the Fermat value of `b`.

`Iter` is the iterated average position used in the test and `Iterated Average` is the value (after difference) used in the `gcd`.

``` julia
@time VF(798607,5)
Actual diff=35  true diff=3903.0  Iter=4  Iterated Average=748713  FACTORS=101 7907.0
  0.000885 seconds (76 allocations: 3.156 KB)
```
```julia
@time VF(978508015703,5)
Actual diff=50745  true diff=273421.0  Iter=6  Iterated Average=963218792667  FACTORS=752867 1.299709e6
  0.014069 seconds (40.66 k allocations: 637.797 KB)
```


## Julia code
Below is the `julia` code.  I'm sure it can be optimized further!

``` julia
function VF(n,step_size)
min_real=ceil(sqrt(n))
max_im= (n-9)/6

last_digit = mod(n,10)

iterated_average = [3, (min_real+n)/2]

#Create Iterated Average for log2(n) number of points
for i in 1:log2(max_im)
  push!(iterated_average,(iterated_average[end]+n)/2)
end

deleteat!(iterated_average,1) #Drop the 3 from the initial assignment
cut=ceil(log2(length(iterated_average)))  #Use only the first few, relationship breaks at higher iterated average points
cut=convert(Int,cut)
resize!(iterated_average,cut)
likely_point = findmax(mod(iterated_average,1))[2] #Highest decimal iterated average is most likely...
likely_point = max(likely_point,2)  #Can always decrease steps from first iterated average...
ascending_iterations = ceil(iterated_average)
descending_iterations = floor(iterated_average)

j=0
while true

i = likely_point     # use [for i in (cut:-1:1)] to sequentially check over all iterated average points

    d=descending_iterations[i]-j; d=convert(Int,d)
    gcd_d=gcd(d,n)
    if(gcd_d>1)
    return(println("Actual diff=", j, "  true diff=", (gcd_d+n/gcd_d)/2 - gcd_d, "  Iter=", i, "  Iterated Average=", d , "  FACTORS=", gcd_d,  " ",n/gcd_d))
    end
    a=ascending_iterations[i]+j; a=convert(Int,a)
    gcd_a=gcd(a,n)
    if(gcd_a>1)
    return(println("Actual diff=", j, "  true diff=", (gcd_a+n/gcd_a)/2 - gcd_a, "  Iter=", i, "  Iterated Average=", a , "  FACTORS=", gcd_a, " ",n/gcd_a))
    end
    gcd_d1=gcd(d-1,n)
    if(gcd_d1>1)
      return(println("Actual diff=", j, "  true diff=", (gcd_d1+n/gcd_d1)/2 - gcd_d1, "  Iter=", i, "  Iterated Average=", d-1 , "  FACTORS=", gcd_d1, " ", n/gcd_d1))
    end
    gcd_a1=gcd(a+1,n)
    if(gcd_a1>1)
      return(println("Actual diff=", j, "  true diff=", (gcd_a1+n/gcd_a1)/2 - gcd_a1, "  Iter=", i, "  Iterated Average=", a+1 , "  FACTORS=", gcd_a1, " ", n/gcd_a1))
    end
    if last_digit==1 | last_digit==9
    gcd_d2=gcd(d-2,n)
    if(gcd_d2>1)
      return(println("Actual diff=", j, "  true diff=", (gcd_d2+n/gcd_d2)/2 - gcd_d2, "  Iter=", i, "  Iterated Average=", d-2 , "  FACTORS=", gcd_d2, " ",n/gcd_d2))
    end
    gcd_a2=gcd(a+2,n)
    if(gcd_a2>1)
      return(println("Actual diff=", j, "  true diff=", (gcd_a2+n/gcd_a2)/2 - gcd_a2, "  Iter=", i, "  Iterated Average=", a+2 , "  FACTORS=", gcd_a2, " ",n/gcd_a2))
    end
    end #last digit Fermat sieve for additional +- 2 gcd check

j += step_size  # Deterministic step_size setting should be j+= 5...Haven't run into many issues with j+= 10...
end #while
end #function
```
