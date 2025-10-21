# Naive Grouping
The naive grouping routine is as stated, *a naive grouping method for factorization.*
It is deterministic and can be distributed along any number of "seed groups".
Like factors, one of 3 *known* seed groups resides below `0.5*sqrt(n)`.

This method is loosely related to the [Ramanujan partition function](https://en.wikipedia.org/wiki/Partition_(number_theory)) whereby the relevant partitions of `n` are those whose summands differ by 1.  Those relevant partitions are then further partitioned to identify factors of `n`.

## The method is best described with an elementary example:

1. For `n = 187` if we start with a seed group of **2**, we would have **1** group of **93** and **1** group of **94**.
    > * The summands (e.g., 93 and 94 generated from `187/2 = 93.5`) are to have a maximum difference of 1 (hence naive grouping!).
    > * Checking our total: `1 * 93 + 1 * 94 = 187`
2. Neither of the groups will evenly parse the other summands so we double the seed group and halve the summands.
3. Now we have **4** groups in total, **3** groups of **47** and **1** groups of **46**.
    > * Checking our total: `3 * 47 + 1 * 46 = 187`
4. The **46** cannot be distributed evenly among the **3** groups, so we repeat.
5. Now we have **8** groups in total, **5** groups of **23** and **3** groups of **24**.
    > * Checking our total: `5 * 23 + 3 * 24 = 187`
6. The **24** cannot be distributed evenly among the **5** groups, nor can the **23** be evenly distributed among the **3** groups, so we repeat.
7. Now we have **16** groups in total, **11** groups of **12** and **5** groups of **11**.
    > * Checking our total: `11 * 12 + 5 * 11 = 187`
8. The **5** groups of **11** can be evenly distributed to the **11** groups of **12**, leaving us with **11** groups of **17**, our factors for `n = 187`.


## And further described with another example:

1. For `n = 798,607` if we start with a seed group of **59**, we would have **17** groups of **13,535** and **42** groups of **13,536**.
    > * The summands (e.g., 13,535 and 13,536 generated from `798607/59 = 13535.71`) are to have a maximum difference of 1 (hence naive grouping!).
    > * Checking our total: `17 * 13535 + 42 * 13536 = 798607`
2. Neither the **17** nor the **42** will evenly parse the other summands so we double the seed group and halve the summands.
    > * `mod(13535,42) > 0` and `mod(13536,17) > 0`.  
3. Now we have **118** groups in total, **17** groups of **6,767** and **101** groups of **6,768**.
    > * Checking our total: `17 * 6767 + 101 * 6768 = 798607`
4. The **6,767** can be distributed among the **101** groups, resulting in factors of **101** groups of **7,907**.
    > * `mod(6767,101) == 0`. 

There are minimal division operations (one per seed group) performed, rather elementary doubling and halving rules resulting in significantly reduced `mod` checks versus the original number `n`.

## Examples
``` julia
@time NG(798607)
  0.000025 seconds (5 allocations: 192 bytes)
  (101,7907)
```

``` julia
@time NG(978508015703)
  0.008523 seconds (5 allocations: 192 bytes)
  (752867,1299709)
```

## Julia code
Below is the `julia` code.  I'm sure it can be optimized further / parallelized!
``` julia
function NG(n)
r = ceil(sqrt(n))
limit = floor(Int, r/2)
  if (limit & 1)==0
    limit = limit - 1
  end
rr = floor(Int, limit/2)
  if (rr & 1)==0
    rr = rr -1
  end
i = rr   # Possible distribution / parallelization, using multiple "Seed groups".  "rr" is worst case point for single "Seed group".
k = i+2

while true
# i seed group decreasing per iteration
  summand_floor = div(n,i)
  if summand_floor*i==n return (i, summand_floor) end
  summand_ceiling = summand_floor+1
  group_ceiling = mod(n,i)
  group_floor = i - group_ceiling
  
# k seed group increasing per iteration
  summand_floor2 = div(n,k)
  if summand_floor2*k==n return (k, summand_floor2) end
  summand_ceiling2 = summand_floor2+1
  group_ceiling2 = mod(n,k)
  group_floor2 = k - group_ceiling2

  while true
  if (summand_floor & 1)==0
    group_floor=(2*group_floor)+group_ceiling
  end

  if (summand_ceiling & 1)==0
    group_ceiling=(2*group_ceiling)+group_floor
  end
  
  if summand_floor < summand_ceiling
    summand_floor=div(summand_floor,2) 
    summand_ceiling=cld(summand_ceiling,2)
  else
    summand_floor=cld(summand_floor,2)
    summand_ceiling=div(summand_ceiling,2)
  end

  if ((summand_floor < group_ceiling) | (summand_ceiling < group_floor)) 
      break
  end

  if (summand_floor2 & 1)==0
    group_floor2=(2*group_floor2)+group_ceiling2
  end
  if (summand_ceiling2 & 1)==0
    group_ceiling2=(2*group_ceiling2)+group_floor2
  end
  if summand_floor2< summand_ceiling2
    summand_floor2=div(summand_floor2,2)
    summand_ceiling2=cld(summand_ceiling2,2)
  else
    summand_floor2=cld(summand_floor2,2)
    summand_ceiling2=div(summand_ceiling2,2)
  end

  if(group_ceiling>1 && group_floor>1)
    if(div(summand_floor,group_ceiling)*group_ceiling==summand_floor)
        return (group_ceiling , div(n,group_ceiling))
        end

        if(div(summand_ceiling,group_floor)*group_floor==summand_ceiling)
          return (group_floor , div(n,group_floor))
          end
  end

  if(group_ceiling2>1 && group_floor2>1)
    if(div(summand_floor2,group_ceiling2)*group_ceiling2==summand_floor2)
        return (group_ceiling2 , div(n,group_ceiling2))
      end
      if(div(summand_ceiling2,group_floor2)*group_floor2==summand_ceiling2)
        return (group_floor2 , div(n,group_floor2))
      end
  end

end #inner while for i, k expansion
i = i - 2
k = k + 2
  if i<3
    return ("PRIME")
  end
  if k>limit
    return ("PRIME")
  end
end #outer while for i, k sequence
end
```
