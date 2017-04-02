# Naive Grouping
The naive grouping routine is as stated, *a naive grouping method for factorization.*
It is deterministic and can be distributed along any number of "seed groups".
One of 3 known seed groups resides below `0.5*sqrt(n)`.

The method is best described with an example:

1. For `n = 798,607` if we start with a seed group of **59**, we would have **17** groups of **13,535** and **42** groups of **13,536**.
    > * The group counts (e.g., 13,535 and 13,536 generated from `798607/59 = 13535.71`) are to have a maximum difference of 1 (hence naive grouping!).
    > * Checking our total: `17 * 13535 + 42 * 13536 = 798607`
2. Neither the **17** nor the **42** will evenly parse the other grouping counts so we double the seed group and halve the counts.
    > * `mod(13535,42) > 1` and `mod(13536,17) > 1`.  
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
Below is the `julia` code.  I'm sure it can be optimized further!
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
i = rr
k = i+2

while true
  count_floor = div(n,i); count_ceiling = count_floor+1
  group_ceiling = mod(n,i)
  group_floor = i - group_ceiling
  count_floor2 = div(n,k); count_ceiling2 = count_floor2+1
  group_ceiling2 = mod(n,k)
  group_floor2 = k - group_ceiling2

  while true
  if (count_floor & 1)==0
    group_floor=(2*group_floor)+group_ceiling
  end

  if (count_ceiling & 1)==0
    group_ceiling=(2*group_ceiling)+group_floor
  end
  
  if count_floor < count_ceiling
    count_floor=div(count_floor,2) 
    count_ceiling=cld(count_ceiling,2)
  else
    count_floor=cld(count_floor,2)
    count_ceiling=div(count_ceiling,2)
  end

  if (count_floor2 & 1)==0
    group_floor2=(2*group_floor2)+group_ceiling2
  end
  if (count_ceiling2 & 1)==0
    group_ceiling2=(2*group_ceiling2)+group_floor2
  end
  if count_floor2 < count_ceiling2
    count_floor2=div(count_floor2,2)
    count_ceiling2=cld(count_ceiling2,2)
  else
    count_floor2=cld(count_floor2,2)
    count_ceiling2=div(count_ceiling2,2)
  end

if ((count_floor < group_ceiling) | (count_ceiling < group_floor)) & ((count_floor2 < group_ceiling2) | (count_ceiling2 < group_floor2))
  break
end

  if(group_ceiling>1 && group_floor>1)
    if(div(count_floor,group_ceiling)*group_ceiling==count_floor)
        return (group_ceiling , div(n,group_ceiling))
        end
        if(div(count_ceiling,group_floor)*group_floor==count_ceiling)
          return (group_floor , div(n,group_floor))
          end
  end

  if(group_ceiling2>1 && group_floor2>1)
    if(div(count_floor2,group_ceiling2)*group_ceiling2==count_floor2)
        return (group_ceiling2 , div(n,group_ceiling2))
        end
        if(div(count_ceiling2,group_floor2)*group_floor2==count_ceiling2)
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
    return("PRIME")
  end
end #outer while for i, k values
end
```
