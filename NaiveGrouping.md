# Naive Grouping
The naive grouping routine is as stated, *a naive grouping method for factorization.*
It is deterministic and can be distributed along any number of "seed groups".
One of 3 known seed groups resides below `0.5*sqrt(n)`.

The method is best described with an example:

1. For `n = 798,607` if we start with a seed group of **59**, we would have **17** groups of **13,535** and **42** groups of **13,536**.
    > * The group counts (e.g., 13,535 and 13,536 generated from `798607/59 = 13535.71`) are to have a maximum difference of 1 (hence naive grouping!).  
2. Neither the **17** nor the **42** will evenly parse the other grouping counts so we double the seed group and halve the counts.
    > * `mod(13535,42) > 1` and `mod(13536,17) > 1`.  
3. Now we have **118** groups in total, **17** groups of **6,767** and **101** groups of **6,768**. 
4. The **6,767** can be distributed among the **101** groups, resulting in factors of **101** groups of **7,907**.
    > * `mod(6767,101) == 0`. 

There are minimal division operations (one per seed group) performed, rather elementary doubling and halving rules resulting in significantly reduced `mod` checks versus the original number `n`.
``` julia
@time NG(798607)
1×2 Array{Float64,2}:
 101.0  7907.0
0.000044 seconds (7 allocations: 288 bytes)
```

``` julia
@time NG(978508015703)
1×2 Array{Float64,2}:
 752867.0  1.29971e6
0.030094 seconds (7 allocations: 288 bytes)
```

Below is the `julia` code:
``` julia
function NG(n)
r = ceil(sqrt(n))
limit = floor(r/2)
  if mod(limit,2)==0
    limit = limit - 1
  end
rr = floor(limit/2)
  if mod(rr,2)==0
    rr = rr -1
  end
i = rr
k = i+2

while true
  z= n/i
  z2= n/k
  count_ceiling = ceil(z)
  count_floor = floor(z)
  group_ceiling = mod(n,i)
  group_floor = i - group_ceiling
  count_ceiling2 = ceil(z2)
  count_floor2 = floor(z2)
  group_ceiling2 = mod(n,k)
  group_floor2 = k - group_ceiling2

  while true
      if count_ceiling==count_floor
        return ([i n/i])
        end
      if count_ceiling2==count_floor2
        return ([k n/k])
        end

  if mod(count_floor,2)==0
    group_floor=(2*group_floor)+group_ceiling
  end
  if mod(count_ceiling,2)==0
    group_ceiling=(2*group_ceiling)+group_floor
  end
  if count_floor < count_ceiling
    count_floor=floor(count_floor/2)
    count_ceiling=ceil(count_ceiling/2)
  else
    count_floor=ceil(count_floor/2)
    count_ceiling=floor(count_ceiling/2)
  end

  if mod(count_floor2,2)==0
    group_floor2=(2*group_floor2)+group_ceiling2
  end
  if mod(count_ceiling2,2)==0
    group_ceiling2=(2*group_ceiling2)+group_floor2
  end
  if count_floor2< count_ceiling2
    count_floor2=floor(count_floor2/2)
    count_ceiling2=ceil(count_ceiling2/2)
  else
    count_floor2=ceil(count_floor2/2)
    count_ceiling2=floor(count_ceiling2/2)
  end

if ((count_floor < group_ceiling) | (count_ceiling < group_floor)) &((count_floor2 < group_ceiling2) | (count_ceiling2 < group_floor2))
  break
end

  if(group_ceiling>1 && group_floor>1)
  if(group_ceiling>1)
    if(mod(count_floor,group_ceiling)==0)
        return [group_ceiling n/group_ceiling]
        end
        end
  if(group_floor>1)
        if(mod(count_ceiling,group_floor)==0)
          return [group_floor n/group_floor]
          end
          end
  end

  if(group_ceiling2>1 && group_floor2>1)
    if(mod(count_floor2,group_ceiling2)==0)
      return [group_ceiling2 n/group_ceiling2]
      end
      if(mod(count_ceiling2,group_floor2)==0)
        return [group_floor2 n/group_floor2]
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
