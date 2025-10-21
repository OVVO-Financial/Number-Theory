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
  summand_floor = div(n,i)
  if summand_floor*i==n return (i, summand_floor) end
  summand_ceiling = summand_floor+1
  group_ceiling = mod(n,i)
  group_floor = i - group_ceiling
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
