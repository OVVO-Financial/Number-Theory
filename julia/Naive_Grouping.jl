function NG(n)
r = Newton_sqrt(n) + 1
limit = div(r,2)
  if (limit & 1)==0
    limit = limit - 1
  end
rr = div(limit,2)
  if (rr & 1)==0
    rr = rr -1
  end
i = rr
k = i+2

while true
  count_floor = div(n,i)
  if count_floor*i==n return (i, count_floor) end
  count_ceiling = count_floor+1
  group_ceiling = mod(n,i)
  group_floor = i - group_ceiling
  count_floor2 = div(n,k)
  if count_floor2*k==n return (k, count_floor) end
  count_ceiling2 = count_floor2+1
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

  if ((count_floor < group_ceiling+1) | (count_ceiling < group_floor+1))
      break
  end

  if (count_floor2 & 1)==0
    group_floor2=(2*group_floor2)+group_ceiling2
  end
  if (count_ceiling2 & 1)==0
    group_ceiling2=(2*group_ceiling2)+group_floor2
  end
  if count_floor2< count_ceiling2
    count_floor2=div(count_floor2,2)
    count_ceiling2=cld(count_ceiling2,2)
  else
    count_floor2=cld(count_floor2,2)
    count_ceiling2=div(count_ceiling2,2)
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
    return ("PRIME")
  end
end #outer while for i, k values
end
