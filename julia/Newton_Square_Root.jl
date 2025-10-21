function Newton_sqrt{T<:Union{Int64,UInt64,Int128,UInt128,BigInt}}(N::T)
    a = 2
    b = N
    while abs(a-b) > 1
        b = div(N,a)
        a = div((a + b), 2)
    end
    return a
end
