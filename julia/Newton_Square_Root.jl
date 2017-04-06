function Newton_sqrt(N)
    a = 1
    b = N
    while abs(a-b) > 1
        b = div(N,a)
        a = div((a + b), 2)
    end
    return a
end
