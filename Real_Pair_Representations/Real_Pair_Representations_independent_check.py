"""Independent Python arithmetic cross-check of the formulas in the R paper.
This does NOT execute the companion base-R test suite.
"""
import cmath
import math
import random


def phi(z):
    return (z.real - z.imag, z.real + z.imag)


def phi_inv(p):
    return complex((p[0] + p[1]) / 2, (p[1] - p[0]) / 2)


def multiply(p, q):
    s = (p[0] * q[1] + p[1] * q[0]) / 2
    t = (p[1] * q[1] - p[0] * q[0]) / 2
    return (s - t, s + t)


def inverse(p):
    m = max(map(abs, p))
    if m == 0:
        raise ZeroDivisionError
    u = tuple(v / m for v in p)
    scale = 2 / (m * sum(v * v for v in u))
    return (scale * u[1], scale * u[0])


def psi(q):
    a, b, c, d = q
    return (a - b, a + b, c - d, c + d)


def qmul(p, q):
    p1, p2 = p[:2], p[2:]
    q1, q2 = q[:2], q[2:]
    a = multiply(p1, q1)
    b = multiply(p2, q2[::-1])
    c = multiply(p1, q2)
    d = multiply(p2, q1[::-1])
    return (a[0] - b[0], a[1] - b[1], c[0] + d[0], c[1] + d[1])


def qconj(p):
    return (p[1], p[0], -p[2], -p[3])


def qinverse(p):
    m = max(map(abs, p))
    if m == 0:
        raise ZeroDivisionError
    u = tuple(v / m for v in p)
    scale = 2 / (m * sum(v * v for v in u))
    return tuple(scale * v for v in qconj(u))


def qnative(v, w):
    a, b, c, d = v
    e, f, g, h = w
    return (a*e - b*f - c*g - d*h,
            a*f + b*e + c*h - d*g,
            a*g - b*h + c*e + d*f,
            a*h + b*g - c*f + d*e)


def qexp(p):
    a = (p[0] + p[1]) / 2
    b = (p[1] - p[0]) / 2
    c = (p[2] + p[3]) / 2
    d = (p[3] - p[2]) / 2
    r = math.sqrt(b*b+c*c+d*d)
    if r < 1e-4:
        r2 = r*r
        sinc = 1-r2/6+r2*r2/120
    else:
        sinc = math.sin(r)/r
    return psi((math.exp(a)*math.cos(r),
                math.exp(a)*sinc*b,
                math.exp(a)*sinc*c,
                math.exp(a)*sinc*d))


def close(a, b, tol=1e-10):
    return all(abs(x-y) <= tol*(1+abs(y)) for x,y in zip(a,b))


ONE = psi((1,0,0,0))
I, J, K = (psi(q) for q in ((0,1,0,0),(0,0,1,0),(0,0,0,1)))
assert close(qmul(I,J),K)
assert close(qmul(J,I),tuple(-x for x in K))
assert close(qinverse(J),tuple(-x for x in J))
assert close(qinverse(K),tuple(-x for x in K))
assert close(qmul(J, qinverse(J)),ONE)
assert close(qmul(K, qinverse(I)),tuple(-x for x in J))
assert close(qexp((0,0,0,0)),ONE)
assert close(qexp(tuple(math.pi*x for x in J)),tuple(-x for x in ONE))
assert close(qexp(tuple(1e-10*x for x in J)),psi((math.cos(1e-10),0,math.sin(1e-10),0)))
for z in (complex(0.5,0.3),complex(-3.1,0),complex(2,-7)):
    assert abs(phi_inv(phi(z))-z) < 1e-12
    assert close(multiply(phi(z), phi(2+3j)), phi(z*(2+3j)))
    assert close(multiply(phi(z),inverse(phi(2+3j))),phi(z/(2+3j)))
    assert close((math.exp(z.real)*(math.cos(z.imag)-math.sin(z.imag)), math.exp(z.real)*(math.cos(z.imag)+math.sin(z.imag))), phi(cmath.exp(z)))

rng=random.Random(20261104)
for _ in range(10000):
    v=[rng.gauss(0,1) for _ in range(4)]
    w=[rng.gauss(0,1) for _ in range(4)]
    p,q=psi(v),psi(w)
    assert close(qmul(p,q),psi(qnative(v,w)))
    assert close(qmul(p,qinverse(p)),ONE)
    assert close(qmul(qinverse(p),p),ONE)
    assert close(qconj(p),psi((v[0],-v[1],-v[2],-v[3])))
    assert close((sum(x*x for x in p)/2,), (sum(x*x for x in v),))
    assert close((sum(x*x for x in qmul(p,q))/2,),
                 ((sum(x*x for x in v))*(sum(x*x for x in w)),))

print('PASS: complex identities, quaternion basis, inverse identities,')
print('      near-zero exponential, and 10000 randomized quaternion cases.')
print('NOTE: Python cross-check only; companion R code requires Rscript.')
