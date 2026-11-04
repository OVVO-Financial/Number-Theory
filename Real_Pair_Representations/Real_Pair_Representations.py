"""Real-pair complex and real-quadruple quaternion arithmetic.

Companion to Fred Viole, *Real-Pair Representations of Complex Numbers
and Quaternions*, November 2026.

The arithmetic operators use only real floating-point tuples. ``phi`` and
``phi_inv`` are convenience adapters to Python's built-in complex type.
The quaternion convention is q = a + b*i + c*j + d*k with i*j = k.
Quaternion division is RIGHT division: p / q := p * inverse(q).

Requires Python 3.9+; no third-party packages.
"""

import math
from typing import Sequence, Tuple

Pair = Tuple[float, float]
Quat = Tuple[float, float, float, float]


def _values(v: Sequence[float], length: int, label: str) -> tuple:
    """Validate a finite real tuple and normalize numeric types."""
    try:
        out = tuple(float(x) for x in v)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a finite real vector") from exc
    if len(out) != length or not all(math.isfinite(x) for x in out):
        raise ValueError(f"{label} must be a finite real vector of length {length}")
    return out


def phi(z: complex) -> Pair:
    """a + b*i -> (a-b, a+b)."""
    z = complex(z)
    if not (math.isfinite(z.real) and math.isfinite(z.imag)):
        raise ValueError("z must be a finite complex scalar")
    return (z.real - z.imag, z.real + z.imag)


def phi_inv(p: Sequence[float]) -> complex:
    """(x,y) -> (x+y)/2 + i*(y-x)/2."""
    x, y = _values(p, 2, "p")
    return complex((x + y) / 2, (y - x) / 2)


def real_pair_add(p: Sequence[float], q: Sequence[float]) -> Pair:
    x, y = _values(p, 2, "p")
    u, v = _values(q, 2, "q")
    return (x + u, y + v)


def real_pair_subtract(p: Sequence[float], q: Sequence[float]) -> Pair:
    x, y = _values(p, 2, "p")
    u, v = _values(q, 2, "q")
    return (x - u, y - v)


def real_pair_multiply(p: Sequence[float], q: Sequence[float]) -> Pair:
    """Transport of complex multiplication, NOT componentwise product."""
    x, y = _values(p, 2, "p")
    u, v = _values(q, 2, "q")
    s = (x * v + y * u) / 2
    t = (y * v - x * u) / 2
    return (s - t, s + t)


def real_pair_conjugate(p: Sequence[float]) -> Pair:
    x, y = _values(p, 2, "p")
    return (y, x)


def real_pair_norm_sq(p: Sequence[float]) -> float:
    x, y = _values(p, 2, "p")
    return (x * x + y * y) / 2


def real_pair_inverse(p: Sequence[float]) -> Pair:
    x, y = _values(p, 2, "p")
    m = max(abs(x), abs(y))
    if m == 0:
        raise ZeroDivisionError("zero complex value has no inverse")
    u, v = x / m, y / m
    scale = (2.0 / (u * u + v * v)) / m
    result = (scale * v, scale * u)
    if not all(math.isfinite(t) for t in result):
        raise OverflowError("complex inverse is outside finite float range")
    return result


def real_pair_divide(p: Sequence[float], q: Sequence[float]) -> Pair:
    return real_pair_multiply(p, real_pair_inverse(q))


def real_pair_exp(p: Sequence[float]) -> Pair:
    x, y = _values(p, 2, "p")
    a, b = (x + y) / 2, (y - x) / 2
    t = math.exp(a)
    return (t * (math.cos(b) - math.sin(b)),
            t * (math.cos(b) + math.sin(b)))


def psi(q: Sequence[float]) -> Quat:
    """Quaternion coefficients (a,b,c,d) -> two real pairs."""
    a, b, c, d = _values(q, 4, "q")
    return (a - b, a + b, c - d, c + d)


def psi_inv(p: Sequence[float]) -> Quat:
    x, y, u, v = _values(p, 4, "p")
    return ((x + y) / 2, (y - x) / 2,
            (u + v) / 2, (v - u) / 2)


def quat_add(p: Sequence[float], q: Sequence[float]) -> Quat:
    a = _values(p, 4, "p")
    b = _values(q, 4, "q")
    return tuple(x + y for x, y in zip(a, b))


def quat_subtract(p: Sequence[float], q: Sequence[float]) -> Quat:
    a = _values(p, 4, "p")
    b = _values(q, 4, "q")
    return tuple(x - y for x, y in zip(a, b))


def quat_multiply(p: Sequence[float], q: Sequence[float]) -> Quat:
    """Cayley-Dickson product transported through the pair mapping."""
    a = _values(p, 4, "p")
    b = _values(q, 4, "q")
    p1, p2, q1, q2 = a[:2], a[2:], b[:2], b[2:]
    aa = real_pair_multiply(p1, q1)
    bb = real_pair_multiply(p2, real_pair_conjugate(q2))
    cc = real_pair_multiply(p1, q2)
    dd = real_pair_multiply(p2, real_pair_conjugate(q1))
    return (aa[0] - bb[0], aa[1] - bb[1],
            cc[0] + dd[0], cc[1] + dd[1])


def quat_conjugate(p: Sequence[float]) -> Quat:
    """Quaternion conjugation: swap first pair, NEGATE second pair."""
    x, y, u, v = _values(p, 4, "p")
    return (y, x, -u, -v)


def quat_norm_sq(p: Sequence[float]) -> float:
    x, y, u, v = _values(p, 4, "p")
    return (x * x + y * y + u * u + v * v) / 2


def quat_inverse(p: Sequence[float]) -> Quat:
    a = _values(p, 4, "p")
    m = max(abs(x) for x in a)
    if m == 0:
        raise ZeroDivisionError("zero quaternion has no inverse")
    u = tuple(x / m for x in a)
    scale = (2.0 / sum(x * x for x in u)) / m
    result = tuple(scale * x for x in quat_conjugate(u))
    if not all(math.isfinite(x) for x in result):
        raise OverflowError("quaternion inverse is outside finite float range")
    return result


def quat_divide(p: Sequence[float], q: Sequence[float]) -> Quat:
    """Right division p/q := p * inverse(q). Order matters."""
    return quat_multiply(p, quat_inverse(q))


def quat_exp(p: Sequence[float]) -> Quat:
    """Quaternion exponential using sinc(r) near the real axis."""
    a, b, c, d = psi_inv(p)
    r = math.hypot(b, c, d)
    if r < 1e-4:
        r2 = r * r
        sinc = 1 - r2 / 6 + r2 * r2 / 120
    else:
        sinc = math.sin(r) / r
    ea = math.exp(a)
    return psi((ea * math.cos(r), ea * sinc * b,
                ea * sinc * c, ea * sinc * d))
