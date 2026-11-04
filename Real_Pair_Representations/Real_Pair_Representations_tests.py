"""Deterministic and randomized tests of the companion Python implementation.

Run: python3 Real_Pair_Representations_tests.py
No third-party packages required.
"""

import cmath
import math
import random
import unittest

from Real_Pair_Representations import (
    phi, phi_inv, psi, psi_inv,
    real_pair_add, real_pair_subtract, real_pair_multiply,
    real_pair_conjugate, real_pair_norm_sq, real_pair_inverse,
    real_pair_divide, real_pair_exp,
    quat_add, quat_subtract, quat_multiply, quat_conjugate,
    quat_norm_sq, quat_inverse, quat_divide, quat_exp,
)


def native_quaternion_multiply(v, w):
    """Independent coefficient-by-coefficient reference formula."""
    a, b, c, d = v
    e, f, g, h = w
    return (a*e - b*f - c*g - d*h,
            a*f + b*e + c*h - d*g,
            a*g - b*h + c*e + d*f,
            a*h + b*g - c*f + d*e)


class TestRealPair(unittest.TestCase):
    def assertClose(self, actual, expected, tolerance=1e-10):
        self.assertEqual(len(actual), len(expected))
        for x, y in zip(actual, expected):
            self.assertTrue(math.isclose(x, y, rel_tol=tolerance,
                                          abs_tol=tolerance), (x, y))

    def test_complex_round_trip(self):
        for z in (0j, 1+0j, 2+3j, -5-8j, 0.5+1e-10j):
            self.assertClose((phi_inv(phi(z)).real, phi_inv(phi(z)).imag),
                             (z.real, z.imag))

    def test_complex_operators(self):
        z, w = 2+3j, -4+2j
        p, q = phi(z), phi(w)
        self.assertClose(real_pair_add(p, q), phi(z+w))
        self.assertClose(real_pair_subtract(p, q), phi(z-w))
        self.assertClose(real_pair_multiply(p, q), phi(z*w))
        self.assertClose(real_pair_conjugate(p), phi(z.conjugate()))
        self.assertClose((real_pair_norm_sq(p),), (abs(z)**2,))
        self.assertClose(real_pair_inverse(q), phi(1/w))
        self.assertClose(real_pair_divide(p, q), phi(z/w))
        self.assertClose(real_pair_exp(p), phi(cmath.exp(z)))

    def test_complex_zero_inverse(self):
        with self.assertRaises(ZeroDivisionError):
            real_pair_inverse((0, 0))

    def test_quaternion_basis(self):
        one = psi((1, 0, 0, 0))
        i = psi((0, 1, 0, 0))
        j = psi((0, 0, 1, 0))
        k = psi((0, 0, 0, 1))
        self.assertClose(quat_multiply(i, j), k)
        self.assertClose(quat_multiply(j, i), tuple(-x for x in k))
        self.assertClose(quat_multiply(j, j), tuple(-x for x in one))
        self.assertClose(quat_inverse(j), tuple(-x for x in j))
        self.assertClose(quat_inverse(k), tuple(-x for x in k))
        self.assertClose(quat_multiply(j, quat_inverse(j)), one)
        self.assertClose(quat_divide(k, i), tuple(-x for x in j))

    def test_quaternion_add_subtract(self):
        p = psi((1, 2, 3, 4))
        q = psi((-2, 0, -1, 8))
        self.assertClose(quat_add(p, q), psi((-1, 2, 2, 12)))
        self.assertClose(quat_subtract(p, q), psi((3, 2, 4, -4)))

    def test_quaternion_exp(self):
        one = psi((1, 0, 0, 0))
        j = psi((0, 0, 1, 0))
        zero = psi((0, 0, 0, 0))
        self.assertClose(quat_exp(zero), one)
        self.assertClose(quat_exp(tuple(math.pi*x for x in j)),
                         tuple(-x for x in one))
        r = 1e-10
        self.assertClose(quat_exp(tuple(r*x for x in j)),
                         psi((math.cos(r), 0, math.sin(r), 0)))

    def test_quaternion_zero_inverse(self):
        with self.assertRaises(ZeroDivisionError):
            quat_inverse(psi((0, 0, 0, 0)))

    def test_small_large_inverse(self):
        one = psi((1, 0, 0, 0))
        for scale in (1e-200, 1e200):
            p = psi(tuple(scale * x for x in (1, 2, -3, 4)))
            self.assertClose(quat_multiply(p, quat_inverse(p)), one)
            self.assertClose(quat_multiply(quat_inverse(p), p), one)
            c = tuple(scale*x for x in phi(1+2j))
            self.assertClose(real_pair_multiply(c, real_pair_inverse(c)),
                             (1, 1))

    def test_gaussian_integers_and_fermat(self):
        for a, b in ((3, 5), (-2, 7), (0, -4)):
            x, y = phi(complex(a, b))
            self.assertEqual((int(x) - int(y)) % 2, 0)
            self.assertEqual((x*x+y*y)/2, a*a+b*b)
        for c, q in ((3, 5), (7, 9), (0, 12), (-3, 5)):
            z = complex((c+q)/2, (q-c)/2)
            self.assertClose(phi(z), (c, q))
            self.assertClose(((z*z).real, (z*z).imag),
                             (c*q, (q*q-c*c)/2))

    def test_zeta_reflection_coordinates(self):
        rho = 0.75+8j
        u, v = phi(rho)
        self.assertClose(phi(1-rho.conjugate()), (1-v, 1-u))
        self.assertClose(((u+v)/2,), (rho.real,))

    def test_randomized_quaternion_properties(self):
        rng = random.Random(20261104)
        one = psi((1, 0, 0, 0))
        for _ in range(10000):
            v = tuple(rng.gauss(0, 1) for _ in range(4))
            w = tuple(rng.gauss(0, 1) for _ in range(4))
            p, q = psi(v), psi(w)
            self.assertClose(psi_inv(p), v)
            self.assertClose(quat_multiply(p, q),
                             psi(native_quaternion_multiply(v, w)))
            self.assertClose(quat_multiply(p, quat_inverse(p)), one)
            self.assertClose(quat_multiply(quat_inverse(p), p), one)
            self.assertClose(quat_conjugate(p),
                             psi((v[0], -v[1], -v[2], -v[3])))
            self.assertClose((quat_norm_sq(p),), (sum(x*x for x in v),))
            self.assertClose((quat_norm_sq(quat_multiply(p, q)),),
                             (sum(x*x for x in v)*sum(x*x for x in w),))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            psi((1, 2, 3))
        with self.assertRaises(ValueError):
            quat_inverse((1, 2, 3, float("nan")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
