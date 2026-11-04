"""Reproducible arithmetic examples from Fred Viole's real-pair paper.
Run: python3 Real_Pair_Representations_examples.py
"""

import math
from Real_Pair_Representations import (
    phi, phi_inv, psi, real_pair_multiply, real_pair_inverse,
    quat_multiply, quat_inverse, quat_exp,
)


def main():
    print("Complex real-pair representation:")
    z = 3+4j
    p = phi(z)
    print(f"phi({z}) = {p}; phi_inv = {phi_inv(p)}")
    print(f"phi(z)*phi(z^-1) = {real_pair_multiply(p, real_pair_inverse(p))}")

    print("\nFermat factor-pair coordinates:")
    c, q = 7, 13
    z = complex((c+q)/2, (q-c)/2)
    print(f"factors ({c},{q}) -> phi(z) = {phi(z)}, Re(z^2) = {(z*z).real:g}")
    assert phi(z) == (c, q) and (z*z).real == c*q

    print("\nQuaternion identities:")
    one, i, j, k = [psi(t) for t in ((1,0,0,0), (0,1,0,0),
                                    (0,0,1,0), (0,0,0,1))]
    print(f"ij = {quat_multiply(i,j)} (expected k = {k})")
    print(f"inverse(j) = {quat_inverse(j)} (expected -j)")
    print(f"exp(pi*j) = {quat_exp(tuple(math.pi*x for x in j))} (expected -1)")
    assert quat_multiply(i, j) == k
    assert quat_inverse(j) == tuple(-x for x in j)
    assert all(abs(x-y)<1e-12 for x,y in zip(
        quat_exp(tuple(math.pi*x for x in j)), tuple(-x for x in one)))

    print("\nZeta-coordinate reflection (geometric identity only):")
    rho = 0.75 + 8j
    u, v = phi(rho)
    print(f"phi(rho) = {(u,v)}; phi(1-conj(rho)) = {phi(1-rho.conjugate())}")
    assert phi(1-rho.conjugate()) == (1-v, 1-u)
    print("All examples passed.")


if __name__ == "__main__":
    main()
