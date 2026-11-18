# Squared complex-plane Newton lattice and relation-generation experiment

**Scope:** A mathematically exact Fermat/CRT factor-point region sieve and GMP Newton-lattice certificate, followed by a controlled smooth-relation comparison. **Not an implemented GNFS**, and not a GNFS-scale performance claim. Original source: `Complex space.md`; mathematical proofs in `math/FORMALIZATION.md`.

## Verified insight: the vertical green line

For an exact factor point, $z=a+bi$, $a^2-b^2=N$, $Y=2ab$, $H=a^2+b^2$, so $H^2=N^2+Y^2$, $(H-Y)(H+Y)=N^2$, and $H-Y=p^2$, $H+Y=q^2$. The curve $H(Y)=\sqrt{N^2+Y^2}$ is increasing and convex; its **rounded integer residual is not monotone**. For any off-curve CTM point with $E=a^2-b^2-N$, the squared residual is exactly $E(2N+E)$; a smooth $E$ does not create a new smooth square relation by itself.

The factor-point wheel uses unknown factor residues $p$ modulo $2m$ only as **a finite modular enumeration**, not as factor knowledge. For $m=420$ and $N=978508015703$, only six paired $(Y,H)$ residues modulo 420 are possible. The exact rational Newton bounds and Euclidean floor sum are proved in the mathematical appendix. The wheel **must not** be applied as an exclusion sieve to general nonzero-residual relations.

## Exact regional geometry, 80 deterministic windows per case

| Input | Bits | Empty / 80 with integer ceilings | Empty / 80 after 1 Newton step | Surviving integer envelope points after Newton |
|---|---:|---:|---:|---:|
| 40b | 40 | 8 | 79 | 1 |
| 64b | 63 | 10 | 80 | 0 |
| 128b | 128 | 9 | 80 | 0 |

Windows were 100,000 $Y$ units wide, distributed by an $N$-only fixed schedule on $0\le Y\le 3N/4$ (valid as a full factor-domain bound only when $p<q<2p$). The observed *factor-bearing* 128-bit validation window retains its exact integer point: 477 candidate lattice points under integer ceilings, **1 after Newton**. Local exclusion is not a full-domain speedup: there can be an immense number of separate windows.

## Fixed-work relation generation: B = 1021, K = 32768

All modes use the same GMP prime-base smoothness test, the required one- or two-sided exponent matrix, full sign tracking, duplicate-residual suppression, exact mod-$N$ congruences and gcd verification. The factor inputs were withheld from the generator binaries. `CTM` is the nearest-boundary Fermat projection; `CTMWalk` is a 4-lane increment-only CTM staircase. `GroupingBand` uses only the specified balance constraint. `SquareY` samples the vertical line uniformly in the balanced range.

| Actual bits | QS valid / deps | CTM valid / deps | CTMWalk valid / deps | Grouping-band valid / deps | Vertical Y valid |
|---:|---:|---:|---:|---:|---:|
| 28 | 1021 / 939 | 1750 / 1408 | 1661 / 1322 | 198 / 0 | 13 |
| 32 | 446 / 363 | 1096 / 756 | 1046 / 715 | 309 / 10 | 0 |
| 36 | 1117 / 1030 | 675 / 351 | 665 / 337 | 409 / 90 | 0 |
| 40 | 370 / 283 | 378 / 58 | 311 / 3 | 327 / 6 | 0 |
| 44 | 209 / 127 | 207 / 0 | 138 / 0 | 214 / 0 | 0 |
| 47 | 221 / 138 | 121 / 0 | 51 / 0 | 114 / 0 | 0 |
| 55 | 127 / 43 | 25 / 0 | 4 / 0 | 20 / 0 | 0 |
| 63 | 23 / 0 | 3 / 0 | 0 / 0 | 0 / 0 | 0 |

In the 40-bit manuscript example, CTM produced **11,041 distinct smooth residual hits**, but only **378 valid two-sided smooth relations** (rank 320). `CTMWalk` produced 311 valid relations, grouping-band 327, and QS 370. Unlike the others, QS starts with a square on one side; counting raw smooth residuals as QS-quality relations would be misleading.

## Time to first independently verified factor under the same maximum candidate budget

Seven of eight semiprimes were resolved by the simple QS relation generator, four by each CTM strategy, three by the balanced grouping-band approach, none by direct vertical-Y sampling. `—` means no proper factor within 32,768 candidates. Timings are three-run medians of the *C++ computation* including candidate construction, smooth factoring, matrix/dependency work and gcd verification, excluding process start.

| Actual bits | QS | CTM projection | CTM walk | Grouping band |
|---:|---:|---:|---:|---:|
| 28 | 1.02 ms | 1.44 ms | 6.00 ms | — |
| 32 | 3.63 ms | 7.96 ms | 11.32 ms | 12.21 ms |
| 36 | 0.96 ms | 16.38 ms | 28.09 ms | 26.43 ms |
| 40 | 6.16 ms | 59.90 ms | 94.32 ms | 64.19 ms |
| 44 | 11.63 ms | — | — | — |
| 47 | 13.56 ms | — | — | — |
| 55 | 30.39 ms | — | — | — |
| 63 | — | — | — | — |

## Smoothness sensitivity

At 40 bits, with the same 32,768 candidates:

| Smoothness B | QS valid / deps | CTM valid / deps | CTM walk valid / deps | Grouping band valid / deps |
|---:|---:|---:|---:|---:|
| 257 | 38 / 12 | 16 / 0 | 20 / 0 | 13 / 0 |
| 1021 | 370 / 283 | 378 / 58 | 311 / 3 | 327 / 6 |
| 4093 | 1403 / 1130 | 2070 / 1040 | 2006 / 965 | 1912 / 849 |

CTM and grouping can produce genuine factors by combining smooth relations, especially at larger B. But they require both sides smooth, and at B=4093 CTM projection needed 345 ms for its full-work 40-bit sample compared with 187 ms for QS, with more columns in the linear algebra. Results do not establish a consistent efficiency gain.

## Valid relations versus smooth-looking residuals

1. **CTM:** $L=(a-b)(a+b)=N+E$, so $L\equiv E\pmod N$. Two-sided smoothness of $L$ and $E$ supports a proper parity matrix and a possible square congruence. Different CTM paths with the same $E$ give the **same** arithmetic relation and are deduplicated.
2. **Grouping:** $N=g\lfloor N/g\rfloor+r$, with $E=-r$ and $L=N+E$. Its relations collapse to the same smooth-neighbor form; multiplicity of successful seed chains does not automatically imply the same multiplicity of independent congruences.
3. **Squared vertical line:** $H^2-Y^2=N^2+W$. A generic smooth $W$ requires the left-hand product to be smooth for our two-sided relation protocol. For a CTM-derived point, $W=E(2N+E)$ and $H^2-Y^2=(N+E)^2$; modulo $N$ this is a tautological square relation, not independent data.
4. **QS baseline:** $x^2-N=f(x)$; only $f(x)$ needs to be smooth, and a dependency of its exponent vectors produces a genuine square congruence.

## Verdict and limitations

The rational-Newton vertical sieve provides a strong, rigorously certified **local** geometry filter. The separate relation experiment confirms that CTM and grouping can produce valid dependencies and actual factors, but this implementation does not outperform the conventional QS relation generator on the tested hard end, and it gives no GNFS-level scaling. The earlier uploaded `gnfs_gmp.cpp` cannot be a GNFS control until its algebraic norm and other NFS machinery are corrected; there was no CADO-NFS installation in the environment. The tests cover 28–63 bits with fixed K and factor-base bounds, not RSA-100/110 or GNFS-relevant input sizes.

The remaining mathematical opportunity is a new **one-sided square relation** from the complex/grouping coordinates, or a truly cheaper way to identify useful $Y$ intervals without visiting a square-root number of them. Those would change the relevant scaling; more smooth *nonzero* residuals alone do not.

## Files and reproducibility

`math/FORMALIZATION.md` is the proof reference. `src/vertical_lattice.cpp` and `src/relation_bench.cpp` are the standalone C++17/GMP implementations. `tests/test_vertical.py`, `tests/test_newton.py`, and `tests/test_relations.py` provide independent mathematical checks. Raw JSON with every test case, full- and early-stop metrics, smoothness sensitivity, and environment notes is in `results/`. All commands are in `README.md`.
