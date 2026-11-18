# Integration and verification status

## Implemented families

- Shared GMP product-tree, VF, CRT/Fermat residue search, and regional workers (`vf-factor-lab/`).
- Deterministic naive grouping and odd-seed chains (`naive_grouping_gmp/`).
- Exact integer Newton/convex region certificates and one/two-sided smooth-relation generators (`vertical_relation_lab/`).
- Direct `Y=2ab` Fermat wheel and batch global factor scans (`global_vertical_locator/`).
- Three genuine nontrivial witnesses on `H^2-Y^2=N^2`, monotone integer staircase, cube-root cutoff (`vertical_witness_bridge/`).
- Iterated-average dyadic GCD progressions and signed witnesses (`iterated_average_vertical_bridge/`).
- N-only residue covers, modular Newton squaring, Gaussian norm-one group factoring, independent Gaussian norm residual, and bounded hybrid controls (`dyadic_global_attack/`).

## Boundaries of claims

These are executable experiments with mathematical proofs for stated identities and
finite checks for selected numerical properties. There is no demonstrated
random/general RSA-scale factorization advantage over QS/GNFS, no CUDA timing
for these new experiments, and no improved general factoring complexity bound.
The original GNFS-like source is retained under `experimental/` and is **not**
a validated GNFS implementation. The `N mod 10 == 5` semiprime case is
short-circuited by a factor-5 precheck in the updated original C++ engine.

## Reproduction

```bash
cd 'Prime Factorization/hybrid_research'
make -C naive_grouping_gmp check
make -C global_vertical_locator test
make -C vertical_relation_lab test
make -C iterated_average_vertical_bridge test
make -C dyadic_global_attack test
make -C vf-factor-lab -f Makefile.gmp selftest
```

The `dyadic_global_attack` sample-set tests consume included JSON oracle files;
these oracles are for verification and are never passed as hidden factor inputs to
any C++ factoring command.
