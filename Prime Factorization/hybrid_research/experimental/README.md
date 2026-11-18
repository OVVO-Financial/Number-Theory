# Experimental GNFS-like lattice code (not validated)

`gnfs_gmp_unverified.cpp` is included only for reproducibility. It is **not** a correct or benchmarked General Number Field Sieve. With `alpha^2=-c`, the algebraic norm of `a+b*alpha` is `a^2+c*b^2`; the uploaded code substitutes `a^2+b^2` in its sieving stage and has additional matrix and algebraic-square shortcomings. Do not use it as a factorization baseline until repaired.
