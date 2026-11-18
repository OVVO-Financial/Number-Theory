# Dyadic iterated-average / vertical-green-line bridge

Standalone, exact add-on to the `Number-Theory` hybrid research PR-ready package.

Files:
- `DYADIC_VERTICAL_THEORY.md`: proofs of the GCD-equivalent dyadic residual, `q mod 2^k` multiplier channels, real monotonic error on the original `Y=2ab` line, and all three integer Pythagorean witness lifts.
- `dyadic_scan.cpp`: GNU GMP C++17 bounded arithmetic-progression batch-GCD worker; N-only input; no factor oracle; exact verification. Does not assert primality on exhaustion.
- `test_dyadic_vertical.py`: 5,000-semiprime exact arithmetic tests and actual worker controls.
- `benchmark_retrospective.py`: 600-case factor-informed rank census, explicitly **not** a real-time factor benchmark.
- `results.json` and `runtime_batch.json`: raw results and host-local kernel sample timings.

Build and run:

```
make
make test
python3 benchmark_retrospective.py
./dyadic_scan 978508015703 8 12686 256
```

`dyadic_scan N level_k max_offset [batch_size=128]` searches the plus/minus integer offsets from the selected rounded iterated average. It uses only N and known arithmetic, not the true p or q. One can distribute levels and bounded offset blocks, but parallel speedup cannot by itself beat the global arithmetic search horizon.

**No claim of a global region-locator breakthrough, subexponential complexity, or GNFS competitiveness.**
