# Additional one-sided relation benchmark

`relation_bench` was extended with three N-only candidate-ordering mechanisms which **all** generate the same valid one-sided relation `x² == x²-N (mod N)`. No second independent smooth side is required. These methods are tested with the same prime bound B=1021 and raw maximum K=32768; their unique-x coverage differs and is reported.

For the 40-bit example `N=978508015703`, `--stop-on-factor`, median of 3 CPU kernel runs:

| candidate method | candidates evaluated | unique x | smooth, valid QS relations | extracted factor | median time |
|---|---:|---:|---:|---:|---:|
| QS consecutive x | 2,621 | 2,621 | 75 | 752867 | 5.19 ms |
| CTMWalkOneSide | 32,768 | 628 | 14 | none | 8.24 ms |
| GroupingBandOneSide | 18,666 | 2,621 | 75 | 752867 | 12.88 ms |
| VerticalOneSide | 6,381 | 2,621 | 75 | 752867 | 8.13 ms |

See `results/one_sided_benchmark.json` for all five inputs and the actual median samples. The CTM duplicate collapse is essential: 32,768 states represent only 628 distinct midpoint trials in this run. The grouping and vertical schedules eventually reached the same 2,621 QS midpoints at higher arithmetic/setup cost. **The new one-sided mechanism is correct but did not outperform QS in this experiment.** The experiment is not a GNFS performance claim.

## Exact prime-power valuation mass over intervals

The additional theorem and Python implementation count `sum_R v_p(R²-N)` via prime-power roots and arithmetic-progression floors, then certify zero B-smooth values when `product p^{mass_p}<min(R²-N)`. This uses only N and the interval endpoints. An independent oracle passed **504** window/multiplier cases, and the mass certificate proved **276** of those windows empty of B-smooth QS relations. The fixed-window performance experiment (`results/prime_power_region_benchmark.json`) used 120 windows per width: at B=7 and width=16, the bound certified 119 of 120 empty; at B=101 and width=16 it certified **none** despite 119 being actually empty. Certificate computation was also slower than directly checking these short windows. This is a legitimate interval inference, but is not yet a competitive relation-generation accelerator.
