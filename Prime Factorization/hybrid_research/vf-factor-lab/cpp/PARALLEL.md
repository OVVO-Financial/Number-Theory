# Parallel regional factorization with C++ and GMP

The parallel executable combines the Fermat wheel, divisor-product search,
and optional VF probes in one local worker pool. Workers explore different
regions and return exact exclusion intervals to a shared coordinator. This
lets completed work from either search reduce the remaining work for both.

There are two meanings of “geography” here. Search-space regions are built
into the scheduler. Separate physical machines can run manually assigned
candidate bands, each with its own local worker pool. The executable does not
provide network transport, a remote queue, remote cancellation, or automatic
aggregation of results from different machines.

The implementation is in [gmp_parallel.hpp](gmp_parallel.hpp), with the CLI
in [gmp_parallel.cpp](gmp_parallel.cpp), range arithmetic in
[gmp_regions.hpp](gmp_regions.hpp), and product jobs in
[gmp_product_worker.hpp](gmp_product_worker.hpp). The prior sequential port
is documented in [README.md](README.md).

## Build and quick start

Run these commands from the project root on Linux or WSL. The build requires
GMP development headers and libraries, a C++17 compiler, and threading
support. The retained native kernel uses GCC/Clang `unsigned __int128`.

```bash
make -f Makefile.parallel
./build/vf_factor_gmp_parallel 978508015703 --threads 4 --timeout 10 --trace
./build/vf_factor_gmp_parallel --metadata
```

The equivalent direct build is:

```bash
mkdir -p build
g++ -O3 -std=c++17 -pthread \
    cpp/gmp_parallel.cpp cpp/region_native.cpp \
    -lgmpxx -lgmp -o build/vf_factor_gmp_parallel
```

`CXX`, `CPPFLAGS`, `CXXFLAGS`, and `LDFLAGS` can be supplied to Make.
`region_native.cpp` includes the existing `vf_factor/fusion_kernel.cpp`
and adds the seek operation. Do not also compile the original kernel into
this executable: that would define its symbols twice.

The binary has no Python runtime dependency. Decimal inputs and arithmetic
values in its JSON output remain exact; large integers are serialized as
strings.

## CLI controls

```bash
./build/vf_factor_gmp_parallel N [options]
./build/vf_factor_gmp_parallel --batch [options] < inputs.txt
```

| Option | Meaning and default |
|---|---|
| `--threads N` | Local worker count; 4, accepted range 1–256 |
| `--timeout SECONDS` | One cooperative wall-clock deadline; 10 seconds |
| `--max-steps U64` | One shared represented-work allowance across all workers |
| `--full` | Recursively split cofactors under the same overall budgets |
| `--trace` | Include job ranges, results, work, and exclusion certificates |
| `--factor-lo INTEGER` | Inclusive lower boundary of the candidate band |
| `--factor-hi INTEGER` | Inclusive upper boundary of the candidate band |
| `--degree N` | Product block degree; 128, accepted range 1–65536 |
| `--fermat-span U64` | Raw midpoint positions per Fermat tile; 655360 |
| `--job-seconds SECONDS` | Cooperative time quantum for a job; 0.025 seconds |
| `--fermat-ratio NUMBER` | Target Fermat/product job-time ratio; 2 |
| `--initial-fermat-tests U64` | Serial initial admissible midpoint tests; 64 |
| `--no-vf` | Disable auxiliary VF progression probes |
| `--no-fermat` | Search product regions only |
| `--no-sharing` | Use separate polynomial constructions for different steps |
| `--gmp-fermat` | Force the GMP cursor for inputs below $2^{64}$ too |
| `--no-probable-prime-check` | Continue exact search on unchecked large primes |
| `--metadata` | Report compiler, GMP, and search-policy metadata |
| `--help` | Print the command reference |

The regional scheduler uses its own fixed `--degree`. The sequential port's
adaptive-degree and `--no-feedback` CLI options are not parallel CLI options.
Exact shared coverage is part of the regional scheduler.

Batch mode uses an independent budget for each line. Invalid lines emit an
error object, subsequent lines still run, and the final exit code is nonzero
if any line is invalid. Full mode processes pending cofactors one at a time,
using the worker pool for each search, and does not start another cofactor
after the common deadline. `--full` cannot be combined with factor-band
options.

## Which regions workers search

After the common prechecks, the possible smaller-divisor domain is

$$
[L,U]=[\max(7,\text{factor-lo}),
        \min(\lfloor\sqrt N\rfloor,\text{factor-hi})],
$$

with omitted boundaries defaulting to 7 and $\lfloor\sqrt N\rfloor$.
The prechecks already handle factors 2, 3, and 5. A requested band with no
intersection with this domain is rejected when a search is needed.

| Search | Job coordinate | What a failed completed job certifies |
|---|---|---|
| Fermat wheel | Integer midpoint $R$ in a closed interval | A converted interval of possible smaller divisors |
| Divisor products | Integer candidate $d$ in a closed interval, with the modulo-30 wheel | No divisor in the completed candidate strip |
| Auxiliary VF probes | Centrally assigned progression blocks evaluated by product workers | No general divisor-strip certificate from failure alone |

A Fermat worker seeks directly to its assigned midpoint. It does not have
to scan from $\lceil\sqrt N\rceil$ first. A product worker owns its arithmetic
state and normalized-polynomial caches, reusing those caches across jobs.
The coordinator assigns distinct new tiles within each search queue and
requeues an unfinished tile from the first uncertified position.

These queues cooperate through the ledger, but in-flight product and Fermat
jobs can still overlap in the factor regions they ultimately cover. Work
already running cannot always benefit from a certificate another worker has
just produced.

Before creating threads, the default configuration performs up to 64 cheap
Fermat square tests. Very close factor pairs can therefore finish without
thread startup. The remaining scheduler estimates the time already assigned
to each search from completed job durations plus average-cost estimates for
in-flight jobs. It prefers Fermat while its assigned time is at most
`fermat_ratio` times the product time. The default is a soft 2:1 target;
it does not fix two cores to Fermat or prescribe a 2:1 job-count ratio.
Available jobs, discoveries, and budgets can change the realized allocation.

## Exact certificates despite out-of-order completion

For odd $N$, an actual divisor $d\le\sqrt N$ has integer complementary
factor $N/d$ and integer midpoint

$$
R(d)=\frac{d+N/d}{2}.
$$

This map decreases on the smaller-factor domain. Thus a candidate band
$[l,u]$ maps to the integer midpoint interval

$$
\left[
\left\lceil\frac{u^2+N}{2u}\right\rceil,
\left\lfloor\frac{l^2+N}{2l}\right\rfloor
\right].
$$

`midpoint_interval_for_factors` calculates these bounds using integer
arithmetic. The Fermat cursor skips midpoints rejected by the exact
congruence filters and tests the remaining $R^2-N$ values for squareness.
Rejected residues count as excluded positions even though no square test is
needed there.

Suppose a job has excluded every integer midpoint from $A$ through $B$.
Because an actual divisor's midpoint is an integer, this also excludes
possible divisors satisfying

$$
A-1<R(d)<B+1.
$$

Writing $f(x)=x-\sqrt{x^2-N}$, inversion on the decreasing branch gives
$d>f(B+1)$ and, when applicable, $d<f(A-1)$. Rounding these strict bounds
to integer divisor coordinates yields the certificate

$$
L_B=B+2-\left\lceil\sqrt{(B+1)^2-N}\right\rceil,
$$

$$
U_A=
\begin{cases}
\lfloor\sqrt N\rfloor,
    & A=\lceil\sqrt N\rceil,\\
A-2-\left\lfloor\sqrt{(A-1)^2-N}\right\rfloor,
    & A>\lceil\sqrt N\rceil.
\end{cases}
$$

Here $A,B$ are clipped to the mathematical midpoint domain, and the final
certificate $[L_B,U_A]$ is clipped to the requested factor band. These bounds
are implemented by `excluded_integer_midpoint_interval`, called through
`RegionFermatCursor::excluded_factors()`.

The proof uses the integer midpoint of an **actual factor pair**. It does
not assert that every integer candidate inside the certificate has its real
midpoint in $[A,B]$. The strict inequalities retain possible hits at $A-1$
and $B+1$. This integer-aware rounding lets adjacent completed midpoint
windows combine without artificial holes between them. A more conservative
real-interval helper remains available for comparison, but the coordinator
uses the integer-aware helper above.

The cursor's `covered_through` certifies only its own interval
`[start, covered_through]`. It is not a prefix beginning at
$\lceil\sqrt N\rceil$. Applying the old global Fermat upper-bound formula
directly to a distant worker's endpoint would incorrectly discard unsearched
regions. A successful midpoint is never included in a failed certificate;
an earlier completed prefix can still be committed.

`CoverageLedger` stores the union of certified closed factor intervals and
merges overlaps and adjacent intervals. It preserves every uncovered gap,
regardless of job completion order. Its first and last uncovered positions
constrain future jobs. A completed middle interval alone cannot advance a
frontier across an unsearched hole. Only actual certificates are inserted;
an unfinished polynomial strip does not become excluded merely because a
worker attempted it.

## Shared budgets and cancellation

The represented-work allowance counts admissible Fermat square tests,
scalar divisor tests, and positions represented by successfully evaluated
product blocks, including auxiliary VF blocks. Recovery revisits already
represented positions without charging them twice. Congruence-rejected
midpoint gaps do not consume square-test units. These units describe search
coverage, not a fixed number of CPU instructions; costs differ between
methods and input sizes.

Before starting a job, the coordinator reserves part of the remaining work
allowance. At completion it charges actual work and refunds the unused
reservation. Under the coordinator mutex it maintains

$$
\text{spent}+\text{reserved}\le\text{global allowance}.
$$

`--max-steps` is therefore a shared cap, not a separate allowance for each
thread. Small remaining quotas can use scalar product prefixes so that a
job can make progress without requiring a whole polynomial block.

Each job observes the common deadline and its local time quantum. Fermat
workers check cancellation between chunks of at most 16384 admissible tests;
the cursor also checks its time budget internally. Product workers check
between bounded operations and scalar tails. A factor is verified by exact
division before the coordinator requests cancellation of the other workers.
All threads are joined before a result is returned, and no work reservations
remain outstanding.

Cancellation is cooperative. An individual GMP operation already in progress
is allowed to finish, so the deadline and job quantum are not hard real-time
limits. Only the completed prefix of a cancelled or timed-out job is
certified; its remaining range can be retried.

Mutable GMP arithmetic objects and polynomial caches belong to individual
workers. Shared coordination state is mutex-protected, and the common input
integer is read-only. Worker caches improve reuse but also mean that memory
use and some setup work grow with worker count.

## Interpreting results

| Single-search status | Meaning |
|---|---|
| `factor` | A verified proper divisor was found; `cofactor` is its exact quotient |
| `prime` | Primality was proved by a precheck or complete coverage of the full smaller-divisor domain |
| `probable_prime` | The primality precheck reported an unproved probable prime |
| `region_exhausted` | The requested restricted band was completely certified; this alone makes no global primality claim |
| `budget` | Search stopped with some requested coverage still missing |

`details.coverage_intervals` contains the merged certificates.
`details.lo` and `details.hi` are the first and last uncovered candidates,
not a claim that every position between them is unsearched. Read the interval
list to distinguish internal covered regions from gaps. With `--trace`,
each job records its requested coordinate range, committed factor certificate,
consumed work, worker identity, and the resulting outer uncovered bounds.
For Fermat jobs the requested coordinates are midpoints; the certificate
coordinates are smaller-divisor candidates.

`elapsed_seconds` measures wall time. `process_cpu_seconds` sums process CPU
use and can exceed wall time when threads execute concurrently. Counters
include actual work, jobs by search, cancellations, peak active workers, and
reservation high-water marks. Precheck-only results do not contain a searched
band or coverage certificates; consult `details.phase` before interpreting
the default range fields.

Full-factorization output separates `prime_factors`,
`probable_prime_factors`, and `unresolved` cofactors. Their product, with
multiplicity, is the input. Its statuses are `complete`, `probable_complete`,
and `partial`; individual search results appear in `runs`.

## Separate machines with explicit bands

Choose inclusive, adjacent candidate bands with the next lower boundary
equal to the preceding upper boundary plus one. For example,
$\lfloor\sqrt{1009000007063}\rfloor=1004489$, so three machines can cover
the domain with these independent commands:

```bash
# Machine A
./build/vf_factor_gmp_parallel 1009000007063 \
    --factor-lo 7 --factor-hi 1000 \
    --threads 4 --timeout 60 --no-vf > region-a.json

# Machine B
./build/vf_factor_gmp_parallel 1009000007063 \
    --factor-lo 1001 --factor-hi 10000 \
    --threads 4 --timeout 60 --no-vf > region-b.json

# Machine C
./build/vf_factor_gmp_parallel 1009000007063 \
    --factor-lo 10001 --factor-hi 1004489 \
    --threads 4 --timeout 60 --no-vf > region-c.json
```

These bands partition candidate/witness values. They do **not** guarantee
that a returned factor lies in its assigned band or that only one machine
can find a particular factor. A product GCD can recover a small factor from
a multiple within a higher band: in this example, the candidate
$11099=11\cdot1009$ in machine C's band can reveal factor 1009. The shared
prechecks also run globally, and optional VF probes use their own global
progressions. `--no-vf` removes those auxiliary progressions, but does not
disable prechecks or change the product-GCD behavior.

This does not weaken a failed-strip certificate: after a strip is fully
certified, there is no divisor in that checked interval. Any returned proper
factor is useful globally, regardless of which band produced its witness.

Each machine has its own deadline and work cap. An external operator or
controller must collect results and stop the remaining processes after a
factor is found. `region_exhausted` certifies its whole requested band;
`budget` certifies only the reported completed intervals. A global exhaustion
argument requires the union of certificates to cover the full domain with
the common prechecks accounted for. This executable does not perform that
cross-machine merge or resume from remote certificate files.

Physical distance does not change the factorization mathematics. The useful
resources are additional compute capacity and coordinated candidate ranges;
network latency matters to whatever external controller supplies the jobs.

## C++ API

Include `gmp_parallel.hpp` and link as shown in the build section:

```cpp
#include "gmp_parallel.hpp"

vf::ParallelOptions options;
options.threads = 4;
options.degree = 128;
options.algorithm.timeout = 10.0;
options.algorithm.max_steps = 1000000;
options.algorithm.trace = true;

const vf::Z n("978508015703");
vf::ParallelResult result = vf::factor_parallel(n, options);
```

`factor_parallel_full(n, options)` returns `ParallelFullResult` under one
overall deadline and work allowance. For an existing outer deadline,
`parallel_find_factor(n, options, global_time, max_work)` accepts a `Budget&`
and a separate optional shared work allowance. Its deadline must be finite;
it manages its own synchronized work reservations rather than using the
outer budget's work counter.

`ParallelOptions` exposes `threads`, `degree`, `fermat_span`, `job_seconds`,
optional `factor_lo`/`factor_hi`, and the common `algorithm` options. The
parallel scheduler honors the common toggles for VF, Fermat, polynomial
sharing, native arithmetic, tracing, and probable-prime checks, plus the
initial Fermat tests and time ratio. Its fixed degree and exact ledger replace
the sequential adaptive-degree and feedback controls.

For custom dispatchers, `RegionFermatCursor(n, start, end, native64)` supports
`advance(max_candidates, seconds)`, `certified_interval()`, and
`excluded_factors()`. Its default backend retains native64 arithmetic where
applicable and uses GMP for larger inputs. `CoverageLedger` provides interval
insertion, completeness, uncovered bounds, gaps, and exact total certified
length; callers must synchronize access when sharing a ledger.

## Validation and performance measurement

Build and run the C++ correctness suites, then the Python integration tests:

```bash
make -f Makefile.parallel selftest
python3 -m pip install -r requirements.txt
python3 -m pytest -q tests/test_parallel_cpp.py
```

The suites cover rounded midpoint boundaries, arbitrary-size seeks,
native/GMP agreement, out-of-order ledger gaps, regional product recovery,
partial-job certificates, shared budget limits, cancellation, and CLI/result
semantics. Mathematical certificates are checked against direct small-input
oracles as well as larger arithmetic cases.

For the frozen scaling experiment, also build the sequential comparison:

```bash
make -f Makefile.gmp
python3 parallel_benchmark.py \
    --timeout 1 --repeats 3 --outdir results/parallel
```

That harness requires at least four exposed logical CPUs and compares the
same regional policy with one, two, and four workers. The older sequential
GMP algorithm is a separate algorithm control. It records wall time, process
CPU, represented work, environment metadata, and unsolved cases. Compare
worker counts under the same policy to assess parallel scaling; a difference
from the sequential executable also includes changes in scheduling and
regional product strategy.

There is no guaranteed speedup. Close-factor inputs may finish in the serial
seed pass. For other inputs, thread startup, repeated cache construction,
memory pressure, overlapping in-flight search, and available CPU capacity
can outweigh the gain from concurrent jobs. Completion order can change the
amount of work performed before a hit, so report success rates, wall time,
CPU time, and work together. More workers provide additional parallel search
capacity; they do not change the underlying asymptotic factorization method.
