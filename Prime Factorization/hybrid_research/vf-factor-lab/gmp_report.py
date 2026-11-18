"""Generate the GMP extension report from frozen measurements; run no searches."""
from __future__ import annotations
import html
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "results/gmp"
FRESH = json.loads((DATA/"fresh/benchmark-summary.json").read_text())
LARGE = json.loads((DATA/"large/benchmark-summary.json").read_text())
REPLAY = json.loads((DATA/"original-replay/benchmark-summary.json").read_text())
PROTOCOL = json.loads((DATA/"protocol.json").read_text())
ENV = json.loads((DATA/"environment.json").read_text())
AUDIT = json.loads((DATA/"audit.json").read_text())
EXAMPLE = json.loads((DATA/"example.json").read_text())
LABELS = {"gmp_cpp":"C++ / GMP (default)", "python_fused":"Python / FLINT reference",
          "gmp_cpp_all_mpz":"C++ with GMP Fermat at every size"}


def table(headers, rows):
    return ('<div class="scroll"><table><thead><tr>'
            + ''.join('<th>'+str(x)+'</th>' for x in headers)
            + '</tr></thead><tbody>'
            + ''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)
            + '</tbody></table></div>')


def select(rows, method, family=None, bits=None):
    return [r for r in rows if r['algorithm']==method and
            (family is None or r['family']==family) and
            (bits is None or r['target_bits']==bits)]


def paired(rows, family=None):
    cpp = {r['case_id']:r for r in select(rows,'gmp_cpp',family) if r['family']!='prime'}
    py = {r['case_id']:r for r in select(rows,'python_fused',family) if r['family']!='prime'}
    ratios = [py[k]['median_seconds']/v['median_seconds'] for k,v in cpp.items()
              if v['solved'] and py[k]['solved']]
    return dict(pairs=len(ratios),cpp_faster=sum(r>1 for r in ratios),
                median_ratio=statistics.median(ratios),geometric_ratio=statistics.geometric_mean(ratios))


def timing_cell(rows):
    done = [r for r in rows if r['solved']]
    value = f"{statistics.median(r['median_seconds'] for r in done)*1000:,.4f} ms" if done else '—'
    return value+f'<small>{len(done)}/{len(rows)} completed</small>'


comparisons = [("Fresh composites",paired(FRESH)),
               ("Fresh balanced",paired(FRESH,'balanced')),
               ("Fresh close factors",paired(FRESH,'close')),
               ("Fresh unbalanced",paired(FRESH,'unbalanced')),
               ("Large completed cases",paired(LARGE)),
               ("Large close factors",paired(LARGE,'close')),
               ("Large small factors",paired(LARGE,'small_factor'))]
comparison_table = table(["Subset","Matched successes","C++ faster","Median Python/C++ ratio","Geometric ratio"],
    [[name,r['pairs'],f"{r['cpp_faster']}/{r['pairs']}",f"{r['median_ratio']:.3f}×",f"{r['geometric_ratio']:.3f}×"]
     for name,r in comparisons])
completion_rows = []
for dataset,rows,methods in [("Fresh",FRESH,list(LABELS)),("Original replay",REPLAY,['gmp_cpp']),
                             ("80–256 bits",LARGE,['gmp_cpp','python_fused'])]:
    for method in methods:
        chosen = select(rows,method)
        composites = [r for r in chosen if r['family']!='prime']
        primes = [r for r in chosen if r['family']=='prime']
        completion_rows.append([dataset,LABELS[method],
            f"{sum(r['solved'] for r in composites)}/{len(composites)}",
            f"{sum(r['solved'] for r in primes)}/{len(primes)}" if primes else '—'])
completions = table(["Dataset","Implementation","Composites completed","Prime controls"],completion_rows)
balanced = table(["Bits"]+list(LABELS.values()),
    [[bits]+[timing_cell(select(FRESH,method,'balanced',bits)) for method in LABELS]
     for bits in (24,32,40,48,56,64)])
large_times = table(["Bits","Family","C++ / GMP","Python / FLINT"],
    [[bits,family]+[timing_cell(select(LARGE,method,family,bits)) for method in ('gmp_cpp','python_fused')]
     for bits in (80,96,128,192,256) for family in ('close','small_factor','balanced')])
files = table(["File","Role"],[
    ['cpp/gmp_factor.cpp','Standalone JSON CLI, batch mode, and metadata'],
    ['cpp/gmp_factor.hpp','Unified coordinator, exact averages, recursive shared budget, reusable C++ API'],
    ['cpp/gmp_common.hpp','GMP integer type, counters, and deadline/work accounting'],
    ['cpp/gmp_arithmetic.hpp','Exact roots, prime-status handling, reciprocal bounds, persistent Fermat cursor'],
    ['cpp/fermat_native_bridge.hpp','C++ bridge to the unchanged validated 64-bit kernel'],
    ['cpp/gmp_products.hpp','Normalized block products, GMP polynomial engine, and aggregate GCD recovery'],
    ['Makefile.gmp','Build the executable and independent C++ self-tests without Python'],
    ['tests/test_gmp_cpp.py','Cross-language integration and resource-limit tests'],
    ['gmp_benchmark.py','Frozen comparisons, corpora, validation, and raw measurements'],
])
test_text=(ROOT/'results/gmp-full-test-results.txt').read_text()
test_summary=next(line for line in reversed(test_text.splitlines()) if 'passed' in line)
build_text=(ROOT/'results/gmp-build-results.txt').read_text()
cpp_checks=[json.loads(line) for line in build_text.splitlines() if line.startswith('{')]
assert all(c['status']=='passed' for c in cpp_checks)
assert sum(c['checks'] for c in cpp_checks)==704138
raw_calls=sum(len((DATA/path/'benchmark-runs.jsonl').read_text().splitlines())
              for path in ('fresh','original-replay','large'))
assert raw_calls==1117
example_result=EXAMPLE['result']
assert example_result['status']=='complete' and not example_result['unresolved']
analysis=dict(paired={name:data for name,data in comparisons},cpp_assertions=704138,
              pytest_summary=test_summary,timed_calls=raw_calls,
              conclusion='Standalone GMP port works. Standard-corpus median paired improvement is 2.621x; large small-factor searches are slower than the FLINT reference on 7/10 inputs. Balanced 80–256-bit cases remain unresolved.')
(DATA/'report-analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')

document=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>C++ / GMP unified factorizer — implementation and validation</title>
<style>
:root{{color-scheme:light;--ink:#18364c;--muted:#596c7c;--line:#dce7ef}}*{{box-sizing:border-box}}
body{{margin:0;background:#f3f7fa;color:var(--ink);font:16px/1.65 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1280px;margin:auto;padding:32px 26px 64px}}header{{background:#123e50;color:white;padding:32px;border-radius:16px}}
h1{{font-size:35px;line-height:1.2;margin:12px 0 20px}}header p{{color:#e3eff5;font-size:19px}}.tag{{letter-spacing:.1em;text-transform:uppercase;font-size:12px;color:#a8dcd7}}
h2{{font-size:25px;margin:24px 0 12px}}h3{{font-size:19px}}section{{background:white;padding:10px 24px 24px;border:1px solid var(--line);border-radius:12px;margin:24px 0}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:20px 0}}.card{{background:white;padding:18px;border:1px solid var(--line);border-radius:10px}}.card b{{display:block;font-size:29px;color:#08796d}}.card span{{font-size:14px;color:var(--muted)}}
.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;font-size:14px;margin:16px 0}}td,th{{text-align:left;vertical-align:top;padding:10px;border-bottom:1px solid var(--line)}}th{{background:#eaf2f7}}small{{display:block;color:var(--muted);font-size:12px}}
pre{{background:#102a3b;color:#e9f2f7;border-radius:8px;padding:18px;overflow:auto;font-size:13px}}code{{font-family:ui-monospace,SFMono-Regular,monospace}}.note{{border-left:4px solid #bd892b;background:#fff7e8;padding:14px 18px}}.good{{border-left:4px solid #168071;background:#edf8f3;padding:14px 18px}}.math{{background:#edf4f8;padding:15px 19px;overflow:auto;font:19px/1.8 Georgia,serif}}a{{color:#075f98}}summary{{cursor:pointer;font-weight:650}}
@media(max-width:720px){{main{{padding:16px}}header{{padding:22px}}h1{{font-size:28px}}section{{padding:8px 14px 18px}}.cards{{grid-template-columns:repeat(2,1fr)}}}}
</style></head><body><main>
<header><div class="tag">Number Theory and Factorization · 17 November 2026 · version 0.3.0</div>
<h1>The unified factorizer now runs in C++ with GMP</h1>
<p>A standalone arbitrary-precision implementation of the shared Fermat, VF, and divisor-product search, with its own fast polynomial engine and one recursive resource budget.</p></header>
<div class="cards"><div class="card"><b>75 / 75</b><span>standard fresh composites completed</span></div>
<div class="card"><b>2.621×</b><span>median paired speedup over Python on those 75 inputs</span></div>
<div class="card"><b>20 / 25</b><span>larger composites completed; five hard balanced cases unresolved</span></div>
<div class="card"><b>5,859</b><span>tests passed, plus 704,138 C++ assertions</span></div></div>

<section><h2>Build and run</h2>
<p>The executable needs GMP development headers/libraries and GCC or Clang with C++17 and unsigned 128-bit support. It does not need Python or FLINT at runtime. These commands are for Linux or WSL, inside the extracted project directory.</p>
<pre>make -f Makefile.gmp
./build/vf_factor_gmp 978508015703 --full --timeout 1
./build/vf_factor_gmp --metadata</pre>
<p>To compile directly:</p>
<pre>mkdir -p build
g++ -O3 -std=c++17 cpp/gmp_factor.cpp vf_factor/fusion_kernel.cpp \\
    -lgmpxx -lgmp -o build/vf_factor_gmp</pre>
<h3>A verified 128-bit example</h3>
<pre>{html.escape(EXAMPLE['command'])}</pre>
<div class="math">{example_result['input']}<br>= {example_result['prime_factors'][0]} × {example_result['prime_factors'][1]}</div>
<p>This execution uses the GMP Fermat backend and finishes with two proven 64-bit prime factors and no unresolved cofactor. The factors differ by 42; it is a favorable close-factor example, not evidence about general balanced 128-bit difficulty.</p>
<details><summary>Complete example output</summary><pre>{html.escape(json.dumps(example_result,indent=2))}</pre></details></section>

<section><h2>The complete port</h2>{files}
<p>Input integers, factors, midpoint state, VF averages and cursors, and reciprocal bounds remain GMP integers. Small loop indices are converted only after a configured strip-width bound. Below 2⁶⁴ the default reuses the existing fast native Fermat kernel. Above that threshold it uses GMP differences, perfect-square screening, and exact roots. The <code>--gmp-fermat</code> flag forces GMP at every size for comparison.</p>
<p>The port preserves the exact period-160 sieve, published VF level selection, full divisor strips, normalized polynomial sharing, and exceptional GCD recovery. A strip is excluded only when its blocks, tails, and recovery finish. Failed VF probes do not exclude a divisor interval.</p>
<p>For a remaining smaller divisor d, the shared bounds remain:</p>
<div class="math">ℓ ≤ d ≤ B − ⌈√(B² − N)⌉<br>Fermat midpoint cap = ⌊(ℓ² + N)/(2ℓ)⌋</div>
<p>B is the next admissible untested midpoint, while every possible divisor below ℓ has been excluded. Both bounds and their integer rounding are exact.</p></section>

<section><h2>GMP polynomial acceleration</h2>
<p>The standalone backend builds block polynomials and evaluation-point polynomials with balanced product trees. Larger multiplications pack nonnegative modular coefficients into separated integer slots, call GMP integer multiplication, and unpack the result. A slot width of at least 2b + ⌈log₂t⌉ bits suffices when coefficients have at most b bits and each convolution coefficient has at most t summands. The implementation rounds slots up to bytes and checks allocation bounds.</p>
<p>Multipoint evaluation uses monic product/remainder trees. Larger remainders use reversed-polynomial division and Newton series inversion. The reversed monic divisor has constant coefficient one, so this arithmetic is valid over composite moduli. Small operations use schoolbook arithmetic to avoid packing overhead.</p>
<div class="math">P<sub>m</sub>(X) = ∏<sub>j=0</sub><sup>m−1</sup>(X+j)<br>
∏(α+sj) ≡ s<sup>m</sup>P<sub>m</sub>(αs<sup>−1</sup>) (mod N)</div>
<p>Steps 30 and 5 are units after the common small-factor checks. Their requests use one normalized polynomial per degree and a combined GCD stage. This retains fast block evaluation instead of individually multiplying every represented candidate.</p>
<p class="note">The GMP implementation is independent of FLINT, but independence does not guarantee better performance. FLINT remains faster on most large small-factor cases in this experiment.</p></section>

<section><h2>Matched benchmark results</h2>
<p>All comparisons use a 0.5-second search budget, three repetitions on success, one pinned logical CPU, and one FLINT thread. A failed repetition ends further attempts and makes the aggregate case unsolved. The fresh seed is 926491; the 80–256-bit seed is 926493. Generated semiprimes are separate from previous tuning and validation seeds; published examples and fixed regression controls remain reused.</p>
<p>Primary timings are internally measured search time. Process/import startup is excluded, per-input arithmetic setup is included, and external request wall time is retained. The Python reference is required to have its native 64-bit kernel available. Adaptive time allocation can change search trajectories between implementations and repetitions.</p>
{completions}
<p class="good">The default C++ implementation and Python reference both complete every fresh standard composite and prime control. C++ also completes every original-corpus replay input. On the 75 fresh composites, C++ is faster on 60, with a median paired Python/C++ ratio of 2.621×.</p>
{comparison_table}
<p>Ratios above one favor C++. Ratios are formed separately for each matched successful input before aggregation. The large small-factor median ratio is 0.853×: C++ is slower on seven of ten such inputs. This result is not hidden inside the overall improvement.</p></section>

<section><h2>Balanced standard inputs</h2>
<p>Each cell reports a median over completed case medians plus the completion count. Four inputs are tested per size. Partial-group timings do not represent the missing cases.</p>
{balanced}
<p>The default C++ median is 313.303 ms on the four balanced 64-bit inputs, versus 238.169 ms for Python/FLINT. Thus the overall 2.621× result is not a uniform gain at every size. The forced-GMP cursor completes 70/75 standard composites, versus 75/75 with the retained native64 path, supporting the default backend choice.</p></section>

<section><h2>Larger integers: capabilities and limits</h2>
{large_times}
<p>Both implementations complete ten close-factor and ten small-factor inputs spanning 80, 96, 128, 192, and 256 bits. Both exhaust the budget on all five nonclose balanced inputs. Small-factor cases contain a 16- or 20-bit smaller prime. Bit length alone does not describe their factorization difficulty.</p>
<p>The C++ median paired advantage on the ten close-factor inputs is 4.280×. On a 256-bit input with a 20-bit factor, however, C++ takes 93.897 ms versus 46.338 ms for Python/FLINT. Arbitrary-precision support is established; general balanced 256-bit factorization is not.</p></section>

<section><h2>Correctness, budgets, and primality status</h2>
<p><strong>{html.escape(test_summary)}</strong></p>
<p>The independent C++ suites pass 476,174 arithmetic assertions and 227,964 polynomial assertions. Arithmetic tests include native/GMP cursor parity, large resume/coverage oracles, and selected inputs of about 4,097 bits. Polynomial tests include composite moduli through 512 bits, degrees through 512, duplicate points, signed progressions, coefficient carry boundaries, monic remainders, and GCD-N recovery. These are correctness checks, not benchmarks of difficult factorization.</p>
<p>The new 51-test integration module includes loops over all inputs 2 through 1001, larger factors and squares, status handling, recursive multiplicities, malformed input recovery, and exact fixed-work trace/counter agreement with Python. The polynomial suite also passes AddressSanitizer and UndefinedBehaviorSanitizer; leak detection was disabled because this environment does not support LeakSanitizer under its process tracing.</p>
<p>One budget covers all components and all <code>--full</code> recursion. Represented positions are charged even when completed evaluation is followed by a deadline exception. Recovery revisits do not charge the same positions again. Failed recovery cannot advance a divisor boundary. Per-kind diagnostic attribution can be incomplete when an operation is interrupted; total represented work remains charged.</p>
<p>Primality below 2⁶⁴ uses deterministic bases. Larger GMP probable-prime results are explicitly labeled <code>probable_prime</code>; <code>--full</code> separates <code>prime_factors</code>, <code>probable_prime_factors</code>, and <code>unresolved</code>. The product of those lists, with multiplicities, remains the original input. A probable-prime leaf is not presented as a proof.</p>
<p>Independent audit found no invalid factors, false prime claims, algorithm errors, work-accounting mismatches, or summary errors across {raw_calls:,} timed calls. Source, binary, and corpus hashes match. One forced-GMP control input, <code>balanced-56-0</code>, completed one repeat and timed out on the next; it is correctly counted as unsolved. Default C++ and Python completion statuses are consistent across all repetitions.</p>
<details><summary>C++ self-test records</summary><pre>{html.escape(json.dumps(cpp_checks,indent=2))}</pre></details>
<details><summary>Independent benchmark audit</summary><pre>{html.escape(json.dumps(AUDIT,indent=2))}</pre></details></section>

<section><h2>Reproduce and reuse</h2>
<pre>make -f Makefile.gmp all selftest

# Python/reference tests and matched benchmark:
python3 -m pip install -r requirements.txt
python3 build_fusion.py
python3 build_native.py
python3 build_wheel20.py
python3 -m pytest -q
python3 gmp_benchmark.py --all-mpz-control --outdir results/gmp-reproduction

# Regenerate this report from the bundled original measurements:
python3 gmp_report.py</pre>
<p>The reusable C++ API is <code>vf::factor</code>, <code>vf::factor_full</code>, or <code>vf::find_factor</code> with a caller-owned shared budget. See <code>cpp/README.md</code> for usage. JSON encodes arithmetic integers as decimal strings to preserve exact values across languages.</p>
<p>No settings were changed in response to the benchmark outcomes. The fixed degree cap remains 512 by default. This port changes implementation and numeric range; it establishes no new factorization complexity exponent.</p>
<details><summary>Environment</summary><pre>{html.escape(json.dumps(ENV,indent=2))}</pre></details>
<details><summary>Frozen protocol</summary><pre>{html.escape(json.dumps(PROTOCOL,indent=2))}</pre></details></section>

<section><h2>Primary documentation and source</h2><ul>
<li><a href="https://gmplib.org/manual/C_002b_002b-Interface-General">GMP C++ interface and linking.</a></li>
<li><a href="https://gmplib.org/manual/Number-Theoretic-Functions">GMP GCD, modular inverses, and distinct definite/probable prime statuses.</a></li>
<li><a href="https://gmplib.org/manual/Integer-Roots">GMP exact integer roots and perfect-square tests.</a></li>
<li><a href="https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Complete_Fermat_Sieve_Verified.pdf">Complete Fermat sieve reference.</a></li>
</ul><p>The original project sources supply the iterated-average progressions. The preceding unified implementation and report document the reciprocal bounds and normalization identity. All measurements in this report come from the bundled raw calls.</p></section>
</main></body></html>'''
(ROOT/'results/gmp_report.html').write_text(document)
print(ROOT/'results/gmp_report.html')
