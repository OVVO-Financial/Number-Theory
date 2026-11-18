"""Generate the regional parallel report from saved evidence; run no searches."""
from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "results/parallel"


def read(name):
    return json.loads((DATA / name).read_text())


def esc(value):
    return html.escape(str(value))


def table(headers, rows):
    return ('<div class="scroll"><table><thead><tr>'
            + ''.join('<th scope="col">' + esc(x) + '</th>' for x in headers)
            + '</tr></thead><tbody>'
            + ''.join('<tr>' + ''.join('<td>' + str(x) + '</td>' for x in row)
                      + '</tr>' for row in rows)
            + '</tbody></table></div>')


ANALYSIS = read('analysis.json')
AUDIT = read('audit.json')
PROTOCOL = read('protocol.json')
ENV = read('environment.json')
SUMMARY = read('benchmark-summary.json')
CORPUS = read('benchmark-corpus.json')
DEMO = read('region-demo.json')
SANITIZERS = json.loads((ROOT / 'results/parallel-sanitizer-results.json').read_text())
assert AUDIT['audit_status'] == 'passed'
METHODS = ['regional_t1', 'regional_t2', 'regional_t4', 'sequential_gmp']
LABELS = ['Regional: 1 worker', 'Regional: 2 workers', 'Regional: 4 workers', 'Earlier sequential GMP']
FAMILIES = {
    'balanced': 'Balanced, 48–64 bits',
    'small_factor': 'Small factors, 64/128-bit inputs',
    'close_control': 'Close controls, 64–256 bits',
    'hard_balanced': 'Hard balanced, 80 bits',
}
lookup = {(r['case_id'], r['algorithm']): r for r in SUMMARY}

completion_rows = []
for family, label in FAMILIES.items():
    counts = AUDIT['cohort_completion_counts'][family]
    completion_rows.append([label] + [f"{counts[m]['completed']} / {counts[m]['total']}" for m in METHODS])
completion_rows.append(['<strong>All inputs</strong>'] + [
    f"<strong>{ANALYSIS['outcomes'][m]['completed']} / 26</strong>" for m in METHODS])
completions = table(['Input family'] + LABELS, completion_rows)

scaling_rows = []
for row in ANALYSIS['paired']:
    if row['reference'] != 'regional_t1' or row['family'] is not None:
        continue
    cpu = next(x for x in AUDIT['joint_success_comparisons']
               if x['reference'] == 'regional_t1' and x['comparator'] == row['comparator'] and x['family'] is None)
    scaling_rows.append([
        '2 vs 1' if row['comparator'] == 'regional_t2' else '4 vs 1',
        'Exclude close controls' if row['omit_close_controls'] else 'All joint successes',
        row['matched_successes'], f"{row['comparator_faster']} / {row['matched_successes']}",
        f"{row['median_ratio']:.3f}×", f"{row['geometric_ratio']:.3f}×",
        '—' if row['omit_close_controls'] else f"{cpu['cpu_cost_comparator_over_reference']['median']:.3f}×",
    ])
scaling = table(['Workers', 'Cohort', 'Matched cases', 'Faster cases',
                 'Median T₁/Tₖ', 'Geometric T₁/Tₖ', 'Median CPUₖ/CPU₁'], scaling_rows)

def timing_cell(row):
    if row['solved']:
        return f"{row['median_seconds'] * 1000:,.3f} ms"
    return f'<span class="budget">Budget</span><small>{row["median_seconds"]:.3f} s</small>'

case_rows = []
for case in CORPUS:
    rows = [lookup[case['case_id'], method] for method in METHODS]
    ratio = rows[0]['median_seconds'] / rows[2]['median_seconds'] if rows[0]['solved'] and rows[2]['solved'] else None
    case_rows.append(
        f'<tr data-family="{esc(case["family"])}"><td><code>{esc(case["case_id"])}</code>'
        f'<small>N = {case["n"]}<br>Known factors: {case["p"]} × {case["q"]}</small></td>'
        + ''.join('<td>' + timing_cell(r) + f'<small>CPU: {r["median_cpu_seconds"]*1000:,.3f} ms<br>'
                  f'Work: {r["median_work_units"]:,}</small></td>' for r in rows)
        + '<td>' + (f'{ratio:.3f}×' if ratio is not None else '—') + '</td></tr>')
per_case = ('<div class="scroll"><table id="cases"><thead><tr><th>Input</th>'
            + ''.join('<th>' + esc(x) + '</th>' for x in LABELS)
            + '<th>T₁ / T₄</th></tr></thead><tbody>' + ''.join(case_rows) + '</tbody></table></div>')

build_text = (ROOT / 'results/parallel-build-results.txt').read_text()
cpp_checks = [json.loads(line) for line in build_text.splitlines() if line.startswith('{')]
assert all(row['status'] == 'passed' for row in cpp_checks)
cpp_total = sum(row['checks'] for row in cpp_checks)
test_summary = next(line for line in reversed((ROOT / 'results/parallel-full-test-results.txt').read_text().splitlines()) if 'passed' in line)
audit_checks = sum(AUDIT['checks'].values())
extra = AUDIT['extra_completion_evidence'][0]

files = table(['File', 'Purpose'], [
    ['Makefile.parallel', 'Build the parallel CLI and independent native self-tests'],
    ['cpp/PARALLEL.md', 'Build, API, exact region mapping, budgets, and separate-machine instructions'],
    ['cpp/gmp_parallel.hpp / .cpp', 'Worker pool, scheduling, work reservations, recursive API, and JSON CLI'],
    ['cpp/gmp_regions.hpp / region_native.cpp', 'Exact midpoint seeking, interval certificates, and merged coverage ledger'],
    ['cpp/gmp_product_worker.hpp', 'Worker-local product state, polynomial caches, and bounded region work'],
    ['tests/test_parallel_cpp.py', 'Integration, out-of-order coverage, and concurrent resource-limit checks'],
    ['parallel_benchmark.py / parallel_report.py', 'Reproducible fixed protocol and evidence-only report generation'],
    ['results/parallel/', 'Every timed call, summaries, corpus, protocol, environment, and independent audit'],
])

document = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Parallel regional GMP factorization — implementation and measured scaling</title>
<style>
:root{color-scheme:light;--ink:#163147;--muted:#546c7b;--line:#dbe5eb}*{box-sizing:border-box}
body{margin:0;background:#f1f5f8;color:var(--ink);font:16px/1.65 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1320px;margin:auto;padding:28px 24px 60px}header{background:#163e52;color:white;padding:30px;border-radius:14px}
h1{font-size:35px;line-height:1.2;margin:12px 0 18px}header p{font-size:19px;color:#e2edf4}.tag{text-transform:uppercase;letter-spacing:.09em;font-size:12px;color:#afdcd4}
h2{font-size:25px;margin:24px 0 12px}h3{font-size:19px}section{background:white;border:1px solid var(--line);border-radius:12px;padding:8px 24px 24px;margin:22px 0}
.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:20px 0}.card{padding:18px;border:1px solid var(--line);background:white;border-radius:10px}.card b{display:block;font-size:29px;color:#08776b}.card span{font-size:14px;color:var(--muted)}
.scroll{overflow-x:auto}table{width:100%;border-collapse:collapse;font-size:14px;margin:16px 0}td,th{text-align:left;vertical-align:top;border-bottom:1px solid var(--line);padding:10px}th{background:#eaf1f5}small{display:block;color:var(--muted);font-size:12px}
pre{overflow:auto;white-space:pre;background:#102d3e;color:#edf5fa;border-radius:8px;padding:16px;font-size:13px}code{font-family:ui-monospace,SFMono-Regular,monospace}td code{font-size:12px}a{color:#08668e}
.note{border-left:4px solid #b08325;background:#fff8e9;padding:14px 18px}.good{border-left:4px solid #168071;background:#edf8f4;padding:14px 18px}.math{font:18px/1.8 Georgia,serif;background:#edf4f8;padding:14px;overflow:auto}.budget{color:#8d5020}summary{cursor:pointer;font-weight:650}select{font:inherit;padding:6px 10px;border:1px solid #afc1cd;border-radius:6px;background:white;color:var(--ink)}
@media(max-width:740px){main{padding:14px}header{padding:22px}h1{font-size:28px}section{padding:8px 14px 18px}.cards{grid-template-columns:repeat(2,1fr)}}
@media print{body{background:white}header{color:var(--ink);background:white;border:1px solid var(--line)}header p{color:var(--ink)}section{break-inside:auto}.cards{grid-template-columns:repeat(4,1fr)}select{display:none}}
</style></head><body><main>
'''
document += f'''
<header><div class="tag">Number Theory and Factorization · 17 November 2026 · version 0.4.0</div>
<h1>Searching several factor regions at once</h1>
<p>The C++ / GMP factorizer now has local parallel workers, exact region certificates, and explicit candidate bands for separate machines. The measured speed benefit varies substantially by input.</p></header>
<div class="cards"><div class="card"><b>20 / 26</b><span>completed with four regional workers</span></div>
<div class="card"><b>0.958×</b><span>median T₁/T₄ across 19 joint successes</span></div>
<div class="card"><b>3.84</b><span>measured CPU seconds per wall second, four workers</span></div>
<div class="card"><b>5,926</b><span>tests passed, plus {cpp_total:,} new C++ assertions</span></div></div>

<section><h2>What the implementation supports</h2>
<p>Each local worker receives a bounded Fermat midpoint window or a divisor-product strip. VF progression requests are allocated centrally and batched with product work. Workers own their mutable GMP arithmetic and polynomial caches. The coordinator merges only completed exclusion certificates and keeps all unproved gaps visible.</p>
<p>Local threads share a deadline, a represented-work allowance, a coverage ledger, and a cancellation flag. Explicit <code>--factor-lo</code> and <code>--factor-hi</code> options also make independent region jobs runnable on different machines. Network transport, cross-host cancellation, and a global cross-host budget are not implemented.</p>
<pre>make -f Makefile.parallel
./build/vf_factor_gmp_parallel 17315933965806733211 --threads 4 --timeout 2
./build/vf_factor_gmp_parallel 978508015703 --threads 4 --full --timeout 1
./build/vf_factor_gmp_parallel --metadata</pre>
<p>Build with a C++17 compiler supporting unsigned 128-bit arithmetic and GMP development headers/libraries. The executable has no Python or FLINT runtime dependency. Its Makefile links threads and GMP; <code>region_native.cpp</code> already includes the unchanged native kernel and must not be linked with a second copy of that kernel.</p>
{files}</section>

<section><h2>Measured results: concurrency works; overall speedup is unproven</h2>
<p>All methods received one second of wall time per input. There are 26 fixed-seed composites, three repetitions on success, and no additional repetitions after a failure. Factors were withheld from the executables and used only for independent validation. There were no mixed-status repetition groups.</p>
{completions}
<p class="note">The earlier sequential GMP implementation completed 23 inputs, compared with 20 for four regional workers. Its four additional completions all had 24-bit smaller factors. The regional implementation is a useful parallel experiment, but these results do not justify replacing the earlier implementation as the general default.</p>
<h3>One, two, and four workers under the same regional policy</h3>
{scaling}
<p>T₁/Tₖ above one means the multiworker search finishes sooner. Ratios below one indicate a regression. Ratios are computed per jointly completed input before aggregation; budget failures are reported separately and never treated as successful timings. Excluding cheap close-factor controls does not remove all short searches or all regressions. The filtered median of 1.298× must be read with its geometric mean of 0.737×.</p>
<p>Across the 19 joint one/four-worker successes, median CPU cost rose 2.091× and median represented work rose 1.078×. The scheduler changes region order and feedback timing, so these are observed time-to-factor results, not equal-work strong-scaling measurements.</p>
<h3>One additional completed input</h3>
<div class="math">{extra['n']} = {extra['p']} × {extra['q']}</div>
<p>Four workers succeeded on all three repetitions: 0.831688, 0.841226, and 0.832652 seconds. Median process CPU time was 3.284089 seconds. The one-worker and two-worker regional versions, and the older sequential version, each reached the one-second limit without a factor.</p>
<p>A 128-bit small-factor case improved from 242.687 ms with one regional worker to 117.540 ms with four; the older sequential implementation still took only 80.534 ms. Conversely, two easy 48-bit inputs regressed from 0.148 to 3.152 ms and from 0.202 to 2.369 ms. Saved counters show extra speculative work. Worker startup, polynomial setup, changed scheduling, and joining cancelled jobs are plausible contributors; separating their exact costs requires profiling.</p>
</section>

<section><h2>Why disjoint regions can be searched safely</h2>
<p>For an odd input N = pq with p ≤ q, the factor pair corresponds to a Fermat midpoint R and half-gap Y:</p>
<div class="math">R = (p + q)/2 = (p + N/p)/2<br>Y = (q − p)/2<br>R² − Y² = N</div>
<p>The function R(d) = (d + N/d)/2 decreases on 1 ≤ d ≤ √N. Small-divisor strips and near-root Fermat windows therefore approach the possible factors from opposite ends. A factor band [L,U] maps to the integer midpoint band:</p>
<div class="math">⌈(U² + N)/(2U)⌉ ≤ R ≤ ⌊(L² + N)/(2L)⌋</div>
<p>The worker initializes exact recurrence state directly at its assigned midpoint; it need not scan all preceding windows. The complete wheel-20 condition intersected with QR64 gives the existing period-160 sieve. Skipped midpoint classes cannot contain a square difference. Tested square differences and recovered factors are verified exactly.</p>
<p>A fully checked window [A,B] in which no factor was found certifies the corresponding smaller-factor interval. A factor-hit job certifies only its preceding failed portion. Because an actual odd-factor-pair midpoint is an integer, the certificate uses the open real interval A−1 &lt; R(d) &lt; B+1:</p>
<div class="math">L<sub>cert</sub> = B + 2 − ⌈√((B+1)² − N)⌉<br>
U<sub>cert</sub> = ⌊√N⌋ if A = ⌈√N⌉;<br>
U<sub>cert</sub> = A − 2 − ⌊√((A−1)² − N)⌋ otherwise.</div>
<p>Bounds are clamped to the assigned domain. The strict endpoints preserve possible hits at A−1 and B+1. Adjacent completed windows consequently meet without artificial rounding gaps. All calculations use exact integer roots.</p>
<p class="good">The coverage ledger stores an interval union, not just a single moving frontier. If a worker finishes an interior region while an earlier region is unfinished, that hole remains searchable. A new frontier advances only across a contiguous union of completed certificates. Failed VF probes never certify a divisor interval.</p>
<p>Within each branch, leases allocate separate coordinates, with unfinished tails returned for later work. Different branches can still perform overlapping speculative work before a certificate arrives. A product GCD can also recover the same factor from multiple candidate multiples. The design does not assert unique ownership of a factor.</p>
</section>

<section><h2>Scheduling, budgets, and thread safety</h2>
<p>The frozen policy uses degree 128 product blocks, midpoint windows of up to 655,360 raw positions, and a 25 ms job allocation. It starts with 64 retained Fermat tests before creating workers. The coordinator targets twice as much Fermat job time as product job time, accounting for completed time and estimates for in-flight jobs. This target is a scheduling heuristic, not a guaranteed CPU allocation ratio.</p>
<p>The older sequential control grows its degree from 16 to 512. Its different block sizes and scheduling policy confound comparisons with the regional engine. Only regional one/two/four-worker comparisons isolate the worker-count option, and even these can perform different search work.</p>
<p>Before a local job starts, it reserves part of the remaining work allowance. The invariant is <code>spent + reserved ≤ max_steps</code>. On completion, actual work is charged and unused reservation returned. Small remaining quotas use smaller batches or scalar work so a quota cannot trigger an endless zero-progress loop. Full recursive factorization shares the same global deadline and allowance across cofactors.</p>
<p>Cancellation is cooperative. Fermat checks between chunks of at most 16,384 retained square tests; product workers check between bounded operations. GMP calls are not interrupted internally. The coordinator joins every worker before returning. The largest measured overrun of the one-second deadline was 28.163 ms.</p>
<p>Workers keep mutable GMP objects separate and share N only for reading. The ledger, assignments, reservations, counters, and published result are mutex protected; cancellation is atomic. This follows the GMP manual’s rule permitting simultaneous reads while forbidding concurrent mutation of the same object. See <a href="https://gmplib.org/manual/Reentrancy">GMP: Reentrancy</a>.</p>
<p>Integer input, factors, and bounds are decimal strings in JSON. Above 64 bits, unproved GMP probable primes remain explicitly labeled <code>probable_prime</code>. A budget result proves no primality claim. A completely certified partial band reports <code>region_exhausted</code>.</p>
</section>

<section><h2>Separate machines and physical geographies</h2>
<p>The same binary can run these independent band jobs on different computers. The following commands were validated as three separate processes on one host; they are not evidence of a network deployment.</p>
<pre># Machine A
./build/vf_factor_gmp_parallel 1009000007063 --threads 2 --timeout 10 \\
    --factor-lo 7 --factor-hi 1000 --no-vf
# Machine B
./build/vf_factor_gmp_parallel 1009000007063 --threads 2 --timeout 10 \\
    --factor-lo 1001 --factor-hi 10000 --no-vf
# Machine C
./build/vf_factor_gmp_parallel 1009000007063 --threads 2 --timeout 10 \\
    --factor-lo 10001 --factor-hi 1004489 --no-vf</pre>
<p>A returns <code>region_exhausted</code>. B and C both find factor 1009, with cofactor 1,000,000,007. C may recover 1009 from a multiple of 1009 in its candidate band. Thus band bounds describe search candidates and coverage; they do not constrain the magnitude of a valid factor result. Global prechecks and enabled VF probes can also find factors outside the assigned band.</p>
<p>A production network coordinator would need job IDs, leases, failure/retry handling, central certificate merging, and stop broadcasts. Work units should be large enough to amortize network latency; the local 25 ms allocation has not been tuned for remote dispatch. For many independent N values, assigning whole inputs to different machines avoids most coordination. Physical distance supplies no mathematical advantage; extra compute and coordinated coverage do.</p>
</section>

<section><h2>Validation and reproducibility</h2>
<p><strong>{esc(test_summary)}</strong>, including 67 new integration tests. New C++ self-tests passed {cpp_total:,} assertions: exact mapping and ledger cases, native/GMP seek agreement, regional product recovery and quotas, and coordinator/API edge cases. The result audit independently passed {audit_checks:,} checks over 266 timed calls and 369 final certificate intervals. Source and executable hashes match the frozen run.</p>
<p>AddressSanitizer and UndefinedBehaviorSanitizer passed the regional self-test, with leak detection disabled because of the environment’s ptrace limitation. ThreadSanitizer produced no reported findings on the exercised four-worker searches. That TSan run preceded the final finite-deadline API guard; concurrency logic was unchanged by the guard, and final-source integration/API tests passed afterward. Sanitizers cover exercised paths and are not a proof of race freedom.</p>
<pre>make -f Makefile.parallel selftest
python3 -m pytest -q
python3 parallel_benchmark.py --outdir results/parallel-reproduction
python3 parallel_report.py</pre>
<p>The full Python suite requires the existing test dependencies and native reference kernels, as documented in README.md. The benchmark requires both C++ executables. The report generator only reads saved evidence; it does not rerun searches.</p>
<p>The environment exposed nine logical CPUs. Processes were restricted to one, two, or four of them. Aggregate process CPU/wall ratios, excluding close controls, were 0.998, 1.973, and 3.840. Cgroup quota files were unavailable, so these are observed CPU consumption figures, not claims about physical-core topology. The earlier executable’s CPU time uses 10 ms procfs ticks and is too coarse for its shortest calls.</p>
<p>Compilation used GCC {esc(ENV['cpp_metadata']['compiler'])} and GMP {esc(ENV['cpp_metadata']['gmp_version'])}. Corpus seed: {PROTOCOL['corpus_seed']}; randomized order seed: {PROTOCOL['order_seed']}. Executables were run sequentially with no simultaneous benchmark competitors. Internal per-input search time includes arithmetic setup and thread startup; process startup is excluded.</p>
<details><summary>Frozen protocol and hashes</summary><pre>{esc(json.dumps(PROTOCOL, indent=2))}</pre></details>
<details><summary>Sanitizer evidence</summary><pre>{esc(json.dumps(SANITIZERS, indent=2))}</pre></details>
</section>

<section><h2>All 26 inputs</h2>
<p>Primary values are median internal elapsed time. CPU and represented work are separate medians, not necessarily from the same repetition. Known factors shown here were never given to the algorithms. Four close controls can finish in the serial seed without creating workers.</p>
<label for="family">Show family: </label><select id="family"><option value="">All inputs</option>
{''.join('<option value="'+esc(k)+'">'+esc(v)+'</option>' for k,v in FAMILIES.items())}</select>
{per_case}
</section>
<section><h2>What this changes mathematically</h2>
<p>Parallelism makes several regions searchable at once and lets completed regions shrink the shared uncertainty. It does not identify a favorable region in advance, establish a new factorization exponent, or make arbitrary balanced large inputs easy. Both 80-bit hard balanced cases remained unresolved in every configuration under this one-second budget.</p>
<p>The current result is a tested basis for further scheduling experiments: amortizing very short searches, reducing speculative product work, and choosing block sizes for multiworker throughput. These would need fresh evaluation before being claimed as improvements.</p></section>
<footer><p>Evidence bundle: <code>vf_factor_lab.zip</code> · Report generated by <code>parallel_report.py</code> · Source and raw records retained alongside earlier experiments.</p></footer>
</main><script>
document.getElementById('family').addEventListener('change', function () {{
  document.querySelectorAll('#cases tbody tr').forEach(row => {{
    row.hidden = !!this.value && row.dataset.family !== this.value;
  }});
}});
</script></body></html>
'''

output = ROOT / 'results/parallel_report.html'
output.write_text(document)
(DATA / 'report-analysis.json').write_text(json.dumps({
    'pytest_summary': test_summary,
    'new_cpp_assertions': cpp_total,
    'benchmark_audit_checks': audit_checks,
    'timed_calls': ANALYSIS['timed_calls'],
    'outcomes': ANALYSIS['outcomes'],
    'conclusion': 'Parallel execution is functioning; the frozen experiment shows mixed scaling and no general improvement over the unchanged sequential GMP implementation.',
}, indent=2) + '\n')
print(f'Wrote {output} ({output.stat().st_size:,} bytes)')
