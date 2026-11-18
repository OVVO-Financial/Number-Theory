"""Create the unified-factorizer report from saved measurements, without searches."""
from __future__ import annotations
import collections
import html
import json
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results"
DATA=OUT/"fused"
ROWS=json.loads((DATA/"fresh/benchmark-summary.json").read_text())
REPLAY=json.loads((DATA/"original-replay/benchmark-summary.json").read_text())
STRESS=json.loads((DATA/"stress/benchmark-summary.json").read_text())
RAW=[json.loads(s) for s in (DATA/"fresh/benchmark-runs.jsonl").read_text().splitlines()]
PROTOCOL=json.loads((DATA/"protocol.json").read_text())
ENV=json.loads((DATA/"fresh/environment.json").read_text())
TRACE=json.loads((DATA/"interaction-trace.json").read_text())
LABELS={
    "fused":"Unified algorithm", "fused_no_vf":"Unified, VF probes disabled",
    "fused_no_feedback":"Unified, bound feedback disabled",
    "fused_unshared":"Unified, separate polynomials",
    "adaptive_product":"Adaptive product search alone", "fermat160":"Same sieved Fermat kernel alone",
    "vf_fast":"Earlier fast VF", "product_trial":"Earlier ordinary product baseline",
    "native_pollard_brent":"Brent rho (C++)", "flint_native":"FLINT full factorization",
}
MAIN=["fused","fermat160","adaptive_product","vf_fast","product_trial","native_pollard_brent","flint_native"]


def select(method,family=None,bits=None,composites=True):
    return [r for r in ROWS if r["algorithm"]==method and (not composites or r["family"]!="prime")
            and (family is None or r["family"]==family) and (bits is None or r["target_bits"]==bits)]


def median_ms(rows):
    values=[r["median_seconds"]*1000 for r in rows if r["solved"]]
    return statistics.median(values) if values else None


def cell(rows):
    value=median_ms(rows)
    return ("—" if value is None else f"{value:,.4f} ms")+f'<small>{sum(r["solved"] for r in rows)}/{len(rows)} completed</small>'


def table(headers,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+str(v)+'</th>' for v in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def paired(reference, family=None, active_only=False):
    a={r["case_id"]:r for r in select(reference,family)}
    b={r["case_id"]:r for r in select("fused",family)}
    keys=[k for k in a if a[k]["solved"] and b[k]["solved"] and (not active_only or b[k]["counters"].get("mixed_batches",0)>0)]
    ratios=[a[k]["median_seconds"]/b[k]["median_seconds"] for k in keys]
    return dict(reference=reference,family=family or "all composites",active_only=active_only,
                pairs=len(keys),fused_faster=sum(v>1 for v in ratios),
                median_ratio=statistics.median(ratios),geometric_ratio=statistics.geometric_mean(ratios),
                fused_extra=sum(b[k]["solved"] and not a[k]["solved"] for k in a),
                reference_extra=sum(a[k]["solved"] and not b[k]["solved"] for k in a))


statuses={}
for row in RAW:
    assert row["status"] in {"factor","prime","budget"}
    if row["factor"]:
        assert 1 < row["factor"] < row["n"] and row["n"]%row["factor"]==0
    statuses.setdefault((row["case_id"],row["algorithm"]),set()).add(row["status"])
mixed=[dict(case_id=k[0],algorithm=k[1],statuses=sorted(v)) for k,v in statuses.items() if len(v)>1]
full=select("fused")
paired_rows=[paired(m) for m in LABELS if m!="fused"]
paired_rows += [paired(m,active_only=True) for m in ("fused_no_vf","fused_no_feedback","fused_unshared")]
for family in ("balanced","close","unbalanced"):
    paired_rows.extend(paired(m,family) for m in ("fermat160","product_trial","vf_fast"))
phases=collections.Counter(r["details"].get("phase","precheck") for r in full)
builds={m:sum(r["counters"].get("polynomials_built",0) for r in select(m)) for m in LABELS}
analysis=dict(composite_successes={m:sum(r["solved"] for r in select(m)) for m in LABELS},
              replay_composite_successes=sum(r["solved"] for r in REPLAY if r["family"]!="prime"),
              phases=dict(phases),polynomial_build_totals=builds,paired=paired_rows,mixed_repeat_cases=mixed,
              fresh_rows=len(ROWS),fresh_raw_calls=len(RAW),
              mixed_batch_inputs=sum(r["counters"].get("mixed_batches",0)>0 for r in full),
              max_fused_wall_seconds=max(r["wall_seconds"] for r in RAW if r["algorithm"]=="fused"))
(DATA/"analysis.json").write_text(json.dumps(analysis,indent=2)+'\n')

outcomes=table(["Method","Composites completed","Prime controls correct"],
    [[LABELS[m],f'{sum(r["solved"] for r in select(m))}/75',
      f'{sum(r["solved"] for r in select(m,"prime",composites=False))}/8'] for m in LABELS])
family_counts=table(["Family","Unified","Sieved Fermat alone","Adaptive product alone","Earlier product baseline"],
    [[family]+[f'{sum(r["solved"] for r in select(m,family))}/{len(select(m,family))}'
               for m in ("fused","fermat160","adaptive_product","product_trial")]
     for family in ("balanced","close","unbalanced")])
timings=''
for family in ("balanced","close","unbalanced"):
    timings+='<h3>'+{"balanced":"Balanced, nonclose factors","close":"Close factors","unbalanced":"Unbalanced factors"}[family]+'</h3>'
    timings+=table(["Bits"]+[LABELS[m] for m in MAIN],
        [[bits]+[cell(select(m,family,bits)) for m in MAIN] for bits in (24,32,40,48,56,64)])
pairs_table=table(["Comparator / unified","Subset","Matched successes","Median ratio","Geometric ratio","Unified faster","Additional completions"],
    [[LABELS[p["reference"]],p["family"]+("; mixed-product inputs" if p["active_only"] else ""),p["pairs"],
      f'{p["median_ratio"]:.3f}×',f'{p["geometric_ratio"]:.3f}×',f'{p["fused_faster"]}/{p["pairs"]}',
      f'Unified +{p["fused_extra"]}; comparator +{p["reference_extra"]}'] for p in paired_rows])
trace_table=table(["Event","First possible divisor","Upper divisor bound","Fermat endpoint cap"],
    [[e["event"],f'{e["lo"]:,}',f'{e["hi"]:,}',f'{e["fermat_cap"]:,}' if e["fermat_cap"] else "—"]
     for e in TRACE["details"]["trace"]])
stress_table=table(["Case","Method","Outcome","Observed time","Factor"],
    [[r["case_id"],LABELS.get(r["algorithm"],"Brent rho (Python)"),r["status"],
      f'{r["median_seconds"]*1000:,.3f} ms',r["factor"] or "—"] for r in STRESS])
test_text=(OUT/"fused-full-test-results.txt").read_text()
test_summary=next(line for line in reversed(test_text.splitlines()) if "passed" in line)
validation=json.loads((OUT/"fusion-native-sqrt-validation.json").read_text())

document=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Unified factorizer — implementation and measured results</title>
<style>
:root{{color-scheme:light;--ink:#16324c;--muted:#52677b;--line:#dae4ef}}*{{box-sizing:border-box}}
body{{margin:0;background:#f5f8fc;color:var(--ink);font:16px/1.65 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1360px;margin:auto;padding:32px 28px 64px}}header{{padding:35px;background:#103b4c;color:white;border-radius:15px}}
h1{{font-size:35px;line-height:1.18;margin:12px 0 20px}}header p{{font-size:19px;color:#e4f0f7}}.tag{{font-size:12px;letter-spacing:.1em;text-transform:uppercase;color:#aad9d7}}
h2{{font-size:25px;margin:24px 0 12px}}h3{{font-size:19px}}section{{background:white;border:1px solid var(--line);border-radius:12px;padding:10px 24px 24px;margin:24px 0}}
p{{max-width:1120px}}.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;font-size:14px;margin:16px 0}}td,th{{text-align:left;vertical-align:top;padding:9px 11px;border-bottom:1px solid var(--line)}}th{{background:#eaf2f7}}small{{display:block;color:var(--muted);font-size:12px}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:22px 0}}.card{{background:white;border:1px solid var(--line);border-radius:10px;padding:18px}}.card b{{display:block;font-size:29px;color:#07776d}}.card span{{font-size:14px;color:var(--muted)}}
.note{{border-left:4px solid #c18a25;background:#fff8e9;padding:14px 18px}}.good{{border-left:4px solid #168075;background:#edf8f5;padding:14px 18px}}.math{{font:19px/1.8 Georgia,serif;background:#eef5f9;padding:14px 20px;overflow-x:auto}}
pre{{background:#112a40;color:#e7f0f9;padding:19px;border-radius:9px;overflow:auto;font-size:13px}}code{{font-family:ui-monospace,SFMono-Regular,monospace}}summary{{cursor:pointer;font-weight:650}}a{{color:#075e99}}.fine{{font-size:13px;color:var(--muted)}}
@media(max-width:700px){{main{{padding:16px}}header{{padding:23px}}h1{{font-size:28px}}section{{padding:8px 14px 18px}}.cards{{grid-template-columns:repeat(2,1fr)}}}}
</style></head><body><main>
<header><div class="tag">Number Theory and Factorization · 17 November 2026</div>
<h1>One unified factorizer, implemented and tested</h1>
<p>The complete Fermat sieve, iterated-average witnesses, and divisor products now share persistent search state, reciprocal factor bounds, a common polynomial engine, aggregate GCD extraction, and one global deadline.</p></header>
<div class="cards"><div class="card"><b>75 / 75</b><span>fresh composites factored within one second</span></div><div class="card"><b>75 / 75</b><span>original composites factored on replay</span></div><div class="card"><b>5,808</b><span>correctness tests passed</span></div><div class="card"><b>112 vs 220</b><span>polynomial builds: shared versus unshared</span></div></div>

<section><h2>The implemented synthesis</h2>
<p>The entry point is <code>vf_factor.fused.factor</code>. It receives N and resource settings only. Shared checks handle small factors, squares, and deterministic primality below 2⁶⁴ once. It then keeps one possible-smaller-factor interval, a persistent complete-sieve Fermat cursor, one cursor for each retained VF progression, and reusable polynomial engines by degree.</p>
<p>Wheel 20 and QR64 are intersected into a midpoint wheel of period 160. The source's complete paired residue condition is preserved. Congruence-rejected midpoints are skipped directly; each tested point uses an exact square check. State is retained across time slices.</p>
<p>The ordinary divisor wheel uses the eight classes coprime to 30. Its block requests and step-five VF requests enter one normalized multipoint batch. A usual batch needs one aggregate GCD. GCD equal to N triggers recovery that continues past trivial VF-only hits. No unfinished divisor strip is marked certified.</p>
<p>One deadline covers prechecks, Fermat work, polynomial construction and evaluation, and recovery. The optional work cap counts scalar divisor positions, represented polynomial positions, and actual Fermat square tests. Recovery revisits already represented terms. The <code>--full</code> command shares the same deadline and work allowance across recursive cofactor splits.</p></section>

<section><h2>How the searches constrain one another</h2>
<p>Let ℓ be the first integer beyond the certified divisor prefix, and let B be the next admissible untested Fermat midpoint. An unknown smaller divisor d has midpoint R(d) = (d + N/d)/2, which decreases for d ≤ √N.</p>
<div class="math">Remaining smaller factor: &nbsp; ℓ ≤ d ≤ U<br>
Fermat → divisor bound: &nbsp; U = B − ⌈√(B² − N)⌉<br>
Divisor prefix → Fermat cap: &nbsp; R<sub>max</sub> = ⌊(ℓ² + N)/(2ℓ)⌋</div>
<p>Every midpoint below B has failed an exact square test or a necessary congruence. An actual factor midpoint must therefore be at least B, giving the upper divisor bound. Conversely, after excluding divisors below ℓ, no remaining factor can have a midpoint above R<sub>max</sub>. When ℓ exceeds U, the remaining interval is empty.</p>
<p>The implementation uses U when sizing and clipping divisor strips and uses R<sub>max</sub> on each Fermat advance. VF witness misses do not change either factor bound; the source's half-gap restrictions are never misapplied to arbitrary GCD offsets.</p>
<h3>An actual execution trace</h3><p>For N = 978508015703, the following trace records both directions of information exchange. The final factor is 752867, with cofactor 1299709. Trace timing is illustrative; benchmark medians are reported separately.</p>
{trace_table}</section>

<section><h2>One polynomial for both searches</h2>
<div class="math">P<sub>m</sub>(X) = ∏<sub>j=0</sub><sup>m−1</sup>(X + j)<br>
∏<sub>j=0</sub><sup>m−1</sup>(α + sj) ≡ s<sup>m</sup>P<sub>m</sub>(αs<sup>−1</sup>) &nbsp; (mod N)<br>
gcd(∏(α + sj), N) = gcd(P<sub>m</sub>(αs<sup>−1</sup>), N), &nbsp; when gcd(s,N)=1.</div>
<p>After the small-factor checks, steps 30 and 5 are units modulo N. Divisor and VF blocks therefore differ only in their normalized evaluation points. Their common polynomial is built once per degree and reused. Monic remainder trees remain valid over composite moduli.</p>
<p>In the fresh validation, 40 inputs entered mixed VF/divisor batches. Their representative runs used 112 polynomial constructions in total, versus 220 for the matched separate-polynomial control. On those 40 inputs, the unified shared version was faster on 27, with a median paired speedup of about 1.114×. Timing-based scheduling also changes the amount of subsequent search work, so this is an implementation comparison rather than a pure polynomial microbenchmark.</p></section>

<section><h2>Fresh validation results</h2>
<p>The final validation uses seed 926483, separate from tuning seed 826391 and the first engineering run's seed 926477. Generated semiprimes are new; published examples, anchor regressions, and prime cubes remain fixed controls. Each search gets one second on one pinned logical CPU. FLINT uses one thread. Successful cases receive three timed repeats; a failed repeat ends further repetitions and makes the aggregate a failure.</p>
{outcomes}
{family_counts}
<p class="good">The unified algorithm completed every fresh composite and all eight prime controls, and also completed all 75 original composites plus eight controls on replay. The standalone sieved Fermat kernel misses nine fresh unbalanced inputs; the earlier product baseline misses seven balanced or close 64-bit inputs. The combined algorithm covers both groups.</p>
<p>Fresh runs: {len(ROWS):,} case/method summaries and {len(RAW):,} timed calls. There are no invalid factors, false prime claims, or algorithm errors. One adaptive-product control case had mixed completion across repeats and is conservatively counted as a failure. All unified-algorithm repeats completed consistently.</p></section>

<section><h2>Timing by factor geometry</h2>
<p>Each cell gives the median over completed case medians and the completion count. Partial-success times are not timings for the entire group. The standalone Fermat control uses the same optimized kernel and wrapper as the unified method, so the shared schedule is compared against the same arithmetic.</p>
{timings}
<p>On the fresh 64-bit balanced group, the unified median is 214.042 ms, standalone sieved Fermat is 130.449 ms, native Brent rho is 1.951 ms, and FLINT is 0.384 ms. The unified method spends time covering the small-factor side as well; it does not win every individual geometry. The earlier VF and product methods complete none of these four inputs under the budget.</p></section>

<section><h2>What each component contributed</h2>
<p>The successful unified phases were {html.escape(str(dict(phases)))}. All 24 balanced cases and all 18 close-factor cases finished through Fermat. Seventeen of the 18 unbalanced inputs finished through the mixed-product engine; the other finished through Fermat. A mixed-product result does not by itself identify whether a VF block or an ordinary-divisor block supplied the factor.</p>
<p>The no-VF, no-feedback, and separate-polynomial controls also complete all 75 fresh composites. VF probes do not demonstrate a net speed advantage here: the all-composite median comparator/unified ratio is about 0.917. The no-feedback median ratio is about 0.929. Ratios below one favor the ablation.</p>
<p>The bounds are exact and update in actual runs, but this experiment does not establish a net work or speed saving attributable to bound feedback. The final clipping counter is zero; that counter alone does not measure earlier degree/strip sizing. The evidence supports a useful unified search and actual arithmetic sharing, while leaving the extra benefit of VF probes and reciprocal bounds unproved on this corpus.</p>
<p class="note">This is a working integrated algorithm, not a claim that each constituent heuristic is independently beneficial. Native Brent and FLINT remain faster general-purpose references on most inputs. No new complexity exponent follows from these finite experiments.</p>
<details><summary>Complete paired comparisons</summary><p>Ratios are comparator time / unified time, computed per matched successful input before aggregation. Values above one favor the unified method. Additional completions are recorded separately.</p>{pairs_table}</details></section>

<section><h2>Arithmetic and budget verification</h2>
<p><strong>{html.escape(test_summary)}</strong></p>
<p>The new tests check resumable coverage, native/Python cursor agreement, temporary Fermat caps, all small odd inputs, normalized mixed-step products over composite rings, negative starts and steps, trivial GCD-N recovery, both bound directions against actual divisors, trace invariants, arbitrary-precision fallback, and one shared work limit. Every returned factor is checked by exact division.</p>
<p>The first engineering run identified the integer square-root loop as a bottleneck. The native uint64 routine now uses a hardware estimate only to seed exact integer correction loops. The final root satisfies the integer floor-root inequalities, and all squares are compared exactly. Arguments larger than uint64 keep the integer-Newton fallback. An independent 45,389-case comparison against Python's integer square root passed; the entire suite passed again afterward.</p>
<p>The first engineering data and its source snapshot are retained under <code>results/fused-development-integer-root</code>. Scheduling settings were selected on the earlier pilot and not retuned after the square-root optimization. The final fresh seed was chosen separately.</p>
<details><summary>Independent square-root validation record</summary><pre>{html.escape(json.dumps(validation,indent=2))}</pre></details></section>

<section><h2>Arbitrary-precision stress tests</h2>
<p>Inputs above 64 bits use the exact Python Fermat cursor and FLINT's arbitrary-precision polynomial backend. The small-factor 80-, 96-, and 128-bit cases complete. The balanced 80-bit case remains unresolved after two seconds, while Python Brent rho completes it. Total input bit length alone is therefore not a measure of demonstrated difficulty.</p>
{stress_table}</section>

<section><h2>Run the single algorithm</h2>
<pre>python3 -m pip install -r requirements.txt
python3 build_fusion.py
python3 -m vf_factor.fused 978508015703 --full --timeout 1

# Optional comparison kernels and the complete test suite:
python3 build_native.py
python3 build_wheel20.py
python3 -m pytest -q

# Reproduce the final comparison in a separate directory:
python3 fused_benchmark.py --outdir results/fused-reproduction
python3 fused_report.py</pre>
<p>After package installation, <code>vf-factor</code> is the same unified entry point; <code>vf-factor-compare</code> retains access to the older individual methods. The report generator reads the bundled <code>results/fused</code> measurements. The reproduction command writes a separate directory. The native kernel requires GCC or Clang with C++17 and unsigned 128-bit integer support. Core factoring remains available through the Python cursor if the native kernel is absent.</p>
<p>Timeouts are cooperative: an in-flight polynomial operation can overshoot slightly. <code>status: budget</code> means unresolved. The fixed degree cap of 512 bounds resources and should not be read as an implementation of an unlimited asymptotic degree choice.</p>
<details><summary>Environment and frozen protocol</summary><pre>{html.escape(json.dumps(ENV,indent=2))}</pre><pre>{html.escape(json.dumps(PROTOCOL,indent=2))}</pre></details></section>

<section><h2>Sources</h2><ul>
<li>Fred Viole, supplied <em>Complex Space Factorization</em> PDF and <em>Iterated Averages.md</em>.</li>
<li><a href="https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Complete_Fermat_Sieve_Verified.pdf">The Complete Fermat Sieve by Terminal Digit: A corrected and verified reference.</a></li>
<li><a href="https://mathexp.eu/bostan/publications/BoGaSc04.pdf">Bostan, Gaudry and Schost: polynomial product and multipoint evaluation methods.</a></li>
<li><a href="https://python-flint.readthedocs.io/en/latest/fmpz_mod_poly.html">python-flint polynomial documentation.</a></li>
</ul><p>The reciprocal-bound and normalization identities are derived explicitly above and independently checked in the supplied tests.</p></section>
</main></body></html>'''
(OUT/"fused_report.html").write_text(document)
print(OUT/"fused_report.html")

