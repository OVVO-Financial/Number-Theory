"""Build the self-contained report from saved benchmark data; no algorithms run."""
from __future__ import annotations
import base64
import collections
import html
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
ROWS = json.loads((OUT / "main/benchmark-summary.json").read_text())
CASES = json.loads((OUT / "main/benchmark-corpus.json").read_text())
ENV = json.loads((OUT / "main/environment.json").read_text())
STRESS = json.loads((OUT / "stress/benchmark-summary.json").read_text())

LABELS = {
    "published_vf":"Original VF (Python)", "native_published_vf":"Original VF (C++)",
    "vf_fast":"Fast VF", "vf_fast_multi":"Fast VF, all levels", "vf_hybrid":"Fermat + fast VF",
    "vf_first":"Fast VF, level 1", "product_trial":"Ordinary product search",
    "trial_division":"Trial wheel (Python)", "native_trial_division":"Trial wheel (C++)",
    "fermat":"Fermat (Python)", "native_fermat":"Fermat (C++)",
    "pollard_brent":"Brent rho (Python)", "native_pollard_brent":"Brent rho (C++)",
    "flint_native":"FLINT full factorization",
}
ORDER = list(LABELS)


def group(method, family=None, bits=None, composites=True):
    return [r for r in ROWS if r["algorithm"] == method
            and (not composites or r["family"] != "prime")
            and (family is None or r["family"] == family)
            and (bits is None or r["target_bits"] == bits)]


def median_ms(rows):
    values = [r["median_seconds"] * 1000 for r in rows if r["solved"]]
    return statistics.median(values) if values else None


def cell(rows):
    successes = sum(r["solved"] for r in rows)
    value = median_ms(rows)
    text = "—" if value is None else f"{value:,.3f} ms"
    return f'{text}<small>{successes}/{len(rows)} completed</small>'


def paired(reference, contender, family=None):
    left = {r["case_id"]:r for r in group(reference, family)}
    right = {r["case_id"]:r for r in group(contender, family)}
    common = [k for k in left if left[k]["solved"] and right[k]["solved"]]
    ratios = [left[k]["median_seconds"]/right[k]["median_seconds"] for k in common]
    return dict(completed_pairs=len(common), median_ratio=statistics.median(ratios),
                geometric_ratio=statistics.geometric_mean(ratios),
                contender_faster=sum(v>1 for v in ratios),
                contender_extra=sum(right[k]["solved"] and not left[k]["solved"] for k in left),
                reference_extra=sum(left[k]["solved"] and not right[k]["solved"] for k in left))


plt.rcParams.update({"font.size":10.5,"axes.spines.top":False,"axes.spines.right":False,
                     "font.family":"DejaVu Sans"})
fig, (ax, bar) = plt.subplots(1,2,figsize=(14.7,5.7),gridspec_kw={"width_ratios":[1.25,1]})
fig.patch.set_facecolor("#f7f9fc")
colors = {"native_published_vf":"#697586", "vf_fast":"#176bca", "vf_first":"#8a4ec4",
          "product_trial":"#d87920", "native_pollard_brent":"#0d7c65", "flint_native":"#cc4267"}
bits_list=[24,32,40,48,56,64]
for offset,(method,color) in enumerate(colors.items()):
    values=[]
    for bits in bits_list:
        rows=group(method,"balanced",bits)
        completed=all(r["solved"] for r in rows)
        values.append(median_ms(rows) if completed else np.nan)
        if not completed:
            ax.scatter(bits+(offset-2.5)*.18,1000,marker="x",color=color,s=55,linewidths=2,zorder=5)
    ax.plot(bits_list,values,"o-",color=color,label=LABELS[method],linewidth=2,markersize=4)
ax.axhline(1000,color="#9ba8b9",linestyle=":",linewidth=1)
ax.set_yscale("log")
ax.set_ylim(.002,2000)
ax.set_xticks(bits_list)
ax.set_xlabel("Input size (bits); four nonclose balanced semiprimes per size")
ax.set_ylabel("Median time (milliseconds, logarithmic scale)")
ax.set_title("Larger balanced inputs benefit from batching",loc="left",fontweight="bold",pad=14)
ax.grid(axis="y",alpha=.18)
ax.legend(loc="upper left",fontsize=8.4,frameon=False)

bar_methods=["native_published_vf","vf_fast","vf_fast_multi","vf_hybrid","vf_first",
             "product_trial","native_pollard_brent","flint_native"]
counts=[sum(r["solved"] for r in group(m)) for m in bar_methods]
y=np.arange(len(bar_methods))
bar.barh(y,[75]*len(y),color="#e5ebf2",height=.6)
bar.barh(y,counts,color=[colors.get(m,"#77a7dd") for m in bar_methods],height=.6)
for idx,count in enumerate(counts): bar.text(count+.7,idx,f"{count}/75",va="center",fontsize=10)
bar.set_yticks(y,[LABELS[m] for m in bar_methods])
bar.invert_yaxis();bar.set_xlim(0,84);bar.set_xticks([0,25,50,75])
bar.set_xlabel("Composite inputs factored within the one-second budget")
bar.set_title("Established methods completed the entire corpus",loc="left",fontweight="bold",pad=14)
bar.grid(axis="x",alpha=.12);bar.set_axisbelow(True)
fig.suptitle("Iterated-average factorization: measured capability",x=.055,ha="left",fontsize=17,fontweight="bold")
fig.text(.055,.015,"Left: lines are shown only when all four cases completed; × at 1,000 ms marks a group with at least one timeout.\nRight: includes all 75 composites, including easy cases; eight prime controls are excluded. Times include per-input setup.",fontsize=9,color="#4a586a")
fig.subplots_adjust(left=.065,right=.97,top=.84,bottom=.20,wspace=.48)
plot=OUT/"benchmark_comparison.png"
fig.savefig(plot,dpi=160,facecolor=fig.get_facecolor())
plt.close(fig)

paired_rows=[]
for reference,contender,family in [
    ("native_published_vf","vf_fast","balanced"),
    ("native_published_vf","vf_fast",None),
    ("vf_first","vf_fast","balanced"),
    ("vf_first","vf_fast",None),
    ("product_trial","vf_fast","balanced"),
    ("product_trial","vf_fast","unbalanced"),
    ("native_pollard_brent","vf_fast","balanced")]:
    p=paired(reference,contender,family)
    paired_rows.append(dict(reference=reference,contender=contender,family=family or "all composites",**p))
(OUT/"analysis.json").write_text(json.dumps(dict(paired=paired_rows,
    composite_successes={m:sum(r["solved"] for r in group(m)) for m in ORDER}),indent=2))

success_table="".join(f'<tr><td>{LABELS[m]}</td><td>{sum(r["solved"] for r in group(m))}/75</td><td>{sum(not r["solved"] for r in group(m))}</td><td>{sum(r["solved"] for r in group(m,"prime",composites=False))}/8</td></tr>' for m in ORDER)
balanced_methods=["native_published_vf","vf_fast","vf_first","product_trial","native_pollard_brent","flint_native"]
balanced_table="".join('<tr><td>'+str(bits)+'</td>'+''.join('<td>'+cell(group(m,"balanced",bits))+'</td>' for m in balanced_methods)+'</tr>' for bits in bits_list)
pair_table="".join(f'<tr><td>{html.escape(LABELS[p["reference"]])}</td><td>{p["family"]}</td><td>{p["completed_pairs"]}</td><td>{p["median_ratio"]:.3f}×</td><td>{p["geometric_ratio"]:.3f}×</td><td>{p["contender_faster"]}/{p["completed_pairs"]}</td><td>VF +{p["contender_extra"]}; comparator +{p["reference_extra"]}</td></tr>' for p in paired_rows)
stress_table="".join(f'<tr><td>{r["case_id"]}</td><td>{LABELS[r["algorithm"]]}</td><td>{"factor" if r["solved"] else "budget exhausted"}</td><td>{r["median_seconds"]*1000:.3f} ms</td><td>{r["factor"] or "—"}</td></tr>' for r in STRESS)
families=collections.Counter(c["family"] for c in CASES)
phase=collections.Counter(r["details"].get("phase","precheck") for r in group("vf_hybrid") if r["solved"])
source_case=[r for r in ROWS if r["case_id"]=="published-978508015703"]
example_table="".join(f'<tr><td>{LABELS[r["algorithm"]]}</td><td>{r["factor"]}</td><td>{r["median_seconds"]*1000:.3f} ms</td><td>{r["counters"].get("gcd_calls","—")}</td></tr>' for r in source_case if r["algorithm"] in ["published_vf","native_published_vf","vf_fast","product_trial","native_pollard_brent","flint_native"])
png=base64.b64encode(plot.read_bytes()).decode()
tests=(OUT/"test-results.txt").read_text()
test_summary=next((line for line in reversed(tests.splitlines()) if "passed" in line),"See test-results.txt")
matrix=[]
for case in CASES:
    cells=[]
    for m in balanced_methods+["vf_hybrid"]:
        row=next(r for r in ROWS if r["case_id"]==case["case_id"] and r["algorithm"]==m)
        cells.append(f'<td>{row["median_seconds"]*1000:.3f}' + (' ms' if row["solved"] else ' ms · TIMEOUT')+'</td>')
    matrix.append(f'<tr><td>{case["case_id"]}<small>N={case["n"]}</small></td>'+''.join(cells)+'</tr>')

document=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>VF Factor Lab — build and benchmark report</title>
<style>
:root{{color-scheme:light;--ink:#172b43;--muted:#58687c;--line:#dce5ef;--blue:#176bca;}}
*{{box-sizing:border-box}}body{{margin:0;background:#f5f8fc;color:var(--ink);font:16px/1.6 system-ui,-apple-system,Segoe UI,sans-serif}}main{{max-width:1240px;margin:auto;padding:38px 28px 70px}}header{{padding:38px;background:#102d4d;color:white;border-radius:16px}}h1{{font-size:36px;line-height:1.15;margin:8px 0 18px}}h2{{font-size:25px;margin:40px 0 14px}}h3{{font-size:19px;margin-top:25px}}p{{max-width:1000px}}.eyebrow{{font-size:12px;letter-spacing:.12em;text-transform:uppercase;color:#a9d4ff}}.lead{{font-size:19px;color:#e2edfa;max-width:1050px}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:24px 0}}.card{{background:white;border:1px solid var(--line);border-radius:12px;padding:19px}}.card b{{display:block;font-size:30px;color:var(--blue)}}.card span{{color:var(--muted);font-size:14px}}section{{background:white;border:1px solid var(--line);border-radius:13px;padding:10px 25px 25px;margin:22px 0}}section h2{{margin-top:18px}}.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;font-size:14px;margin:14px 0}}th{{text-align:left;background:#eef3f9;font-weight:650}}td,th{{padding:10px 12px;border-bottom:1px solid var(--line);vertical-align:top}}small{{display:block;color:var(--muted);font-size:12px}}.note{{padding:15px 18px;border-left:4px solid #e1a947;background:#fff8e9}}.good{{border-left-color:#176bca;background:#eef6ff}}pre{{background:#10243b;color:#e4effb;padding:20px;overflow:auto;border-radius:8px;font-size:14px}}code{{font-family:ui-monospace,SFMono-Regular,monospace}}img{{max-width:100%;height:auto}}summary{{cursor:pointer;font-weight:650}}a{{color:#075ca7}}.fine{{font-size:13px;color:var(--muted)}}@media(max-width:700px){{main{{padding:16px}}header{{padding:24px}}h1{{font-size:29px}}.cards{{grid-template-columns:repeat(2,1fr)}}section{{padding:8px 15px 20px}}}}
</style></head><body><main>
<header><div class="eyebrow">Number Theory and Factorization · 17 November 2026</div><h1>Iterated averages, built and tested</h1><p class="lead">Fast polynomial block searches improve the original VF routine on larger balanced examples. The current implementation does not outperform established Brent rho or FLINT, and the published choice of averaging level shows no consistent advantage over the simpler controls.</p></header>
<div class="cards"><div class="card"><b>83</b><span>fixed main inputs: 75 composites + 8 primes</span></div><div class="card"><b>14</b><span>implementations and controlled variants</span></div><div class="card"><b>70 / 75</b><span>fast VF completions within one second</span></div><div class="card"><b>0</b><span>invalid factors or false prime claims</span></div></div>
<img src="data:image/png;base64,{png}" alt="Balanced-input timing curves and composite completion counts">
<section><h2>What was actually implemented</h2><p>The project includes an exact reproduction of the later Julia <code>VF(n,5)</code> routine, scalar and native C++ baselines, a reusable fast polynomial-product engine, one-level and multiple-level schedules, a first-level ablation, and an explicit Fermat/VF hybrid. All factor outputs are checked by exact division.</p><p>The two published examples reproduce the same selected levels and search positions. For <code>978508015703</code>, the scalar routine selects level 6 and returns 752867 at distance 50745 using GCD argument 963218792667. The later Julia routine also succeeds on the midpoint-parity regressions; it is distinct from importing the PDF's Fermat digit classes unchanged into a fixed-seed offset.</p><div class="scroll"><table><thead><tr><th>Method on 978508015703</th><th>Factor</th><th>Median time</th><th>GCD calls</th></tr></thead><tbody>{example_table}</tbody></table></div><p class="fine">Times include setup. GCD counts are not interchangeable operation counts; polynomial products represent many terms collectively.</p></section>
<section><h2>Main outcomes</h2><div class="scroll"><table><thead><tr><th>Method</th><th>Composites completed</th><th>Budget failures</th><th>Prime controls correct</th></tr></thead><tbody>{success_table}</tbody></table></div><p class="note">A one-second limit is cooperative. A polynomial operation already running can finish after the nominal deadline; these runs remain budget failures. No failed run is counted as a fast factorization.</p></section>
<section><h2>Balanced semiprimes: complete timing picture</h2><p>Four nonclose balanced semiprimes were generated at each exact bit size. Each cell shows the median among completed cases and explicitly identifies the completion count. Partial-success medians must not be compared as if they describe all four cases.</p><div class="scroll"><table><thead><tr><th>Bits</th>{''.join('<th>'+LABELS[m]+'</th>' for m in balanced_methods)}</tr></thead><tbody>{balanced_table}</tbody></table></div><p>At 48 bits, all four cases completed for both original C++ VF and fast VF. Their family medians were 131.507 ms and 18.210 ms, respectively. The median of the four paired speedup ratios was 8.651×. At 56 bits fast VF completed all four while the original C++ routine completed two. At 64 bits fast VF and ordinary product search completed none within the fixed limit; Brent and FLINT completed all four.</p></section>
<section><h2>Does the averaging level add value?</h2><p>The cleanest schedule comparison is <code>vf_fast</code> against <code>vf_first</code>: the same polynomial engine, prefix and stream rules, but the latter always uses level one. Level one completed 72/75 composites, while the selected-level version completed 70/75. On the 20 balanced inputs both completed, selected-level VF was faster on 11; level one additionally completed two 64-bit balanced cases. These data do not establish a robust benefit from the published level-selection rule.</p><p>Ordinary product search beat selected-level VF on 16 of 20 jointly completed balanced cases. VF was faster on 13 of 18 unbalanced cases, but that comparison also changes adaptive block growth versus larger upfront polynomial construction. The unbalanced gain cannot be assigned solely to the averages.</p><div class="scroll"><table><thead><tr><th>Comparator against fast VF</th><th>Family</th><th>Matched successes</th><th>Median comparator/VF</th><th>Geometric ratio</th><th>VF faster</th><th>Unmatched completions</th></tr></thead><tbody>{pair_table}</tbody></table></div><p class="fine">Ratios above 1 favor fast VF; below 1 favor the comparator. Every ratio uses only inputs both completed. Completion differences are displayed separately.</p></section>
<section><h2>Hybrid benefits and attribution</h2><p>The hybrid completed 71/75 composites. All 18 close-semiprime successes were found by its preliminary Fermat pass, so their speed is credited to Fermat. The successful hybrid phases were: {html.escape(str(dict(phase)))}. All 20 completed balanced fast-VF cases used polynomial searches; their improvement over scalar VF is not a hidden Fermat fallback.</p><p>The shared square and small-factor prechecks are identical across methods. Python/native implementations are labelled separately to avoid presenting a language-runtime difference as a mathematical gain. Native Brent uses a different fixed RNG stream from Python Brent; its timing comparison includes both implementation and trajectory effects.</p></section>
<section><h2>Beyond 64-bit arithmetic</h2><p>A separate stress set used 80-, 96- and 128-bit inputs with small factors, plus one balanced 80-bit semiprime. This validates actual factor searches with the arbitrary-precision polynomial backend. A large total input with a small factor is easier than a balanced input of the same size.</p><div class="scroll"><table><thead><tr><th>Case</th><th>Method</th><th>Outcome</th><th>Observed time</th><th>Factor</th></tr></thead><tbody>{stress_table}</tbody></table></div><p>The 80-bit balanced case remained unresolved by the VF/product variants after the two-second budget; Python Brent completed it. These stress tests had one timing run per method and are not used to infer general scaling.</p></section>
<section><h2>Why the accelerated search is correct</h2><p>For a progression <code>alpha + s*j</code>, the block polynomial <code>F(X)=product(X+s*i, i=0..m-1)</code> evaluated at <code>alpha+s*m*r</code> is exactly the product of block r. The implementation builds products and evaluates many points using polynomial trees, rather than multiplying all represented scalar terms individually.</p><p>Monic polynomial divisors permit the computation over composite moduli. A proper block GCD returns a factor. A block GCD equal to N is inspected, and searching continues if it contains only trivial hits.</p><p>After removing 5, any step-five progression encounters a p-divisible term within p positions. If that term is divisible by all of N, the term p positions later has GCD p because its difference is 5p. Thus 2*floor(sqrt(N))+2 indices provide complete coverage for composite N. Finite time budgets do not inherit that completion guarantee; their failures are explicitly reported.</p><p>This is established Strassen-style collective product-search machinery applied to VF progressions. It supplies a real acceleration but does not establish a new exponent. The fixed degree cap 4096 controls resources in this prototype.</p></section>
<section><h2>Verification and reproducibility</h2><p><strong>{html.escape(test_summary)}</strong></p><p>The suite checks small integers exhaustively, composite-ring polynomial products against independent direct multiplication, negative starts, GCD=N recovery, tails, large moduli, source examples, midpoint regressions, prime powers, and budget behavior. Native checks additionally cover exact source search order and uint64 square-root boundaries. Main benchmark reconciliation found 1,162 summary rows and 3,292 raw timed runs, with no errors, invalid divisors, false prime claims, or mixed completion outcomes across repeated runs.</p><p>Corpus composition: {html.escape(str(dict(families)))}. Known factors are never passed to a timed algorithm. A separate tuning seed was used for the small degree-cap pilot, and the initial 4096 cap was retained before main timings.</p><pre>python3 -m pip install -r requirements.txt
python3 build_native.py
python3 -m vf_factor.cli 978508015703 --method vf_fast --timeout 10
python3 -m pytest -q
python3 -m vf_factor.benchmark --timeout 1 --repeats 3 --outdir results/reproduction
python3 report.py</pre><p class="fine">The benchmark command writes a separate reproduction directory. The report script regenerates this report from the original saved <code>results/main</code> and <code>results/stress</code> measurements. Core factoring supports Python arbitrary-precision integers; optional C++ comparison kernels are restricted to unsigned 64-bit inputs. Native sources and their build script are included. The archive contains every corpus, raw run, summary, configuration, test and report-generation script.</p><details><summary>Environment and measurement protocol</summary><pre>{html.escape(json.dumps(ENV,indent=2))}</pre></details></section>
<section><h2>Assessment</h2><p><strong>The implementation is a working factorization research prototype, and fast collective searches make the original approach more capable on larger balanced examples.</strong> The current evidence does not support calling it competitive with established general-purpose methods. The useful open engineering question is whether better interval/level selection or more work sharing can improve the fast-product baseline; the present maximum-fraction rule has not demonstrated a consistent advantage.</p><p>Very small reported times vary substantially between repetitions. Differences of a few microseconds should not be treated as stable rankings. The much larger gaps between the polynomial VF variants and native Brent/FLINT remain clear in this experiment.</p></section>
<section><h2>Selected per-input results</h2><details><summary>Show the complete 83-input matrix for seven principal methods</summary><div class="scroll"><table><thead><tr><th>Input</th>{''.join('<th>'+LABELS[m]+'</th>' for m in balanced_methods+['vf_hybrid'])}</tr></thead><tbody>{''.join(matrix)}</tbody></table></div></details></section>
<section><h2>Sources</h2><ul><li>Fred Viole, <em>Complex Space Factorization</em>, supplied PDF, sections 5.1–5.2 and 7.3–7.4.</li><li>Fred Viole, supplied <em>Iterated Averages.md</em>, Julia VF routine.</li><li><a href="https://mathexp.eu/bostan/publications/BoGaSc04.pdf">Bostan, Gaudry and Schost: polynomial product and multipoint-evaluation methods.</a></li><li><a href="https://doi.org/10.1090/mcom/3623">Hittmeir: Strassen-style products and deterministic factorization.</a></li><li><a href="https://python-flint.readthedocs.io/en/latest/fmpz_mod_poly.html">python-flint documentation.</a></li></ul></section>
</main></body></html>'''
(OUT/"benchmark_report.html").write_text(document)
print(OUT/"benchmark_report.html")
print(plot)
