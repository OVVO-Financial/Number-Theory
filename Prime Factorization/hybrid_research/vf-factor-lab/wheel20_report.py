"""Generate the focused wheel-20 report from saved measurements only."""
from __future__ import annotations

import html
import json
import statistics
from collections import Counter
from math import isqrt
from pathlib import Path

from vf_factor.wheel20 import admissible_pairs, midpoint_residues

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
DATA = OUT / "wheel20"
ROWS = json.loads((DATA / "benchmark-summary.json").read_text())
CASES = json.loads((DATA / "benchmark-corpus.json").read_text())
RAW = [json.loads(line) for line in (DATA / "benchmark-runs.jsonl").read_text().splitlines()]
ENV = json.loads((DATA / "environment.json").read_text())
PROTOCOL = json.loads((DATA / "protocol.json").read_text())
LABELS = {
    "native_fermat_parity": "Fermat: parity",
    "native_fermat_wheel20": "Fermat: wheel 20",
    "native_fermat_parity_qr64": "Fermat: parity + QR64",
    "native_fermat_wheel20_qr64": "Fermat: wheel 20 + QR64",
    "product_trial": "Ordinary product search",
    "native_pollard_brent": "Brent rho (C++)",
}
WHEEL = "native_fermat_wheel20_qr64"


def selected(method, family=None, bits=None, composites=True, main=True):
    return [r for r in ROWS if r["algorithm"] == method
            and (not main or r["family"] != "wheel20_supplement")
            and (not composites or r["family"] != "prime")
            and (family is None or r["family"] == family)
            and (bits is None or r["target_bits"] == bits)]


def median_ms(rows):
    values = [r["median_seconds"] * 1000 for r in rows if r["solved"]]
    return statistics.median(values) if values else None


def cell(rows):
    solved = sum(r["solved"] for r in rows)
    ms = median_ms(rows)
    return ("—" if ms is None else f"{ms:,.4f} ms") + f"<small>{solved}/{len(rows)} completed</small>"


def table(headers, rows):
    return '<div class="scroll"><table><thead><tr>' + ''.join('<th>'+h+'</th>' for h in headers) + '</tr></thead><tbody>' + ''.join('<tr>'+''.join('<td>'+str(c)+'</td>' for c in row)+'</tr>' for row in rows) + '</tbody></table></div>'


def paired(reference, contender, family=None, minimum_reference_ms=0):
    a = {r["case_id"]: r for r in selected(reference, family)}
    b = {r["case_id"]: r for r in selected(contender, family)}
    keys = [k for k in a if a[k]["solved"] and b[k]["solved"]
            and a[k]["median_seconds"]*1000 >= minimum_reference_ms]
    ratios = [a[k]["median_seconds"] / b[k]["median_seconds"] for k in keys]
    return dict(reference=reference, contender=contender, family=family or "all composites",
                minimum_reference_ms=minimum_reference_ms, pairs=len(keys),
                contender_faster=sum(v > 1 for v in ratios),
                median_ratio=statistics.median(ratios) if ratios else None,
                geometric_ratio=statistics.geometric_mean(ratios) if ratios else None,
                contender_extra=sum(b[k]["solved"] and not a[k]["solved"] for k in a),
                reference_extra=sum(a[k]["solved"] and not b[k]["solved"] for k in a))


assert len(ROWS) == 86 * 6
assert all(r["status"] in {"factor", "prime", "budget"} for r in RAW)
assert all(1 < r["factor"] < r["n"] and r["n"] % r["factor"] == 0 for r in RAW if r["factor"])
statuses = {}
for row in RAW:
    statuses.setdefault((row["case_id"], row["algorithm"]), set()).add(row["status"])
assert all(len(s) == 1 for s in statuses.values())

success_table = table(["Method", "Original composites", "Original prime controls"],
    [[LABELS[m], f'{sum(r["solved"] for r in selected(m))}/75',
      f'{sum(r["solved"] for r in selected(m,"prime",composites=False))}/8'] for m in LABELS])

family_tables = ''
for family in ("balanced", "close", "unbalanced"):
    family_tables += '<h3>'+{"balanced":"Balanced, nonclose factors", "close":"Close factors", "unbalanced":"Unbalanced factors"}[family]+'</h3>'
    family_tables += table(["Bits"] + list(LABELS.values()),
        [[bits] + [cell(selected(m, family, bits)) for m in LABELS] for bits in (24,32,40,48,56,64)])

pairs = []
for family in ("balanced", "close", "unbalanced", None):
    for contender in ("native_fermat_wheel20", WHEEL):
        pairs.append(paired("product_trial", contender, family))
for ref, contender in (("native_fermat_parity", "native_fermat_wheel20"),
                       ("native_fermat_parity_qr64", WHEEL)):
    pairs.append(paired(ref, contender, "balanced"))
    pairs.append(paired(ref, contender, None, minimum_reference_ms=1))
pair_table = table(["Comparator", "Contender", "Family", "Matched inputs", "Median speedup", "Geometric speedup", "Contender faster", "Unmatched completions"],
    [[LABELS[p["reference"]], LABELS[p["contender"]], p["family"] + ("; comparator ≥1 ms" if p["minimum_reference_ms"] else ""),
      p["pairs"], f'{p["median_ratio"]:.3f}×', f'{p["geometric_ratio"]:.3f}×',
      f'{p["contender_faster"]}/{p["pairs"]}', f'Contender +{p["contender_extra"]}; comparator +{p["reference_extra"]}'] for p in pairs])

lanes = []
for n in (1,3,7,9,11,13,17,19):
    rs = midpoint_residues(n)
    heights = sorted({(4*a*b) % 40 for a,b in admissible_pairs(n)})
    lanes.append([n, ', '.join(f'({a},{b})' for a,b in admissible_pairs(n)),
                  ', '.join(map(str,rs)), f'{len(rs)}/10', ', '.join(map(str, heights))])
lane_table = table(["N mod 20", "Complete (R mod 10, i mod 10) pairs", "R mod 10", "Midpoint density", "4Ri mod 40"], lanes)

supplement = table(["Input", "Method", "Outcome", "Median", "Midpoints visited", "Square roots"],
    [[r["case_id"], LABELS[r["algorithm"]], r["status"], f'{r["median_seconds"]*1000:,.5f} ms',
      r["counters"].get("candidates", "—"), r["counters"].get("square_tests", "—")]
     for r in ROWS if r["family"] == "wheel20_supplement"])

tests=(OUT/"wheel20-test-results.txt").read_text()
test_summary=next(line for line in reversed(tests.splitlines()) if "passed" in line)
stats = {"main_composite_successes": {m:sum(r["solved"] for r in selected(m)) for m in LABELS},
         "paired_comparisons": pairs, "raw_runs":len(RAW), "summary_rows":len(ROWS),
         "source":PROTOCOL["source_url"], "test_summary":test_summary}
(DATA/"analysis.json").write_text(json.dumps(stats,indent=2)+'\n')
source = html.escape(PROTOCOL["source_url"], quote=True)

document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Complete Fermat wheel 20 — measured comparison</title>
<style>
:root{{color-scheme:light;--ink:#16304b;--muted:#53677d;--line:#dce5ef}}*{{box-sizing:border-box}}
body{{margin:0;background:#f5f8fc;color:var(--ink);font:16px/1.65 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1320px;margin:auto;padding:32px 28px 60px}}header{{background:#113956;color:white;padding:34px;border-radius:14px}}
h1{{font-size:34px;line-height:1.2;margin:10px 0 20px}}h2{{font-size:25px;margin:24px 0 12px}}h3{{font-size:19px}}
header p{{font-size:19px;color:#e1eef9}}.tag{{font-size:12px;letter-spacing:.1em;text-transform:uppercase;color:#aed9f6}}
section{{background:white;border:1px solid var(--line);border-radius:12px;padding:10px 24px 24px;margin:24px 0}}
p{{max-width:1080px}}.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:16px 0}}
td,th{{text-align:left;vertical-align:top;padding:9px 12px;border-bottom:1px solid var(--line)}}th{{background:#eaf1f8}}
small{{display:block;color:var(--muted);font-size:12px}}.note{{border-left:4px solid #c99126;background:#fff8e8;padding:14px 18px}}
.good{{border-left:4px solid #178278;background:#edf8f5;padding:14px 18px}}pre{{padding:18px;background:#10283f;color:#e6f0fa;border-radius:8px;overflow:auto}}
code{{font-family:ui-monospace,SFMono-Regular,monospace}}.math{{font-family:Georgia,serif;font-size:19px;padding:12px 18px;background:#f0f5fa;overflow-x:auto}}
a{{color:#075f9f}}summary{{cursor:pointer;font-weight:650}}.fine{{font-size:13px;color:var(--muted)}}
@media(max-width:700px){{main{{padding:16px}}header{{padding:23px}}h1{{font-size:28px}}section{{padding:8px 15px 20px}}}}
</style></head><body><main>
<header><div class="tag">Number Theory and Factorization · 17 November 2026</div>
<h1>The complete Fermat wheel 20 does help — in the right search region</h1>
<p>It beats this ordinary product-search implementation on every jointly completed balanced or close-factor input. The wheel also gives a measurable improvement over matched Fermat controls. Its exact coupled residue restrictions remain periodic: they do not by themselves identify an absolute region containing a factor.</p></header>

<section><h2>The exact specification tested</h2>
<p>The source is <a href="{source}">The Complete Fermat Sieve by Terminal Digit: A corrected and verified reference</a>. Pages 1 and 4–5 distinguish a complete digit-only table from a single parity lane. A table indexed by N mod 10 needs eight pairs for each of 1, 3, 7 and 9. Indexing by N mod 20 selects the appropriate four-pair lane for that particular N. Inputs divisible by 2 or 5 are handled separately.</p>
<p>The implementation generates all pairs satisfying R² − i² ≡ N (mod 20), projects onto eligible midpoint residues, and visits those midpoints through a cyclic gap table. An exact square check for R² − N enforces the full paired condition automatically. The regression 551 = 19 × 29 = 24² − 5² retains its required (4,5) lane.</p>
{lane_table}
<p class="good">This is the complete sieve. No parity lane was silently removed, and no additional half-gap or coupled-height condition is missing from a successful exact square check.</p></section>

<section><h2>What was measured</h2>
<p>Four native implementations share the same C++ kernel and exact square-root routine. Two independent switches select wheel 20 and the QR64 square-residue filter. This isolates the wheel's added value instead of comparing it with a different square-root implementation. The original product-search baseline is unchanged, including its degree cap of 4096. Native Brent rho supplies a further reference.</p>
<p>The original 83 inputs comprise 75 composites and eight prime controls. Three supplemental inputs check 551 and the same-lane example below; these are excluded from original-corpus totals. There is a one-second cooperative budget, three timed repeats for completed cases, fixed randomized method order, one pinned logical CPU, and one FLINT thread. Timings include prechecks, per-input wheel or polynomial setup, and wrapper costs. Known factors are used only for validation outside the timed algorithm.</p>
{success_table}
<p class="note">All wheel variants completed 65 of the 75 original composites, while product search completed 68. The distribution matters: wheel Fermat wins on these balanced and close-factor sets, while product search succeeds on more unbalanced inputs. Timeouts are unresolved searches, never primality claims.</p></section>

<section><h2>Timing by factor geometry</h2>
<p>Every cell shows the median among completed inputs and its completion count. There are four balanced inputs per size and three close or unbalanced inputs per size. Do not treat a partial-success median as the timing for all inputs in that cell.</p>
{family_tables}
<p>For the 48-bit balanced set, wheel 20 alone has a 1.700 ms median, compared with 3.129 ms for plain parity Fermat and 14.264 ms for product search. With QR64 retained, wheel Fermat takes 0.758 ms, compared with 1.447 ms for the matched parity + QR64 control. At 56 bits the corresponding medians are 15.535, 37.983, 120.528, 8.608 and 20.754 ms.</p>
<p>At 64 bits both wheel variants complete all four balanced inputs; product search completes none within one second. The pure wheel median is 399.652 ms and wheel + QR64 is 161.658 ms. Native Brent rho remains much faster on this set, with a 0.781 ms median.</p>
<p>Across all original inputs, product search completes all 18 unbalanced semiprimes; each wheel variant completes nine. All 18 close-factor inputs finish at the first eligible Fermat midpoint, so the large advantage over product setup in that family comes from favorable Fermat geometry; wheel filtering is not the source of those first-test successes.</p></section>

<section><h2>Paired comparisons and attribution</h2>
<p>The following speedups divide comparator time by contender time for each matched successful input before aggregation. Values above one favor the contender. Unmatched completions are listed separately. The ≥1 ms rows are descriptive subsets used to reduce the influence of wrapper noise; they were not used for parameter selection.</p>
{pair_table}
<p>Native loop execution versus Python-orchestrated FLINT products is a comparison of these implementations. The matched native controls isolate the sieve itself. Differences of a few microseconds on easy inputs are noisy and should not be read as stable rankings.</p></section>

<section><h2>What the sieve proves about search cost</h2>
<p>Parity-correct Fermat visits half the integers. Wheel 20 visits either two or three midpoint classes per ten integers. Over a long interval, it removes 60% or 40% of the parity-eligible midpoint tests: candidate-count factors of 5/2 or 5/3, with bounded endpoint corrections. Those are not guaranteed wall-clock factors.</p>
<div class="math">R* = (p + q)/2; &nbsp; i* = (q − p)/2<br>
Δ = (p + q)/2 − √(pq) = (q − p)² / [2(√p + √q)²]<br>
Wheel midpoint tests = c<sub>N</sub> Δ + O(1), &nbsp; c<sub>N</sub> ∈ {{1/5, 3/10}}.</div>
<p>For close balanced factors with gap D = q − p, Δ is approximately D²/(8√N). In particular D = O(N<sup>1/4</sup>) gives only O(1) Fermat tests. The wheel makes the admissible portions of that region cheaper to traverse.</p>
<p>The usual uncapped, appropriately scaled polynomial-product construction has a balanced-input work scale of roughly N<sup>1/4</sup>, suppressing logarithmic factors. Comparing search scales places a close-factor crossover around D = N<sup>3/8</sup>, up to constants and logarithmic factors. This is not a measured timing threshold. For a fixed factor ratio q/p &gt; 1, Fermat requires Θ(√N) midpoint range; a fixed wheel does not change that exponent. The fixed 4096 degree cap means the current product baseline also does not realize the uncapped asymptotic indefinitely.</p></section>

<section><h2>Coupled geometry does not supply an absolute location</h2>
<p>The pair table preserves the coupled coordinate Y = 4Ri = q² − p². It yields Y ≡ 0 (mod 40) for N ending in 1 or 9, and Y ≡ 8 or 32 (mod 40) for N ending in 3 or 7. These are valid geometric restrictions. They follow from the complete pair table, so they are not independent additional tests after that table is enforced.</p>
<p>Here are two independently verified semiprimes with almost equal N, the same square-root starting point, the same complete wheel class, the same true pair lane, and the same height residue:</p>
{table(["Quantity", "Near factor midpoint", "Remote factor midpoint"], [
    ["N", "10,093,015,271", "10,093,010,311"],
    ["Factorization", "100459 × 100469", "19 × 531211069"],
    ["ceil(sqrt(N))", "100464", "100464"],
    ["N mod 20", "11", "11"],
    ["True (R mod 10, i mod 10)", "(4,5)", "(4,5)"],
    ["4Ri mod 40", "0", "0"],
    ["True R", "100464", "265605544"],
    ["Wheel midpoints through the factor", "1", "79,651,525"],
])}
<p>The remote case completes in 778.400 ms with wheel 20 alone, or 422.976 ms with wheel + QR64. Product search finds its small factor in 0.761 ms. For the near case, wheel 20 needs one test and about 0.016 ms, while product search takes 0.872 ms. The same residue information does not distinguish these two absolute locations.</p>
<p>In the divisor coordinate, every unit residue u modulo 20 remains compatible with N: choose v ≡ N/u (mod 20), then R = (u+v)/2 and i = (v−u)/2 give an admissible pair modulo ten. Thus the wheel supplies no further rejection of candidate-factor residues beyond the coprimality with 2 and 5 already known.</p>
<p class="good">The measured conclusion is positive and specific: this complete sieve makes Fermat traversal cheaper and beats this product baseline in favorable regions. Completeness of the periodic lanes does not itself tell us which absolute region contains a factor. Additional arithmetic information would be needed for that stronger claim.</p>
<details><summary>All supplemental measurements</summary>{supplement}</details></section>

<section><h2>Verification and reproduction</h2>
<p><strong>{html.escape(test_summary)}</strong></p>
<p>The 128 new tests include all ten odd modulo-20 tables, both parity lanes in digit-only unions, the 551 regression, the coupled height residues, 20,000 small-input searches across four modes, independent first-solution and exact candidate-count checks, native intermediates exceeding 64 bits, and time/work budgets. All {len(RAW):,} timed calls produced valid outcomes, with no invalid factors, false prime claims, algorithm errors, or mixed repeat outcomes.</p>
<pre>python3 -m pip install -r requirements.txt
python3 build_native.py
python3 build_wheel20.py
python3 wheel20_factor.py 10093015271 --qr64 --timeout 10
python3 -m pytest -q
python3 wheel20_benchmark.py --outdir results/wheel20-reproduction
python3 wheel20_report.py</pre>
<p>The report generator reads the bundled original <code>results/wheel20</code> measurements. The reproduction command writes a separate directory. Raw runs, exact inputs, configurations, environment data, source hashes, and this report generator are included. The native wheel implementation accepts integers below 2<sup>64</sup>; unresolved searches return <code>status: budget</code>.</p>
<details><summary>Environment and measurement protocol</summary><pre>{html.escape(json.dumps(ENV,indent=2))}</pre><pre>{html.escape(json.dumps(PROTOCOL,indent=2))}</pre></details></section>
</main></body></html>'''
(OUT/"wheel20_report.html").write_text(document)
print(OUT/"wheel20_report.html")

