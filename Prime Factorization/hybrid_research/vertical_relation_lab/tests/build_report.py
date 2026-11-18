import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]/'results'
base=json.loads((R/'relation_benchmark.json').read_text())
exit_=json.loads((R/'factor_exit_benchmark.json').read_text())
new=json.loads((R/'newton_window_benchmark.json').read_text())
bounds=json.loads((R/'smoothness_bound_sensitivity.json').read_text())
by={(a['case'],a['method']):a for a in base['results']}
byexit={(a['case'],a['method']):a for a in exit_['results']}
methods=('QS','CTM','CTMWalk','Grouping','GroupingBand','SquareY','NearPlus','NearMinus')
s=[]
def add(t=''):s.append(t)
add('# Squared complex-plane Newton lattice and relation-generation experiment')
add('')
add('**Scope:** A mathematically exact Fermat/CRT factor-point region sieve and GMP Newton-lattice certificate, followed by a controlled smooth-relation comparison. **Not an implemented GNFS**, and not a GNFS-scale performance claim. Original source: `Complex space.md`; mathematical proofs in `math/FORMALIZATION.md`.')
add('')
add('## Verified insight: the vertical green line')
add('')
add(r'For an exact factor point, $z=a+bi$, $a^2-b^2=N$, $Y=2ab$, $H=a^2+b^2$, so $H^2=N^2+Y^2$, $(H-Y)(H+Y)=N^2$, and $H-Y=p^2$, $H+Y=q^2$. The curve $H(Y)=\sqrt{N^2+Y^2}$ is increasing and convex; its **rounded integer residual is not monotone**. For any off-curve CTM point with $E=a^2-b^2-N$, the squared residual is exactly $E(2N+E)$; a smooth $E$ does not create a new smooth square relation by itself.')
add('')
add('The factor-point wheel uses unknown factor residues $p$ modulo $2m$ only as **a finite modular enumeration**, not as factor knowledge. For $m=420$ and $N=978508015703$, only six paired $(Y,H)$ residues modulo 420 are possible. The exact rational Newton bounds and Euclidean floor sum are proved in the mathematical appendix. The wheel **must not** be applied as an exclusion sieve to general nonzero-residual relations.')
add('')
add('## Exact regional geometry, 80 deterministic windows per case')
add('')
add('| Input | Bits | Empty / 80 with integer ceilings | Empty / 80 after 1 Newton step | Surviving integer envelope points after Newton |')
add('|---|---:|---:|---:|---:|')
for name in ('40b','64b','128b'):
 x=next(z for z in new['results'] if z['case']==name and z['newton_steps']==0)
 y=next(z for z in new['results'] if z['case']==name and z['newton_steps']==1)
 add(f"| {name} | {x['N_bits']} | {x['empty_tiles']} | {y['empty_tiles']} | {y['candidate_envelope_points']} |")
add('')
add('Windows were 100,000 $Y$ units wide, distributed by an $N$-only fixed schedule on $0\\le Y\\le 3N/4$ (valid as a full factor-domain bound only when $p<q<2p$). The observed *factor-bearing* 128-bit validation window retains its exact integer point: 477 candidate lattice points under integer ceilings, **1 after Newton**. Local exclusion is not a full-domain speedup: there can be an immense number of separate windows.')
add('')
add('## Fixed-work relation generation: B = 1021, K = 32768')
add('')
add('All modes use the same GMP prime-base smoothness test, the required one- or two-sided exponent matrix, full sign tracking, duplicate-residual suppression, exact mod-$N$ congruences and gcd verification. The factor inputs were withheld from the generator binaries. `CTM` is the nearest-boundary Fermat projection; `CTMWalk` is a 4-lane increment-only CTM staircase. `GroupingBand` uses only the specified balance constraint. `SquareY` samples the vertical line uniformly in the balanced range.')
add('')
add('| Actual bits | QS valid / deps | CTM valid / deps | CTMWalk valid / deps | Grouping-band valid / deps | Vertical Y valid |')
add('|---:|---:|---:|---:|---:|---:|')
for c in base['cases']:
 name=c['case'];bits=c['N_bits']
 def ent(m):
  r=by[(name,m)];return f"{r['valid_relations']} / {r['dependencies']}"
 add(f"| {bits} | {ent('QS')} | {ent('CTM')} | {ent('CTMWalk')} | {ent('GroupingBand')} | {by[(name,'SquareY')]['valid_relations']} |")
add('')
add('In the 40-bit manuscript example, CTM produced **11,041 distinct smooth residual hits**, but only **378 valid two-sided smooth relations** (rank 320). `CTMWalk` produced 311 valid relations, grouping-band 327, and QS 370. Unlike the others, QS starts with a square on one side; counting raw smooth residuals as QS-quality relations would be misleading.')
add('')
add('## Time to first independently verified factor under the same maximum candidate budget')
add('')
add('Seven of eight semiprimes were resolved by the simple QS relation generator, four by each CTM strategy, three by the balanced grouping-band approach, none by direct vertical-Y sampling. `—` means no proper factor within 32,768 candidates. Timings are three-run medians of the *C++ computation* including candidate construction, smooth factoring, matrix/dependency work and gcd verification, excluding process start.')
add('')
add('| Actual bits | QS | CTM projection | CTM walk | Grouping band |')
add('|---:|---:|---:|---:|---:|')
for c in base['cases']:
 name=c['case'];bits=c['N_bits']
 def tm(m):
  r=byexit[(name,m)];return f"{r['median_kernel_time_ms']:.2f} ms" if r['factor_found'] else '—'
 add(f"| {bits} | {tm('QS')} | {tm('CTM')} | {tm('CTMWalk')} | {tm('GroupingBand')} |")
add('')
add('## Smoothness sensitivity')
add('')
add('At 40 bits, with the same 32,768 candidates:')
add('')
add('| Smoothness B | QS valid / deps | CTM valid / deps | CTM walk valid / deps | Grouping band valid / deps |')
add('|---:|---:|---:|---:|---:|')
for B in (257,1021,4093):
 def fmt(m):
  r=next(z for z in bounds['results'] if z['case']=='40b_manuscript' and z['B']==B and z['method']==m)
  return f"{r['valid']} / {r['dependencies']}"
 add(f"| {B} | {fmt('QS')} | {fmt('CTM')} | {fmt('CTMWalk')} | {fmt('GroupingBand')} |")
add('')
add('CTM and grouping can produce genuine factors by combining smooth relations, especially at larger B. But they require both sides smooth, and at B=4093 CTM projection needed 345 ms for its full-work 40-bit sample compared with 187 ms for QS, with more columns in the linear algebra. Results do not establish a consistent efficiency gain.')
add('')
add('## Valid relations versus smooth-looking residuals')
add('')
add('1. **CTM:** $L=(a-b)(a+b)=N+E$, so $L\\equiv E\\pmod N$. Two-sided smoothness of $L$ and $E$ supports a proper parity matrix and a possible square congruence. Different CTM paths with the same $E$ give the **same** arithmetic relation and are deduplicated.')
add('2. **Grouping:** $N=g\\lfloor N/g\\rfloor+r$, with $E=-r$ and $L=N+E$. Its relations collapse to the same smooth-neighbor form; multiplicity of successful seed chains does not automatically imply the same multiplicity of independent congruences.')
add('3. **Squared vertical line:** $H^2-Y^2=N^2+W$. A generic smooth $W$ requires the left-hand product to be smooth for our two-sided relation protocol. For a CTM-derived point, $W=E(2N+E)$ and $H^2-Y^2=(N+E)^2$; modulo $N$ this is a tautological square relation, not independent data.')
add('4. **QS baseline:** $x^2-N=f(x)$; only $f(x)$ needs to be smooth, and a dependency of its exponent vectors produces a genuine square congruence.')
add('')
add('## Verdict and limitations')
add('')
add('The rational-Newton vertical sieve provides a strong, rigorously certified **local** geometry filter. The separate relation experiment confirms that CTM and grouping can produce valid dependencies and actual factors, but this implementation does not outperform the conventional QS relation generator on the tested hard end, and it gives no GNFS-level scaling. The earlier uploaded `gnfs_gmp.cpp` cannot be a GNFS control until its algebraic norm and other NFS machinery are corrected; there was no CADO-NFS installation in the environment. The tests cover 28–63 bits with fixed K and factor-base bounds, not RSA-100/110 or GNFS-relevant input sizes.')
add('')
add('The remaining mathematical opportunity is a new **one-sided square relation** from the complex/grouping coordinates, or a truly cheaper way to identify useful $Y$ intervals without visiting a square-root number of them. Those would change the relevant scaling; more smooth *nonzero* residuals alone do not.')
add('')
add('## Files and reproducibility')
add('')
add('`math/FORMALIZATION.md` is the proof reference. `src/vertical_lattice.cpp` and `src/relation_bench.cpp` are the standalone C++17/GMP implementations. `tests/test_vertical.py`, `tests/test_newton.py`, and `tests/test_relations.py` provide independent mathematical checks. Raw JSON with every test case, full- and early-stop metrics, smoothness sensitivity, and environment notes is in `results/`. All commands are in `README.md`.')
(R.parent/'REPORT.md').write_text('\n'.join(s)+'\n')
print('report lines',len(s))
