"""Small independent-seed pilot for two scheduling choices, never factor hints."""
from __future__ import annotations
import json
import statistics
from pathlib import Path

from vf_factor.benchmark import run
from vf_factor.corpus import build_corpus
from vf_factor.fused import factor


def main():
    out=Path("results/fused-tuning")
    corpus=build_corpus(seed=826391,split="heldout")
    cases=[]
    for row in corpus:
        if row["target_bits"] < 40:
            continue
        if row["family"]=="balanced" and row["case_id"].endswith(("-0","-1")):
            cases.append(row)
        elif row["family"]=="close" and row["case_id"].endswith("-0"):
            cases.append(row)
        elif row["family"]=="unbalanced" and row["small_factor_bits"] in (8,16 if row["target_bits"]==40 else 20):
            cases.append(row)
    assert len(cases)==20
    options={f"ratio{ratio}_degree{degree}":dict(fermat_ratio=ratio,max_degree=degree)
             for ratio in (1,2) for degree in (512,1024)}
    registry={name:factor for name in options}
    rows=run(cases,list(registry),timeout=.25,repeats=2,seed=826392,
             outdir=str(out),options=options,algorithm_registry=registry)
    assert not any(r["status"]=="error" for r in rows)
    # Selection rule fixed in advance: maximize completions, then minimize
    # geometric mean time with timeout outcomes assigned the full budget.
    scored=[]
    for name in options:
        rs=[r for r in rows if r["algorithm"]==name]
        score=statistics.geometric_mean(r["median_seconds"] if r["solved"] else .25 for r in rs)
        scored.append(dict(name=name,completed=sum(r["solved"] for r in rs),
                           geometric_seconds=score,options=options[name]))
    scored.sort(key=lambda r:(-r["completed"],r["geometric_seconds"]))
    (out/"selection.json").write_text(json.dumps(dict(rule="max completions, then minimum geometric time; budget failures=0.25s",
                                                     seed=826391,selected=scored[0],ranking=scored),indent=2))
    print(json.dumps(scored,indent=2))


if __name__=="__main__":
    main()
