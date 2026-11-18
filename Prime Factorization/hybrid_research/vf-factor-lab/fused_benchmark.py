"""Frozen unified-algorithm validation on a fresh seed and original-corpus replay."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from vf_factor import fusion_native
from vf_factor.benchmark import run
from vf_factor.corpus import build_corpus, validate_corpus
from vf_factor.fused import factor
from vf_factor.registry import get_algorithms

ROOT=Path(__file__).resolve().parent
FRESH_SEED=926483
CONFIG={"max_degree":512,"fermat_ratio":2}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout",type=float,default=1.0)
    parser.add_argument("--repeats",type=int,default=3)
    parser.add_argument("--outdir",default=str(ROOT/"results/fused"))
    parser.add_argument("--fresh-only",action="store_true")
    args=parser.parse_args()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    registry=get_algorithms()
    names=["fused","fused_no_vf","fused_no_feedback","fused_unshared","adaptive_product",
           "fermat160","vf_fast","product_trial","native_pollard_brent","flint_native"]
    for name in names[:6]:
        registry[name]=factor
    options={
        "fused":dict(CONFIG),
        "fused_no_vf":{**CONFIG,"vf_probes":False},
        "fused_no_feedback":{**CONFIG,"feedback":False},
        "fused_unshared":{**CONFIG,"shared_polynomial":False},
        "adaptive_product":{**CONFIG,"fermat_enabled":False,"vf_probes":False},
        "fermat160":{**CONFIG,"vf_probes":False,"initial_fermat_tests":(1<<63)-1},
    }
    sources=("vf_factor/fused.py","vf_factor/shared_products.py","vf_factor/fusion_native.py",
             "vf_factor/fusion_kernel.cpp","vf_factor/poly_engine.py","vf_factor/benchmark.py",
             "fused_benchmark.py")
    protocol={"fresh_seed":FRESH_SEED,"tuning_seed":826391,
              "selection":json.loads((ROOT/"results/fused-tuning/selection.json").read_text()),
              "frozen_options":options,"timeout_seconds":args.timeout,"repeats":args.repeats,
              "native":fusion_native.backend_info(),
              "synthesis":"Reciprocal exact factor bounds; shared unit-step polynomial for divisor and VF blocks; one aggregate GCD; persistent wheel-160 Fermat state; one deadline",
              "scheduling":"Fermat gets a bounded time slice based on the preceding polynomial batch cost; trajectories may vary across timing repeats",
              "source_sha256":{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              "independence":"Final fresh generated cases use seed 926483, distinct from tuning seed 826391 and integer-root development seed 926477. Published examples, anchor regressions and prime cubes are fixed controls, not new random holdouts."}
    (out/"protocol.json").write_text(json.dumps(protocol,indent=2))
    fresh=build_corpus(seed=FRESH_SEED,split="heldout")
    validate_corpus(fresh)
    tuning_values={c["n"] for c in json.loads((ROOT/"results/fused-tuning/benchmark-corpus.json").read_text())}
    assert not any(c["n"] in tuning_values for c in fresh if c["family"] in {"balanced","close","unbalanced"})
    rows=run(fresh,names,timeout=args.timeout,repeats=args.repeats,seed=926478,
             outdir=str(out/"fresh"),options=options,algorithm_registry=registry)
    assert not any(r["status"]=="error" for r in rows)
    if not args.fresh_only:
        previous=json.loads((ROOT/"results/main/benchmark-corpus.json").read_text())
        rows=run(previous,["fused"],timeout=args.timeout,repeats=args.repeats,seed=926479,
                 outdir=str(out/"original-replay"),options=options,algorithm_registry=registry)
        assert not any(r["status"]=="error" for r in rows)
        stress=json.loads((ROOT/"results/stress-corpus.json").read_text())
        registry["pollard_brent"]=get_algorithms()["pollard_brent"]
        rows=run(stress,["fused","adaptive_product","pollard_brent"],timeout=2.0,repeats=1,seed=926480,
                 outdir=str(out/"stress"),options=options,algorithm_registry=registry)
        assert not any(r["status"]=="error" for r in rows)


if __name__=="__main__":
    main()
