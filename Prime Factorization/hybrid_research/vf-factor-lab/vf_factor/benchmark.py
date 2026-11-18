"""Reproducible single-process, single-thread benchmark with censored failures.

Known factors are used only by validation after a timed search. Algorithm
functions receive N, timeout and fixed public options, never corpus factors.
Completed cases are repeated; budget-limited cases are not rerun needlessly.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import platform
import random
import resource
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import flint

from .corpus import build_corpus
from .registry import get_algorithms


def _json(value):
    return json.dumps(value, sort_keys=True, default=str)


def validate_result(case, result):
    if result.factor is not None:
        assert 1 < result.factor < case["n"] and case["n"] % result.factor == 0
        return True
    if result.status == "prime":
        from sympy import isprime
        assert isprime(case["n"]), f"false prime for {case['n']}"
        return True
    assert result.status in ("budget", "exhausted"), result
    return False


def run(cases, names, *, timeout=1.0, repeats=3, seed=20261117, outdir="results", options=None,
        algorithm_registry=None):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    methods = get_algorithms() if algorithm_registry is None else algorithm_registry
    options = options or {}
    flint.ctx.threads = 1
    affinity = None
    if hasattr(os, "sched_getaffinity"):
        affinity = max(os.sched_getaffinity(0))
        os.sched_setaffinity(0, {affinity})
    corpus_text = _json(cases)
    (out / "benchmark-corpus.json").write_text(json.dumps(cases, indent=2))
    metadata = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version, platform=platform.platform(), cpu_count=os.cpu_count(),
        flint_python=flint.__version__, flint_library=flint.__FLINT_VERSION__,
        flint_threads=flint.ctx.threads, timeout_seconds=timeout,
        pinned_logical_cpu=affinity,
        completed_repeats=repeats, algorithm_order_seed=seed,
        corpus_sha256=hashlib.sha256(corpus_text.encode()).hexdigest(),
        algorithms=names, options=options,
        timing="perf_counter around entire algorithm call, imports prewarmed; includes prechecks and polynomial setup",
        censoring="budget/exhausted rows are failures, never assigned fast successful times; no additional repeats after a budget result",
        fairness="Python-loop baselines, optional C++ baselines, and compiled FLINT reference identified separately")
    (out / "environment.json").write_text(json.dumps(metadata, indent=2))
    # Warm imports and native dispatch, not per-input polynomial state.
    for name in names:
        methods[name](8051, timeout=.2, **options.get(name, {}))
    rng = random.Random(seed)
    raw, summary = [], []
    with (out / "benchmark-runs.jsonl").open("w") as log:
        for case_index, case in enumerate(cases):
            order = list(names)
            rng.shuffle(order)
            for name in order:
                attempts = []
                for repeat in range(repeats):
                    started = perf_counter()
                    try:
                        result = methods[name](case["n"], timeout=timeout, **options.get(name, {}))
                        elapsed = perf_counter() - started
                        solved = validate_result(case, result)
                        row = dict(case_id=case["case_id"], family=case["family"],
                                   bits=case["bits"], target_bits=case["target_bits"],
                                   n=case["n"], algorithm=name, repeat=repeat,
                                   wall_seconds=elapsed, solved=solved,
                                   **{k:v for k,v in result.to_dict().items() if k not in ("n","algorithm")})
                    except Exception as exc:
                        row = dict(case_id=case["case_id"], family=case["family"],
                                   bits=case["bits"], target_bits=case["target_bits"], n=case["n"],
                                   algorithm=name, repeat=repeat, wall_seconds=perf_counter()-started,
                                   solved=False, status="error", factor=None,
                                   error=f"{type(exc).__name__}: {exc}", counters={}, details={})
                    log.write(_json(row) + "\n")
                    log.flush()
                    raw.append(row)
                    attempts.append(row)
                    if not row["solved"]:
                        break
                times = [r["wall_seconds"] for r in attempts]
                successful = all(r["solved"] for r in attempts)
                representative = min(attempts, key=lambda r: abs(r["wall_seconds"]-statistics.median(times)))
                aggregate = dict(case_id=case["case_id"], family=case["family"],
                    bits=case["bits"], target_bits=case["target_bits"], n=case["n"],
                    algorithm=name, solved=successful, status=representative["status"],
                    median_seconds=statistics.median(times), min_seconds=min(times),
                    max_seconds=max(times), attempts=len(attempts), factor=representative["factor"],
                    counters=representative["counters"], details=representative["details"])
                if "error" in representative:
                    aggregate["error"] = representative["error"]
                summary.append(aggregate)
            print(f"{case_index+1}/{len(cases)} {case['case_id']} ({case['bits']} bits)", flush=True)
    (out / "benchmark-summary.json").write_text(json.dumps(summary, indent=2))
    metadata["process_peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    metadata["completed_timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    (out / "environment.json").write_text(json.dumps(metadata, indent=2))
    columns = ["case_id","family","bits","target_bits","n","algorithm","solved","status",
               "median_seconds","min_seconds","max_seconds","attempts","factor","counters","details"]
    with (out / "benchmark-summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in summary:
            writer.writerow({**row, "counters":_json(row["counters"]), "details":_json(row["details"])})
    return summary


def main():
    methods = get_algorithms()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=1.0)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--split", choices=["heldout","tuning"], default="heldout")
    parser.add_argument("--methods", nargs="+", choices=methods, default=list(methods))
    parser.add_argument("--outdir", default="results/main")
    parser.add_argument("--corpus", help="read a fixed JSON corpus instead of generating one")
    args = parser.parse_args()
    cases = json.loads(Path(args.corpus).read_text()) if args.corpus else build_corpus(split=args.split)
    run(cases, args.methods, timeout=args.timeout, repeats=args.repeats, outdir=args.outdir)


if __name__ == "__main__":
    main()
