"""Reproducible C++ GMP versus Python/FLINT validation, with startup excluded.

The executable is kept alive in --batch mode. Primary timings are each
implementation's internally measured search time; external per-request wall
times are retained separately. Both methods include per-input arithmetic setup.
The benchmark passes only N and fixed controls to algorithms. Reference factors
are used for validation after each timed call.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import select
import statistics
import subprocess
import sys
from time import perf_counter

from sympy import isprime

from vf_factor import fusion_native
from vf_factor.corpus import (build_corpus, validate_corpus, _balanced_pair,
                              _close_pair, _random_pair)
from vf_factor.fused import factor


ROOT = Path(__file__).resolve().parent
FRESH_SEED = 926491
LARGE_SEED = 926493
ORDER_SEED = 926497
CONFIG = {"initial_degree": 16, "max_degree": 512, "fermat_ratio": 2.0,
          "initial_fermat_tests": 64}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def large_corpus(seed=LARGE_SEED):
    rng = random.Random(seed)
    cases = []
    for bits in (80, 96, 128, 192, 256):
        for kind, index, pair in (
            ("close", 0, _close_pair(rng, bits)),
            ("close", 1, _close_pair(rng, bits)),
            ("small_factor", 16, _random_pair(rng, bits, 16)),
            ("small_factor", 20, _random_pair(rng, bits, 20)),
            ("balanced", 0, _balanced_pair(rng, bits)),
        ):
            p, q = pair
            cases.append(dict(case_id=f"gmp-{kind}-{bits}-{index}", family=kind,
                              n=p*q, p=p, q=q, bits=(p*q).bit_length(),
                              target_bits=bits, expected_status="factor",
                              seed=seed, split="large", small_factor_bits=p.bit_length(),
                              provenance="fixed-seed prime generation before algorithm calls"))
    validate_corpus(cases)
    return cases


class CppWorker:
    """A persistent process; each input line starts a new independent budget."""
    def __init__(self, binary, timeout, extra_options=()):
        args = [str(binary), "--batch", "--timeout", str(timeout),
                "--initial-degree", str(CONFIG["initial_degree"]),
                "--max-degree", str(CONFIG["max_degree"]),
                "--fermat-ratio", str(CONFIG["fermat_ratio"]),
                "--initial-fermat-tests", str(CONFIG["initial_fermat_tests"]),
                *extra_options]
        self.process = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, text=True, bufsize=1)
        self.response_timeout = max(10.0, 4.0*timeout)
        # Complete one cheap warm-up before the first measured request.
        result, _ = self.call(77)
        assert result["status"] == "factor"

    def call(self, n):
        started = perf_counter()
        self.process.stdin.write(str(n) + "\n")
        self.process.stdin.flush()
        ready, _, _ = select.select([self.process.stdout], [], [], self.response_timeout)
        if not ready:
            self.process.kill()
            raise TimeoutError("GMP batch request exceeded the external safety timeout")
        line = self.process.stdout.readline()
        wall = perf_counter() - started
        if not line:
            raise RuntimeError(f"GMP process stopped: {self.process.stderr.read()}")
        result = json.loads(line)
        if result.get("status") == "error":
            raise RuntimeError(result)
        return result, wall

    def close(self):
        self.process.stdin.close()
        status = self.process.wait(timeout=10)
        error = self.process.stderr.read()
        if status != 0:
            raise RuntimeError(f"GMP process exited {status}: {error}")


def validate_result(case, result):
    n = case["n"]
    assert int(result["n"]) == n, (case, result)
    status = result["status"]
    if status == "factor":
        d = int(result["factor"])
        assert 1 < d < n and n % d == 0, (case, result)
        assert int(result["cofactor"]) == n//d
        assert case["expected_status"] == "factor"
        return True
    if status in ("prime", "probable_prime"):
        assert case["expected_status"] == "prime" and isprime(n), (case, result)
        assert result["factor"] is None
        return status == "prime"  # A probable result is not a prime certificate.
    assert status == "budget", (case, result)
    assert result["factor"] is None
    return False


def pin_and_environment(binary):
    before = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    if before and hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {before[0]})
    import flint
    flint.ctx.threads = 1
    metadata = json.loads(subprocess.check_output([str(binary), "--metadata"], text=True))
    return dict(python=sys.version, platform=platform.platform(), machine=platform.machine(),
                original_affinity=before,
                active_affinity=sorted(os.sched_getaffinity(0)) if before else None,
                python_flint_version=flint.__version__, flint_threads=flint.ctx.threads,
                python_fermat_backend=fusion_native.backend_info(),
                python_native_kernel_sha256=digest(ROOT/"vf_factor/fusion_kernel.so"),
                cpp_metadata=metadata, binary_sha256=digest(binary))


def run_dataset(cases, name, outdir, binary, timeout, repeats, methods, rng):
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir/"benchmark-corpus.json").write_text(json.dumps(cases, indent=2)+"\n")
    ordered = list(cases)
    rng.shuffle(ordered)
    workers = {}
    if "gmp_cpp" in methods:
        workers["gmp_cpp"] = CppWorker(binary, timeout)
    if "gmp_cpp_all_mpz" in methods:
        workers["gmp_cpp_all_mpz"] = CppWorker(binary, timeout, ("--gmp-fermat",))
    # Python imports, dispatch, and shared-library loading are warmed here.
    if "python_fused" in methods:
        factor(77, timeout=timeout, **CONFIG)
        factor(978508015703, timeout=timeout, **CONFIG)
    raw_rows, summaries = [], []
    path = outdir/"benchmark-runs.jsonl"
    try:
        with path.open("w") as raw_file:
            for position, case in enumerate(ordered):
                order = list(methods)
                rng.shuffle(order)
                for method in order:
                    attempts = []
                    for repeat in range(repeats):
                        if method in workers:
                            result, wall = workers[method].call(case["n"])
                        else:
                            tick = perf_counter()
                            result = factor(case["n"], timeout=timeout, **CONFIG).to_dict()
                            wall = perf_counter()-tick
                        solved = validate_result(case, result)
                        row = dict(dataset=name, case_id=case["case_id"], family=case["family"],
                                   bits=case["bits"], target_bits=case["target_bits"], n=case["n"],
                                   algorithm=method, repeat=repeat, timeout_seconds=timeout,
                                   status=result["status"], factor=int(result["factor"]) if result["factor"] else None,
                                   solved=solved, seconds=result["elapsed_seconds"],
                                   wall_seconds=wall, counters=result.get("counters", {}),
                                   details=result.get("details", {}), result=result)
                        assert row["seconds"] >= 0 and wall >= 0
                        raw_file.write(json.dumps(row)+"\n")
                        raw_file.flush()
                        raw_rows.append(row)
                        attempts.append(row)
                        if not solved:
                            break
                    med = statistics.median(r["seconds"] for r in attempts)
                    representative = min(attempts, key=lambda r: abs(r["seconds"]-med))
                    summary = {k: representative[k] for k in
                               ("dataset", "case_id", "family", "bits", "target_bits", "n",
                                "algorithm", "status", "factor", "counters", "details")}
                    summary.update(solved=all(r["solved"] for r in attempts) and len(attempts)==repeats,
                                   attempts=len(attempts), requested_repeats=repeats,
                                   median_seconds=med,
                                   median_wall_seconds=statistics.median(r["wall_seconds"] for r in attempts),
                                   min_seconds=min(r["seconds"] for r in attempts),
                                   max_seconds=max(r["seconds"] for r in attempts),
                                   mixed_status=len({r["status"] for r in attempts}) > 1)
                    summaries.append(summary)
                if (position+1) % 10 == 0 or position+1 == len(ordered):
                    print(json.dumps(dict(dataset=name, completed_inputs=position+1,
                                          total_inputs=len(ordered))), flush=True)
    finally:
        for worker in workers.values():
            worker.close()
    summaries.sort(key=lambda r: (r["case_id"], r["algorithm"]))
    (outdir/"benchmark-summary.json").write_text(json.dumps(summaries, indent=2)+"\n")
    fields = [k for k in summaries[0] if k not in ("counters", "details")]
    with (outdir/"benchmark-summary.csv").open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k:r[k] for k in fields} for r in summaries)
    stats = {}
    for method in methods:
        rows = [r for r in summaries if r["algorithm"]==method]
        composites = [r for r in rows if r["family"]!="prime"]
        stats[method] = dict(composites_completed=sum(r["solved"] for r in composites),
                             composites=len(composites),
                             primes_completed=sum(r["solved"] for r in rows if r["family"]=="prime"),
                             primes=sum(r["family"]=="prime" for r in rows),
                             mixed_repeat_cases=sum(r["mixed_status"] for r in rows))
    analysis = dict(dataset=name, inputs=len(cases), summary_rows=len(summaries),
                    timed_calls=len(raw_rows), outcomes=stats)
    (outdir/"analysis.json").write_text(json.dumps(analysis, indent=2)+"\n")
    return analysis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT/"build/vf_factor_gmp")
    parser.add_argument("--outdir", type=Path, default=ROOT/"results/gmp")
    parser.add_argument("--timeout", type=float, default=0.5)
    parser.add_argument("--large-timeout", type=float, default=0.5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--fresh-only", action="store_true")
    parser.add_argument("--all-mpz-control", action="store_true",
                        help="also measure the GMP Fermat kernel below64 bits")
    args = parser.parse_args()
    if args.timeout <= 0 or args.large_timeout <= 0 or args.repeats < 1:
        parser.error("timeouts and repeats must be positive")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error("build the GMP executable before benchmarking")
    if not fusion_native.is_available():
        parser.error("run python3 build_fusion.py so the Python comparison has the same native64 Fermat kernel")
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    fresh = build_corpus(seed=FRESH_SEED, split="heldout")
    validate_corpus(fresh)
    large = large_corpus()
    replay = json.loads((ROOT/"results/main/benchmark-corpus.json").read_text())
    validate_corpus(replay)
    environment = pin_and_environment(binary)
    source_paths = sorted(list((ROOT/"cpp").glob("*.hpp")) + list((ROOT/"cpp").glob("*.cpp")))
    source_paths += [ROOT/p for p in (
        "gmp_benchmark.py", "build_gmp.py", "Makefile.gmp",
        "vf_factor/fused.py", "vf_factor/shared_products.py", "vf_factor/poly_engine.py",
        "vf_factor/fusion_native.py", "vf_factor/fusion_kernel.cpp", "vf_factor/common.py",
        "vf_factor/vf_algorithm.py", "vf_factor/baselines.py", "vf_factor/corpus.py",
    ) if (ROOT/p).exists()]
    protocol = dict(fresh_seed=FRESH_SEED, large_seed=LARGE_SEED, order_seed=ORDER_SEED,
                    frozen_options=CONFIG, timeout_seconds=args.timeout,
                    large_timeout_seconds=args.large_timeout, repeats=args.repeats,
                    methods=["gmp_cpp", "python_fused"] + (["gmp_cpp_all_mpz"] if args.all_mpz_control else []),
                    timing="Primary: internal search elapsed_seconds. External request wall times retained. Process/import startup and one-time warm-up excluded; per-input setup included.",
                    failure_policy="Stop repetitions after a failure. A case succeeds only when every requested repetition completed.",
                    independence="New generated semiprimes use seeds926491/926493. Published examples, prime cubes and anchor regressions remain fixed controls. No per-input factors enter algorithm calls.",
                    scheduling="Both algorithms use time-dependent Fermat allocation; work trajectories can differ between implementations and repetitions.",
                    probable_primes="The GMP default can return probable_prime above64 bits. Those are never counted as proved primes. Large comparison inputs here are known composites.",
                    source_sha256={str(p.relative_to(ROOT)):digest(p) for p in source_paths},
                    corpus_sha256={"fresh": hashlib.sha256(json.dumps(fresh, sort_keys=True).encode()).hexdigest(),
                                   "large": hashlib.sha256(json.dumps(large, sort_keys=True).encode()).hexdigest(),
                                   "original_replay": hashlib.sha256(json.dumps(replay, sort_keys=True).encode()).hexdigest()})
    (out/"protocol.json").write_text(json.dumps(protocol, indent=2)+"\n")
    (out/"environment.json").write_text(json.dumps(environment, indent=2)+"\n")
    rng = random.Random(ORDER_SEED)
    analyses = [run_dataset(fresh, "fresh", out/"fresh", binary, args.timeout,
                            args.repeats, protocol["methods"], rng)]
    if not args.fresh_only:
        analyses.append(run_dataset(replay, "original_replay", out/"original-replay", binary,
                                    args.timeout, args.repeats, ["gmp_cpp"], rng))
        analyses.append(run_dataset(large, "large", out/"large", binary, args.large_timeout,
                                    args.repeats, ["gmp_cpp", "python_fused"], rng))
    (out/"analysis.json").write_text(json.dumps(analyses, indent=2)+"\n")
    print(json.dumps(analyses, indent=2), flush=True)


if __name__ == "__main__":
    main()
