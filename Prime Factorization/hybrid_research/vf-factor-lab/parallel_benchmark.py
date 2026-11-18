"""Frozen regional-search scaling experiment: one, two, and four workers.

The same parallel search policy is measured at each worker count. The older
sequential GMP factorizer is a separate algorithm control. Per-call internal
wall time, process CPU time, and actual represented work are retained. Success
requires all requested repetitions; timeouts are not reported as speedups.
"""
from __future__ import annotations

import argparse
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

from vf_factor.corpus import _balanced_pair, _close_pair, _random_pair, validate_corpus


ROOT = Path(__file__).resolve().parent
CORPUS_SEED = 926509
ORDER_SEED = 926513
REGIONAL_CONFIG = {"degree": 128, "fermat_span": 655360,
                   "job_seconds": 0.025, "initial_fermat_tests": 64,
                   "fermat_ratio": 2.0}
METHODS = {"regional_t1": 1, "regional_t2": 2, "regional_t4": 4,
           "sequential_gmp": 1}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def corpus(seed=CORPUS_SEED):
    rng = random.Random(seed)
    cases = []

    def add(family, bits, index, pair, **extra):
        p, q = pair
        cases.append(dict(case_id=f"regional-{family}-{bits}-{index}", family=family,
                          n=p*q, p=p, q=q, bits=(p*q).bit_length(), target_bits=bits,
                          expected_status="factor", seed=seed,
                          provenance="fixed-seed prime generation before any timed search",
                          **extra))

    for bits in (48, 56, 60, 64):
        for index in range(3):
            add("balanced", bits, index, _balanced_pair(rng, bits))
    for bits in (64, 128):
        for pbits in (20, 24):
            for index in range(2):
                add("small_factor", bits, f"p{pbits}-{index}",
                    _random_pair(rng, bits, pbits), small_factor_bits=pbits)
    for bits in (64, 128, 192, 256):
        p, q = _close_pair(rng, bits)
        add("close_control", bits, 0, (p, q), factor_gap=q-p)
    for index in range(2):
        add("hard_balanced", 80, index, _balanced_pair(rng, 80))
    validate_corpus(cases)
    assert len(cases) == 26 and len({c["n"] for c in cases}) == 26
    return cases


def cpu_snapshot(pid):
    """Fallback process CPU measurement for the unchanged sequential binary."""
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
        fields = text[text.rfind(")")+2:].split()
        # After comm, index0 is stat field3; fields14/15 are user/system ticks.
        return (int(fields[11])+int(fields[12]))/os.sysconf("SC_CLK_TCK")
    except (OSError, ValueError, IndexError):
        return None


class Worker:
    def __init__(self, binary, timeout, affinity, options=()):
        command = [str(binary), "--batch", "--timeout", str(timeout), *options]
        before_exec = None
        if affinity and hasattr(os, "sched_setaffinity"):
            # This benchmark driver is single-threaded when workers are spawned.
            before_exec = lambda: os.sched_setaffinity(0, set(affinity))
        self.process = subprocess.Popen(command, stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        text=True, bufsize=1, preexec_fn=before_exec)
        self.external_limit = max(15, timeout*5)
        self.command = command
        self.affinity = sorted(os.sched_getaffinity(self.process.pid)) if hasattr(os, "sched_getaffinity") else None
        result, *_ = self.call(77)
        assert result["status"] == "factor", result

    def call(self, n):
        cpu_before = cpu_snapshot(self.process.pid)
        began = perf_counter()
        self.process.stdin.write(str(n)+"\n")
        self.process.stdin.flush()
        ready, _, _ = select.select([self.process.stdout], [], [], self.external_limit)
        if not ready:
            self.process.kill()
            raise TimeoutError("factorizer exceeded the external safety timeout")
        line = self.process.stdout.readline()
        wall = perf_counter()-began
        cpu_after = cpu_snapshot(self.process.pid)
        if not line:
            raise RuntimeError(f"factorizer exited: {self.process.stderr.read()}")
        result = json.loads(line)
        if result.get("status") == "error":
            raise RuntimeError(result)
        if "process_cpu_seconds" in result:
            cpu, cpu_source = result["process_cpu_seconds"], "internal_process_cpu_clock"
        elif cpu_before is not None and cpu_after is not None:
            cpu, cpu_source = max(0, cpu_after-cpu_before), "procfs_user_system_ticks"
        else:
            cpu, cpu_source = None, "unavailable"
        return result, wall, cpu, cpu_source

    def close(self):
        self.process.stdin.close()
        status = self.process.wait(timeout=15)
        error = self.process.stderr.read()
        if status != 0:
            raise RuntimeError(f"factorizer exited {status}: {error}")


def validate(case, result, workers):
    n = case["n"]
    assert int(result["n"]) == n
    assert result["status"] in {"factor", "budget"}, (case, result)
    if result["status"] == "factor":
        factor = int(result["factor"])
        assert 1 < factor < n and n % factor == 0, (case, result)
        assert int(result["cofactor"]) == n//factor
    else:
        assert result["factor"] is None
    if "outstanding_reservations" in result.get("details", {}):
        assert int(result["details"]["outstanding_reservations"]) == 0
        assert int(result["counters"].get("active_workers_high_water", 0)) <= workers
        counts = result["counters"]
        assert int(counts["jobs_finished"]) == int(counts["jobs_started"])
        assert int(counts["work_units"]) == sum(int(counts.get(key,0)) for key in
                                                ("scalar_tests","represented_candidates","fermat_tests"))
    return result["status"] == "factor"


def paired(rows, reference, comparator, family=None, omit_controls=False):
    a = {r["case_id"]:r for r in rows if r["algorithm"]==reference
         and (family is None or r["family"]==family)
         and (not omit_controls or r["family"]!="close_control")}
    b = {r["case_id"]:r for r in rows if r["algorithm"]==comparator
         and (family is None or r["family"]==family)
         and (not omit_controls or r["family"]!="close_control")}
    keys = [k for k in a if a[k]["solved"] and b[k]["solved"]]
    ratios = [a[k]["median_seconds"]/b[k]["median_seconds"] for k in keys]
    works = [a[k]["median_work_units"]/b[k]["median_work_units"] for k in keys
             if b[k]["median_work_units"] > 0]
    return dict(reference=reference, comparator=comparator, family=family,
                omit_close_controls=omit_controls, matched_successes=len(keys),
                comparator_faster=sum(x>1 for x in ratios),
                median_ratio=statistics.median(ratios) if ratios else None,
                geometric_ratio=statistics.geometric_mean(ratios) if ratios else None,
                median_reference_comparator_work_ratio=statistics.median(works) if works else None,
                comparator_extra=sum(b[k]["solved"] and not a[k]["solved"] for k in a),
                reference_extra=sum(a[k]["solved"] and not b[k]["solved"] for k in a))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT/"build/vf_factor_gmp_parallel")
    parser.add_argument("--sequential-binary", type=Path, default=ROOT/"build/vf_factor_gmp")
    parser.add_argument("--outdir", type=Path, default=ROOT/"results/parallel")
    parser.add_argument("--timeout", type=float, default=1.0)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.timeout <= 0 or args.repeats < 1:
        parser.error("timeout and repeats must be positive")
    binary, sequential_binary = args.binary.resolve(), args.sequential_binary.resolve()
    for path in (binary, sequential_binary):
        if not path.is_file():
            parser.error(f"build executable first: {path}")
    available_cpus = sorted(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else list(range(os.cpu_count() or 1))
    if len(available_cpus) < 4:
        parser.error("this 1/2/4-worker experiment requires at least four exposed logical CPUs")
    out = args.outdir
    out.mkdir(parents=True,exist_ok=True)
    cases = corpus()
    (out/"benchmark-corpus.json").write_text(json.dumps(cases,indent=2)+"\n")
    environment = dict(python=sys.version, platform=platform.platform(),
                       logical_cpu_count=os.cpu_count(), original_affinity=available_cpus,
                       process_cpu_clock="Parallel: internally reported summed process CPU. Sequential: procfs ticks; quantized.",
                       procfs_ticks_per_second=os.sysconf("SC_CLK_TCK") if hasattr(os,"sysconf") else None,
                       cpp_metadata=json.loads(subprocess.check_output([str(binary),"--metadata"],text=True)),
                       sequential_metadata=json.loads(subprocess.check_output([str(sequential_binary),"--metadata"],text=True)),
                       binary_sha256=digest(binary),sequential_binary_sha256=digest(sequential_binary))
    for path in ("/sys/fs/cgroup/cpu.max","/sys/fs/cgroup/cpu.stat","/sys/fs/cgroup/cpuset.cpus.effective","/proc/self/cgroup"):
        p = Path(path)
        environment[path] = p.read_text().strip() if p.is_file() else None
    environment["quota_caveat"] = ("No cgroup quota file is exposed; effective concurrency must be inferred from measured process CPU/wall ratios."
        if environment["/sys/fs/cgroup/cpu.max"] is None else "Cgroup quota is recorded above.")
    source_paths = sorted(list((ROOT/"cpp").glob("*.hpp"))+list((ROOT/"cpp").glob("*.cpp")))
    source_paths += [ROOT/p for p in ("parallel_benchmark.py","Makefile.gmp","vf_factor/fusion_kernel.cpp","vf_factor/corpus.py") if (ROOT/p).is_file()]
    protocol = dict(corpus_seed=CORPUS_SEED,order_seed=ORDER_SEED,
                    timeout_seconds=args.timeout,repeats=args.repeats,
                    regional_options=REGIONAL_CONFIG,methods=METHODS,
                    method_distinction="regional_t1/t2/t4 share the same regional policy. sequential_gmp is a separate older algorithm control, not a scaling denominator.",
                    timing="Internal search wall time is primary. External request wall time and summed process CPU time are retained. Process startup is excluded; per-input worker startup, scheduling and arithmetic setup are included.",
                    affinity="Each worker process is restricted to the first k exposed logical CPUs for k workers; sequential baseline uses one CPU. CPU IDs and actual worker affinity are recorded.",
                    failures="Stop repetitions after first failure. Success requires all requested repetitions. Timed-out cases are not ranked as successful speedups.",
                    work_caveat="Changing worker count can change region order and total work to the first factor. Wall-time ratios are empirical scaling, not proof of equal-work linear speedup.",
                    independence="26 fixed-seed cases are generated before timed calls; factors are used only for independent validation, never passed to the executable.",
                    source_sha256={str(p.relative_to(ROOT)):digest(p) for p in source_paths},
                    corpus_sha256=hashlib.sha256(json.dumps(cases,sort_keys=True).encode()).hexdigest())
    workers = {}
    try:
        for name,count in METHODS.items():
            options = [] if name=="sequential_gmp" else ["--threads",str(count),"--degree",str(REGIONAL_CONFIG["degree"]),
                                                       "--fermat-span",str(REGIONAL_CONFIG["fermat_span"]),
                                                       "--job-seconds",str(REGIONAL_CONFIG["job_seconds"]),
                                                       "--initial-fermat-tests",str(REGIONAL_CONFIG["initial_fermat_tests"]),
                                                       "--fermat-ratio",str(REGIONAL_CONFIG["fermat_ratio"])]
            workers[name] = Worker(sequential_binary if name=="sequential_gmp" else binary,
                                   args.timeout,available_cpus[:count],options)
        environment["actual_worker_affinity"] = {name:w.affinity for name,w in workers.items()}
        (out/"environment.json").write_text(json.dumps(environment,indent=2)+"\n")
        (out/"protocol.json").write_text(json.dumps(protocol,indent=2)+"\n")
        rng=random.Random(ORDER_SEED)
        ordered=list(cases);rng.shuffle(ordered)
        raw, summaries = [], []
        with (out/"benchmark-runs.jsonl").open("w") as output:
            for position,case in enumerate(ordered):
                methods=list(METHODS);rng.shuffle(methods)
                for method in methods:
                    attempts=[]
                    for repeat in range(args.repeats):
                        result,wall,cpu,cpu_source=workers[method].call(case["n"])
                        solved=validate(case,result,METHODS[method])
                        row=dict(case_id=case["case_id"],family=case["family"],bits=case["bits"],
                                 target_bits=case["target_bits"],n=case["n"],algorithm=method,
                                 workers=METHODS[method],repeat=repeat,timeout_seconds=args.timeout,
                                 status=result["status"],factor=int(result["factor"]) if result["factor"] else None,
                                 solved=solved,seconds=result["elapsed_seconds"],wall_seconds=wall,
                                 process_cpu_seconds=cpu,process_cpu_source=cpu_source,
                                 work_units=int(result["counters"].get("work_units",0)),
                                 counters=result["counters"],details=result.get("details",{}))
                        assert row["seconds"]>=0 and wall>=0 and (cpu is None or cpu>=0)
                        output.write(json.dumps(row)+"\n");output.flush()
                        raw.append(row);attempts.append(row)
                        if not solved: break
                    med=statistics.median(r["seconds"] for r in attempts)
                    representative=min(attempts,key=lambda r:abs(r["seconds"]-med))
                    summary={k:representative[k] for k in ("case_id","family","bits","target_bits","n","algorithm","workers","status","factor","counters","details")}
                    cpus=[r["process_cpu_seconds"] for r in attempts if r["process_cpu_seconds"] is not None]
                    summary.update(solved=all(r["solved"] for r in attempts) and len(attempts)==args.repeats,
                                   attempts=len(attempts),requested_repeats=args.repeats,
                                   median_seconds=med,min_seconds=min(r["seconds"] for r in attempts),
                                   max_seconds=max(r["seconds"] for r in attempts),
                                   median_wall_seconds=statistics.median(r["wall_seconds"] for r in attempts),
                                   median_cpu_seconds=statistics.median(cpus) if cpus else None,
                                   median_work_units=statistics.median(r["work_units"] for r in attempts),
                                   mixed_status=len({r["status"] for r in attempts})>1)
                    summaries.append(summary)
                print(json.dumps(dict(completed_inputs=position+1,total_inputs=len(ordered))),flush=True)
    finally:
        for worker in workers.values():
            worker.close()
    summaries.sort(key=lambda r:(r["case_id"],r["algorithm"]))
    (out/"benchmark-summary.json").write_text(json.dumps(summaries,indent=2)+"\n")
    columns=[k for k in summaries[0] if k not in ("counters","details")]
    with (out/"benchmark-summary.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        writer.writerows({k:r[k] for k in columns} for r in summaries)
    comparisons=[]
    for target in ("regional_t2","regional_t4"):
        comparisons.append(paired(summaries,"regional_t1",target))
        comparisons.append(paired(summaries,"regional_t1",target,omit_controls=True))
        for family in ("balanced","small_factor","close_control","hard_balanced"):
            comparisons.append(paired(summaries,"regional_t1",target,family))
    comparisons += [paired(summaries,"sequential_gmp",m) for m in ("regional_t1","regional_t2","regional_t4")]
    outcomes={}
    for method in METHODS:
        rs=[r for r in summaries if r["algorithm"]==method]
        actual=[r for r in raw if r["algorithm"]==method and r["family"]!="close_control"]
        cpus=[r for r in actual if r["process_cpu_seconds"] is not None]
        outcomes[method]=dict(completed=sum(r["solved"] for r in rs),total=len(rs),
                              mixed_repeat_cases=sum(r["mixed_status"] for r in rs),
                              measured_cpu_wall_ratio_excluding_close=sum(r["process_cpu_seconds"] for r in cpus)/sum(r["seconds"] for r in cpus) if cpus else None,
                              maximum_active_workers=max(int(r["counters"].get("active_workers_high_water",0)) for r in rs))
    analysis=dict(inputs=len(cases),summary_rows=len(summaries),timed_calls=len(raw),
                  outcomes=outcomes,paired=comparisons)
    (out/"analysis.json").write_text(json.dumps(analysis,indent=2)+"\n")
    print(json.dumps(analysis,indent=2),flush=True)


if __name__=="__main__":
    main()
