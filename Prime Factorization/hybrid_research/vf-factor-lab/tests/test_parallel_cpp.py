"""Independent checks for exact regional parallel factorization.

The binary receives only N and public resource controls. Tests validate factors
and certified coverage outside the search. No timing speedup is a correctness
assertion, and out-of-order results are never required to have a fixed order.
"""
from __future__ import annotations

import json
from math import isqrt, prod
import os
from pathlib import Path
import random
import subprocess

import pytest
from sympy import isprime, nextprime


ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("VF_PARALLEL_BINARY", ROOT/"build/vf_factor_gmp_parallel"))
pytestmark = pytest.mark.skipif(not BINARY.is_file(), reason="regional parallel executable not built")


def invoke(n, *options, check=True):
    p = subprocess.run([str(BINARY), str(n), *map(str, options)],
                       text=True, capture_output=True, timeout=30)
    if not check:
        return p
    assert p.returncode == 0, (p.returncode, p.stdout, p.stderr)
    return json.loads(p.stdout)


def batch(values, *options, check=True):
    values = list(values)
    p = subprocess.run([str(BINARY), "--batch", *map(str, options)],
                       input="".join(f"{n}\n" for n in values),
                       text=True, capture_output=True,
                       timeout=max(30, len(values)*3))
    if check:
        assert p.returncode == 0, (p.returncode, p.stderr, p.stdout[-2000:])
    rows = [json.loads(line) for line in p.stdout.splitlines()]
    assert len(rows) == len(values), (len(rows), len(values), p.stderr)
    return rows, p


def check_result(n, result, threads, completion=False):
    assert isinstance(result["n"], str) and int(result["n"]) == n
    assert result["status"] in {"factor", "prime", "probable_prime", "budget", "region_exhausted"}, result
    if result["status"] == "factor":
        assert isinstance(result["factor"], str) and isinstance(result["cofactor"], str)
        d = int(result["factor"])
        assert 1 < d < n and d*int(result["cofactor"]) == n, result
    elif result["status"] in {"prime", "probable_prime"}:
        assert isprime(n), result
        assert result["factor"] is None
    elif result["status"] == "region_exhausted":
        assert result["factor"] is None
    else:
        assert not completion, result
        assert result["factor"] is None
    assert result["elapsed_seconds"] >= 0
    assert result["process_cpu_seconds"] >= 0
    counters = result["counters"]
    assert int(counters.get("work_units", 0)) >= 0
    assert 0 <= int(counters.get("active_workers_high_water", 0)) <= threads
    assert int(counters.get("jobs_completed", 0)) <= int(counters.get("jobs_started", 0))
    assert int(counters.get("jobs_cancelled", 0)) <= int(counters.get("jobs_started", 0))
    assert int(result["details"]["outstanding_reservations"]) == 0


def check_trace(n, smaller_factor, result):
    """Every accepted certificate must exclude no actual smaller factor."""
    trace = result["details"]["trace"]
    previous_work = 0
    job_ids = []
    domain = None if trace and trace[0]["event"] == "seed" else (7, isqrt(n))
    covered = []
    for event in trace:
        assert event["event"] in {"seed", "job_result"}, event
        work = int(event["work_units"])
        assert previous_work <= work <= int(result["counters"]["work_units"])
        previous_work = work
        if event.get("lo") is not None and event.get("hi") is not None:
            lo, hi = int(event["lo"]), int(event["hi"])
            assert lo <= smaller_factor <= hi, event
            if domain is None:
                # The serial seed is independently validated. Reconstruct all
                # subsequent asynchronous exclusions from their certificates.
                domain = (lo, hi)
        certificate = (event.get("certified_lo"), event.get("certified_hi"))
        if all(x is not None for x in certificate):
            lo, hi = map(int, certificate)
            assert lo <= hi
            assert not lo <= smaller_factor <= hi, event
            if domain is not None and event["event"] == "job_result":
                lo, hi = max(lo, domain[0]), min(hi, domain[1])
                if lo <= hi:
                    covered.append((lo, hi))
                    merged = []
                    for left, right in sorted(covered):
                        if merged and left <= merged[-1][1]+1:
                            merged[-1] = (merged[-1][0], max(right, merged[-1][1]))
                        else:
                            merged.append((left, right))
                    covered = merged
        else:
            assert certificate == (None, None), event
        if event["event"] == "job_result":
            job_ids.append(event["job_id"])
            assert event["kind"] in {"product", "fermat"}
            assert int(event["requested_lo"]) <= int(event["requested_hi"])
            assert int(event["job_work"]) >= 0
            if domain is not None and event.get("lo") is not None and event.get("hi") is not None:
                expected_lo, expected_hi = domain
                for left, right in covered:
                    if left <= expected_lo <= right:
                        expected_lo = right+1
                for left, right in reversed(covered):
                    if left <= expected_hi <= right:
                        expected_hi = left-1
                # Conservative frontiers are allowed. Advancing beyond an
                # uncertified hole is not, even when jobs finish out of order.
                assert int(event["lo"]) <= expected_lo, event
                assert int(event["hi"]) >= expected_hi, event
    assert len(job_ids) == len(set(job_ids)), job_ids
    assert len(job_ids) == int(result["counters"]["jobs_started"])
    if trace:
        assert int(trace[-1]["work_units"]) == int(result["counters"]["work_units"])


def test_metadata_identifies_parallel_backend():
    p = subprocess.run([str(BINARY), "--metadata"], text=True,
                       capture_output=True, check=True, timeout=5)
    metadata = json.loads(p.stdout)
    assert "parallel" in json.dumps(metadata).lower()
    assert "gmp" in json.dumps(metadata).lower()


@pytest.mark.parametrize("threads", (1, 2, 4))
def test_small_input_truth_and_prime_classification(threads):
    values = list(range(2, 502))
    rows, _ = batch(values, "--threads", threads, "--timeout", 1)
    for n, row in zip(values, rows):
        check_result(n, row, threads, completion=True)
        assert (row["status"] == "prime") == bool(isprime(n)), (n, row)


@pytest.mark.parametrize("threads", (1, 2, 4))
@pytest.mark.parametrize("options", [(), ("--no-vf",), ("--gmp-fermat",),
                                     ("--no-fermat", "--no-vf")])
def test_parallel_regions_return_only_verified_factors(threads, options):
    cases = [551, 1351, 8251, 978508015703, 1009*1000000007, 1009**3,
             1009*int(nextprime(1 << 117))]
    if "--no-fermat" in options:
        cases.remove(978508015703)
    rows, _ = batch(cases, "--threads", threads, "--timeout", 2, *options)
    for n, row in zip(cases, rows):
        check_result(n, row, threads, completion=True)


def test_random_interleavings_and_batch_state_are_independent():
    rng = random.Random(926503)
    cases = []
    for _ in range(24):
        p = int(nextprime(rng.randrange(257, 8000)))
        q = int(nextprime(rng.randrange(1 << 30, 1 << 32)))
        cases.append(p*q)
    rows, _ = batch(cases, "--threads", 4, "--timeout", 2)
    for n, row in zip(cases, rows):
        check_result(n, row, 4, completion=True)


@pytest.mark.parametrize("bits", (80, 128, 256, 512))
def test_large_arithmetic_stays_exact(bits):
    p = int(nextprime((3 << (bits//2-2)) + 1433))
    q = int(nextprime(p+42))
    n = p*q
    row = invoke(n, "--threads", 4, "--timeout", 1)
    check_result(n, row, 4, completion=True)
    assert sorted([int(row["factor"]), int(row["cofactor"])]) == [p, q]


@pytest.mark.parametrize("threads", (1, 2, 4))
@pytest.mark.parametrize("limit", (0, 1, 63, 64, 65, 100, 1000, 10000, 100000))
def test_simultaneous_work_reservations_respect_one_cap(threads, limit):
    n = 1000003*int(nextprime(1 << 100))
    row = invoke(n, "--threads", threads, "--timeout", 0.25,
                 "--max-steps", limit, "--trace")
    check_result(n, row, threads)
    assert int(row["counters"]["work_units"]) <= limit
    assert int(row["counters"].get("reservation_high_water", 0)) <= limit
    check_trace(n, 1000003, row)


@pytest.mark.parametrize("options", [(), ("--no-fermat", "--no-vf"), ("--gmp-fermat",)])
def test_out_of_order_certificates_preserve_the_actual_factor(options):
    p = 1000003
    n = p*int(nextprime(1 << 100))
    row = invoke(n, "--threads", 4, "--timeout", 0.2, "--trace", *options)
    check_result(n, row, 4)
    check_trace(n, p, row)
    assert int(row["counters"]["jobs_started"]) > 0
    for event in row["details"]["trace"]:
        if event["event"] == "job_result" and event["kind"] == "product":
            if event.get("certified_lo") is not None:
                assert int(event["requested_lo"]) <= int(event["certified_lo"])
                assert int(event["certified_hi"]) <= int(event["requested_hi"])


def test_multiple_worker_execution_and_cooperative_cancellation():
    n = int(nextprime(1 << 63))*int(nextprime(3 << 62))
    row = invoke(n, "--threads", 4, "--timeout", 0.15, "--trace")
    check_result(n, row, 4)
    assert row["status"] == "budget"
    assert 2 <= int(row["counters"]["active_workers_high_water"]) <= 4
    assert int(row["details"]["outstanding_reservations"]) == 0
    check_trace(n, int(nextprime(1 << 63)), row)


def test_full_factorization_uses_one_budget_and_preserves_multiplicities():
    n = 1009*1013*10007*10009
    row = invoke(n, "--threads", 4, "--full", "--timeout", 2, "--max-steps", 100)
    assert int(row["input"]) == n
    groups = ("prime_factors", "probable_prime_factors", "unresolved")
    assert prod(int(x) for key in groups for x in row[key]) == n
    assert all(isprime(int(x)) for x in row["prime_factors"])
    assert all(isprime(int(x)) for x in row["probable_prime_factors"])
    assert int(row["work_units"]) <= 100
    assert int(row["work_units"]) == sum(int(r["counters"]["work_units"]) for r in row["runs"])
    assert row["process_cpu_seconds"] >= 0
    for run in row["runs"]:
        check_result(int(run["n"]), run, 4)


def test_full_factorization_labels_probable_prime_leaves():
    large_prime = (1 << 127)-1
    n = 2**4*3**2*101**2*large_prime
    row = invoke(n, "--threads", 4, "--full", "--timeout", 2)
    assert row["unresolved"] == []
    assert sorted(map(int, row["prime_factors"])) == [2]*4+[3]*2+[101]*2
    assert list(map(int, row["probable_prime_factors"])) == [large_prime]
    assert prod(int(x) for key in ("prime_factors", "probable_prime_factors") for x in row[key]) == n


def test_zero_deadline_and_full_zero_deadline_remain_unresolved():
    n = int(nextprime(1 << 63))*int(nextprime(3 << 62))
    row = invoke(n, "--threads", 4, "--timeout", 0)
    check_result(n, row, 4)
    assert row["status"] == "budget"
    assert row["counters"]["work_units"] == 0
    full = invoke(n, "--threads", 4, "--full", "--timeout", 0)
    assert full["prime_factors"] == []
    assert full["probable_prime_factors"] == []
    assert list(map(int, full["unresolved"])) == [n]
    assert int(full["work_units"]) == 0


def test_batch_error_and_budget_do_not_poison_next_input():
    hard = int(nextprime(1 << 63))*int(nextprime(3 << 62))
    rows, process = batch([hard, "invalid", 551, 8051], "--threads", 4,
                          "--max-steps", 100, "--timeout", 1, check=False)
    assert process.returncode != 0
    assert rows[0]["status"] == "budget"
    assert rows[1]["status"] == "error"
    for n, row in zip((551, 8051), rows[2:]):
        check_result(n, row, 4, completion=True)
        assert int(row["counters"]["work_units"]) <= 100


def test_completed_partial_region_is_not_a_prime_claim():
    n = 1009*1000000007
    row = invoke(n, "--threads", 4, "--factor-lo", 2000, "--factor-hi", 2200,
                 "--no-fermat", "--no-vf", "--timeout", 2, "--trace")
    check_result(n, row, 4, completion=True)
    assert row["status"] == "region_exhausted"
    assert row["factor"] is None
    for event in row["details"]["trace"]:
        if event.get("certified_lo") is not None:
            assert 2000 <= int(event["certified_lo"]) <= int(event["certified_hi"]) <= 2200


@pytest.mark.parametrize("lower,upper", [(900, 1100), (1009, 1009), (1000, 1009)])
def test_explicit_factor_region_contains_a_factor(lower, upper):
    n = 1009*1000000007
    row = invoke(n, "--threads", 2, "--factor-lo", lower, "--factor-hi", upper,
                 "--no-fermat", "--no-vf", "--timeout", 2)
    check_result(n, row, 2, completion=True)
    assert row["status"] == "factor"
    assert int(row["factor"]) == 1009


def test_partial_region_cannot_be_used_as_full_factorization():
    process = invoke(1009*1000000007, "--full", "--factor-lo", 900,
                     "--factor-hi", 1100, check=False)
    assert process.returncode != 0


def test_inverted_region_is_rejected():
    process = invoke(1009*1000000007, "--factor-lo", 2200,
                     "--factor-hi", 2000, check=False)
    assert process.returncode != 0


@pytest.mark.parametrize("options", [("--threads", "0"), ("--threads", "-1"),
                                     ("--threads", "abc"), ("--timeout", "nan"),
                                     ("--max-steps", "-1")])
def test_invalid_parallel_controls_are_rejected(options):
    process = invoke(77, *options, check=False)
    assert process.returncode != 0
    assert process.stderr.strip()
