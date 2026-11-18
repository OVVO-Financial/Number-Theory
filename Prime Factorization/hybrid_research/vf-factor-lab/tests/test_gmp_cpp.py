"""Independent integration checks for the standalone GMP C++ factorizer.

Build with ``python3 build_gmp.py`` before running this module.  Reference
factors and primality information remain in the test process; the executable
receives only the integer and public resource controls.  Timings are not used
as correctness assertions.
"""
from __future__ import annotations

import json
from math import prod
import os
from pathlib import Path
import random
import subprocess

import pytest
from sympy import isprime, nextprime


ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("VF_GMP_BINARY", ROOT / "build/vf_factor_gmp"))
pytestmark = pytest.mark.skipif(
    not BINARY.is_file(), reason="GMP executable not built; run python3 build_gmp.py"
)


def invoke(n, *options, check=True):
    proc = subprocess.run(
        [str(BINARY), str(n), *map(str, options)], text=True,
        capture_output=True, timeout=20,
    )
    if not check:
        return proc
    assert proc.returncode == 0, (proc.returncode, proc.stdout, proc.stderr)
    return json.loads(proc.stdout)


def batch(values, *options, check=True):
    values = list(values)
    proc = subprocess.run(
        [str(BINARY), "--batch", *map(str, options)],
        input="".join(f"{n}\n" for n in values), text=True,
        capture_output=True, timeout=max(30, len(values) * 2),
    )
    if check:
        assert proc.returncode == 0, (proc.returncode, proc.stderr, proc.stdout[-1000:])
    lines = proc.stdout.splitlines()
    assert len(lines) == len(values), (len(lines), len(values), proc.stderr)
    return [json.loads(line) for line in lines], proc


def check_result(n, result, *, require_completion=False):
    assert isinstance(result["n"], str) and int(result["n"]) == n
    status = result["status"]
    assert status in {"factor", "prime", "probable_prime", "budget"}, result
    if status == "factor":
        assert isinstance(result["factor"], str)
        assert isinstance(result["cofactor"], str)
        factor, cofactor = int(result["factor"]), int(result["cofactor"])
        assert 1 < factor < n and factor * cofactor == n, result
    elif status in {"prime", "probable_prime"}:
        assert isprime(n), result
        assert result["factor"] is None and result["cofactor"] is None
    else:
        assert result["factor"] is None
        assert not require_completion, result
    assert result["elapsed_seconds"] >= 0
    assert result["counters"].get("work_units", 0) >= 0


def check_full(n, result):
    assert isinstance(result["input"], str) and int(result["input"]) == n
    fields = ("prime_factors", "probable_prime_factors", "unresolved")
    for key in fields:
        assert all(isinstance(x, str) for x in result[key])
    factors = [int(x) for key in fields for x in result[key]]
    assert all(x > 1 for x in factors)
    assert prod(factors) == n
    assert all(isprime(int(p)) for p in result["prime_factors"])
    assert all(isprime(int(p)) for p in result["probable_prime_factors"])
    for run in result["runs"]:
        check_result(int(run["n"]), run)
    assert result["work_units"] == sum(r["counters"].get("work_units", 0)
                                        for r in result["runs"])


def test_metadata_identifies_gmp_and_native_build():
    proc = subprocess.run([str(BINARY), "--metadata"], text=True,
                          capture_output=True, check=True, timeout=5)
    meta = json.loads(proc.stdout)
    assert "gmp" in json.dumps(meta).lower()
    assert "compiler" in json.dumps(meta).lower()


def test_all_small_inputs_independently_checked():
    values = list(range(2, 1002))
    results, _ = batch(values, "--timeout", 1)
    for n, result in zip(values, results):
        check_result(n, result, require_completion=True)
        assert (result["status"] == "prime") == bool(isprime(n)), (n, result)


@pytest.mark.parametrize("options", [(), ("--gmp-fermat",), ("--no-vf",),
                                     ("--no-feedback",), ("--no-sharing",),
                                     ("--no-fermat",),
                                     ("--no-fermat", "--no-vf")])
def test_residue_regressions_and_mixed_search_options(options):
    values = [511, 551, 1351, 2291, 8251, 8051, 798607,
              1009 * 1000000007, 1009**3, 4294967279 * 4294967291]
    # The near-2^64 close pair belongs to Fermat tests; its factors are too
    # large for a deliberately Fermat-disabled bounded product search.
    if "--no-fermat" in options:
        values.pop()
    results, _ = batch(values, "--timeout", 2, *options)
    for n, result in zip(values, results):
        check_result(n, result, require_completion=True)


def test_arbitrary_precision_close_pairs_and_small_factors():
    rng = random.Random(246811)
    cases = []
    for bits in (80, 96, 128, 192, 256, 384, 512):
        p = int(nextprime((3 << (bits // 2 - 2)) + rng.randrange(1000, 100000)))
        q = int(nextprime(p + rng.randrange(2, 200)))
        cases.append(p * q)
        qlarge = int(nextprime((1 << (bits - 10)) + rng.randrange(1, 10000)))
        cases.append(1009 * qlarge)
    results, _ = batch(cases, "--timeout", 2)
    for n, result in zip(cases, results):
        assert n > 1 << 64
        check_result(n, result, require_completion=True)


@pytest.mark.parametrize("bits", (65, 80, 128, 256, 512))
def test_large_square_precheck(bits):
    p = int(nextprime((1 << bits) + 1009))
    result = invoke(p*p, "--timeout", 1)
    check_result(p*p, result, require_completion=True)
    assert int(result["factor"]) == p
    assert result["counters"].get("work_units", 0) == 0


def test_proven_and_probable_primes_have_distinct_statuses():
    values = [2, 3, 65537, 18446744073709551557, (1 << 127) - 1, (1 << 521) - 1]
    results, _ = batch(values, "--timeout", 1)
    for n, result in zip(values, results):
        check_result(n, result, require_completion=True)
        if n < 1 << 64:
            assert result["status"] == "prime"
        else:
            assert result["status"] == "probable_prime"


def test_uint64_strong_pseudoprimes_are_not_certified_prime():
    values = [3215031751, 2152302898747, 341550071728321, 3825123056546413051]
    assert all(not isprime(n) for n in values)
    results, _ = batch(values, "--timeout", 0)
    for result in results:
        assert result["status"] == "budget", result


def test_disabling_probable_prime_shortcut_preserves_unresolved_status():
    n = (1 << 127) - 1
    result = invoke(n, "--no-probable-prime-check", "--no-fermat", "--max-steps", 0)
    assert result["status"] == "budget"
    assert result["counters"].get("work_units", 0) == 0


def test_full_factorization_preserves_all_multiplicities_and_probable_leaf():
    large_prime = (1 << 127) - 1
    n = 2**5 * 3**3 * 101**2 * large_prime
    result = invoke(n, "--full", "--timeout", 2)
    check_full(n, result)
    assert result["unresolved"] == []
    assert sorted(map(int, result["prime_factors"])) == [2]*5 + [3]*3 + [101]*2
    assert list(map(int, result["probable_prime_factors"])) == [large_prime]


@pytest.mark.parametrize("limit", (0, 1, 7, 64, 100, 1000, 5000))
def test_work_cap_is_shared_across_search_components(limit):
    n = int(nextprime(1 << 100)) * 1000003
    result = invoke(n, "--timeout", 2, "--max-steps", limit, "--trace")
    check_result(n, result)
    assert result["status"] == "budget"
    assert result["counters"].get("work_units", 0) <= limit
    for event in result["details"].get("trace", []):
        assert int(event["work_units"]) <= limit


def test_full_recursion_uses_one_work_cap():
    n = 1009 * 1013 * 10007 * 10009
    result = invoke(n, "--full", "--timeout", 2, "--max-steps", 100)
    check_full(n, result)
    assert len(result["runs"]) >= 3  # The first Fermat split creates two searches.
    assert result["work_units"] <= 100


def test_zero_timeout_does_not_claim_a_hard_composite_is_prime():
    n = int(nextprime(1 << 100)) * int(nextprime(3 << 99))
    result = invoke(n, "--timeout", 0)
    assert result["status"] == "budget"
    assert result["factor"] is None
    assert result["counters"].get("work_units", 0) == 0


def test_batch_resets_budget_for_each_input():
    values = [int(nextprime(1 << 100)) * 1000003, 551, 8051]
    results, _ = batch(values, "--max-steps", 100, "--timeout", 1)
    assert results[0]["status"] == "budget"
    for n, result in zip(values[1:], results[1:]):
        check_result(n, result, require_completion=True)
    assert all(r["counters"].get("work_units", 0) <= 100 for r in results)


def test_bound_trace_preserves_actual_factor_pair():
    n, p = 978508015703, 752867
    result = invoke(n, "--timeout", 2, "--trace")
    check_result(n, result, require_completion=True)
    trace = result["details"]["trace"]
    assert result["counters"].get("mixed_batches", 0) > 0
    assert result["counters"].get("upper_bound_updates", 0) > 0
    assert result["counters"].get("lower_bound_updates", 0) > 0
    for event in trace:
        lo, hi = int(event["lo"]), int(event["hi"])
        assert lo <= p <= hi, event
        if event.get("fermat_cap") is not None:
            assert int(event["fermat_cap"]) >= (p + n//p)//2, event


def test_deterministic_product_trace_matches_python_under_work_cap():
    from vf_factor.fused import factor
    n = int(nextprime(1 << 100)) * 1000003
    limit = 6000
    cpp = invoke(n, "--timeout", 5, "--max-steps", limit,
                 "--no-fermat", "--no-vf", "--trace")
    py = factor(n, timeout=5, max_steps=limit, fermat_enabled=False,
                vf_probes=False, trace=True).to_dict()
    assert cpp["status"] == py["status"] == "budget"
    for key in ("scalar_tests", "represented_candidates", "divisor_candidates",
                "polynomials_built", "evaluations", "work_units"):
        assert cpp["counters"].get(key, 0) == py["counters"].get(key, 0), key
    assert int(cpp["details"]["lo"]) == py["details"]["lo"]
    assert int(cpp["details"]["hi"]) == py["details"]["hi"]
    def canonical(trace):
        return [(e["event"], int(e["lo"]), int(e["hi"]), int(e["work_units"]))
                for e in trace]
    assert canonical(cpp["details"]["trace"]) == canonical(py["details"]["trace"])


@pytest.mark.parametrize("n", ("0", "1", "-7", "77x", "abc", "0x123", ""))
def test_invalid_decimal_inputs_are_rejected(n):
    result = invoke(n, check=False)
    assert result.returncode != 0
    assert result.stderr.strip()


@pytest.mark.parametrize("options", [
    ("--timeout", "-1"), ("--timeout", "nan"), ("--timeout", "inf"),
    ("--max-steps", "-1"), ("--max-steps", str(1 << 64)),
    ("--initial-degree", "0"), ("--max-degree", "0"),
    ("--initial-degree", "32", "--max-degree", "16"),
    ("--fermat-ratio", "0"), ("--fermat-ratio", "nan"),
    ("--initial-fermat-tests", "-1"), ("--not-an-option",),
])
def test_invalid_controls_are_rejected(options):
    result = invoke(551, *options, check=False)
    assert result.returncode != 0
    assert result.stderr.strip()


def test_batch_errors_are_per_line_and_do_not_drop_later_requests():
    values = ["551", "wrong", "", "1", "8051"]
    results, proc = batch(values, "--timeout", 1, check=False)
    assert proc.returncode != 0
    check_result(551, results[0], require_completion=True)
    for result in results[1:4]:
        assert result["status"] == "error" and result["error"]
    check_result(8051, results[4], require_completion=True)
