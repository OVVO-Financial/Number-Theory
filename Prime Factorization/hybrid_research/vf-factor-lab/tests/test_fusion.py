"""Independent checks for the single shared-bound factorization algorithm.

References enumerate integer factor pairs directly. No known factors are
passed to search functions. The cursor tests verify certificates and resumption,
not just the validity of a returned divisor; timing thresholds are not asserted.
"""
from __future__ import annotations

from math import gcd, isqrt, prod
import random

import pytest
from sympy import isprime, nextprime

from vf_factor import fusion_native


NATIVE = pytest.mark.skipif(
    not fusion_native.is_available(),
    reason="optional fusion kernel not built; run python3 build_fusion.py",
)
QR20 = {r * r % 20 for r in range(20)}
QR64 = {r * r % 64 for r in range(64)}


def ceil_root(n: int) -> int:
    root = isqrt(n)
    return root + (root * root != n)


def eligible_midpoint(n: int, a: int) -> bool:
    """Independent residue-square criterion, without a cyclic gap table."""
    difference = a * a - n
    return difference % 20 in QR20 and difference % 64 in QR64


def next_eligible(n: int, start: int) -> int:
    return next(a for a in range(start, start + 160)
                if eligible_midpoint(n, a))


def direct_fermat(n: int) -> tuple[int | None, int, int]:
    """Scan every integer midpoint, without either wheel or shared bounds."""
    for a in range(ceil_root(n), (n + 1) // 2 + 1):
        b = isqrt(a * a - n)
        if b * b == a * a - n:
            factor = a - b
            if factor > 1:
                return factor, a, b
            return None, a, b
    raise AssertionError("Every odd integer has at least its trivial Fermat pair")


def validate_certificate(n: int, snapshot: dict) -> None:
    a0 = ceil_root(n)
    covered = snapshot["covered_through"]
    # Independent check over all actually excluded integer midpoints.
    for a in range(a0, covered + 1):
        b = isqrt(a * a - n)
        assert b * b != a * a - n or a - b <= 1, (n, snapshot, a)
    assert snapshot["candidates"] == snapshot["square_tests"]
    if snapshot["factor"]:
        assert 1 < snapshot["factor"] < n
        assert n % snapshot["factor"] == 0
        assert covered < snapshot["last_midpoint"]
    elif not snapshot["exhausted"]:
        assert eligible_midpoint(n, snapshot["next_midpoint"])
        assert covered == snapshot["next_midpoint"] - 1


@NATIVE
@pytest.mark.parametrize("n", (3, 5, 7, 9, 15, 25, 35, 49, 77, 101, 143,
                                    343, 551, 8051, 101 * 5003))
def test_cursor_setup_certificate_and_exact_first_candidate(n: int) -> None:
    with fusion_native.FermatCursor(n) as cursor:
        state = cursor.snapshot()
        expected = next_eligible(n, ceil_root(n))
        assert state["status"] == "ready"
        assert state["root_ceiling"] == ceil_root(n)
        assert state["endpoint"] == (n + 1) // 2
        assert state["next_midpoint"] == expected
        assert state["covered_through"] == expected - 1
        assert state["candidates"] == 0
        assert state["first_midpoint"] is None
        assert state["last_midpoint"] is None
        validate_certificate(n, state)


@NATIVE
def test_cursor_all_small_odd_inputs_without_prechecks() -> None:
    for n in range(3, 1002, 2):
        factor, a, b = direct_fermat(n)
        expected_candidates = sum(eligible_midpoint(n, x)
                                  for x in range(ceil_root(n), a + 1))
        with fusion_native.FermatCursor(n) as cursor:
            state = cursor.advance(1000)
        assert state["factor"] == factor, (n, state)
        assert state["last_midpoint"] == a
        assert state["candidates"] == expected_candidates
        assert state["delta_candidates"] == expected_candidates
        if factor:
            assert state["status"] == "factor"
            assert state["half_gap"] == b
        else:
            assert state["status"] == "exhausted"
            assert state["exhausted"]
        validate_certificate(n, state)


@NATIVE
@pytest.mark.parametrize("n", (551, 8051, 101 * 5003, 7**5, 1009 * 3011, 1009))
def test_cursor_resumption_equals_one_continuous_scan(n: int) -> None:
    with fusion_native.FermatCursor(n) as whole:
        expected = whole.advance(1_000_000)
    with fusion_native.FermatCursor(n) as resumed:
        prior_count = prior_coverage = 0
        for iteration in range(1000):
            state = resumed.advance((1, 3, 2, 7, 19)[iteration % 5])
            assert state["delta_candidates"] == state["candidates"] - prior_count
            assert state["covered_through"] >= prior_coverage
            assert state["candidates"] >= prior_count
            validate_certificate(n, state)
            prior_count = state["candidates"]
            prior_coverage = state["covered_through"]
            if state["factor"] or state["exhausted"]:
                break
        else:
            pytest.fail("Resumed cursor did not finish its small test input")
        for key in ("factor", "candidates", "square_tests", "first_midpoint",
                    "last_midpoint", "next_midpoint", "covered_through", "half_gap",
                    "status", "exhausted"):
            assert state[key] == expected[key], (n, key, state, expected)
        repeated = resumed.advance(10)
        assert repeated["delta_candidates"] == 0
        assert repeated["candidates"] == state["candidates"]
        assert repeated["factor"] == state["factor"]


@NATIVE
@pytest.mark.parametrize("kind", ("zero_work", "zero_time"))
def test_cursor_zero_budget_does_not_advance_and_can_resume(kind: str) -> None:
    with fusion_native.FermatCursor(101 * 5003) as cursor:
        initial = cursor.snapshot()
        if kind == "zero_work":
            state = cursor.advance(0)
            assert state["status"] == "candidate_budget"
        else:
            state = cursor.advance(10, timeout_seconds=0)
            assert state["status"] == "time_budget"
        for key in ("factor", "candidates", "square_tests", "next_midpoint",
                    "covered_through", "first_midpoint", "last_midpoint"):
            assert state[key] == initial[key]
        assert state["delta_candidates"] == 0
        continued = cursor.advance(1)
        assert continued["delta_candidates"] == 1
        assert continued["candidates"] == 1


@NATIVE
def test_temporary_midpoint_cap_does_not_destroy_cursor_state() -> None:
    n = 101 * 5003
    factor, solution, _ = direct_fermat(n)
    with fusion_native.FermatCursor(n) as cursor:
        initial = cursor.snapshot()
        before_first = cursor.advance(1000, max_midpoint=initial["next_midpoint"] - 1)
        assert before_first["status"] == "range_limit"
        assert before_first["candidates"] == 0
        assert not before_first["exhausted"]
        assert before_first["next_midpoint"] == initial["next_midpoint"]
        partial = cursor.advance(1000, max_midpoint=solution - 1)
        assert partial["factor"] is None
        assert partial["status"] == "range_limit"
        assert not partial["exhausted"]
        assert partial["last_midpoint"] < solution
        assert partial["covered_through"] == solution - 1
        validate_certificate(n, partial)
        completed = cursor.advance(1000, max_midpoint=solution)
        assert completed["factor"] == factor
        assert completed["last_midpoint"] == solution
        assert completed["delta_candidates"] == 1


@NATIVE
@pytest.mark.parametrize("cap_offset", range(12))
def test_range_cap_never_performs_a_square_check_beyond_cap(cap_offset: int) -> None:
    n = 101 * 5003
    cap = ceil_root(n) + cap_offset
    expected = [a for a in range(ceil_root(n), cap + 1)
                if eligible_midpoint(n, a)]
    with fusion_native.FermatCursor(n) as cursor:
        state = cursor.advance(1000, max_midpoint=cap)
    assert state["factor"] is None
    assert state["candidates"] == len(expected)
    assert state["last_midpoint"] == (expected[-1] if expected else None)
    assert state["next_midpoint"] > cap
    validate_certificate(n, state)


@NATIVE
@pytest.mark.parametrize("n", (4294967291 * 4294967279,
                                    (2**32 + 1)**2 - 1000014**2,
                                    2**64 - 1))
def test_cursor_exact_arithmetic_at_uint64_boundary(n: int) -> None:
    factor, a, b = direct_fermat(n)
    with fusion_native.FermatCursor(n) as cursor:
        state = cursor.advance(100_000)
    assert state["factor"] == factor
    assert state["last_midpoint"] == a
    assert state["half_gap"] == b
    assert state["covered_through"] < a
    assert state["candidates"] == sum(eligible_midpoint(n, x)
                                      for x in range(ceil_root(n), a + 1))


@NATIVE
@pytest.mark.parametrize("n", (-1, 0, 1, 2, 4, True, 77.0, "77", 2**64, 2**64 + 1))
def test_cursor_rejects_out_of_domain_inputs(n) -> None:
    with pytest.raises(ValueError):
        fusion_native.FermatCursor(n)


@NATIVE
@pytest.mark.parametrize("parameter,value", (
    ("max_candidates", -1), ("max_candidates", True), ("max_candidates", 1.5),
    ("max_candidates", 2**64), ("max_midpoint", -1), ("max_midpoint", True),
    ("max_midpoint", 2**64), ("timeout_seconds", -1), ("timeout_seconds", True),
    ("timeout_seconds", float("nan")), ("timeout_seconds", float("inf")),
))
def test_cursor_rejects_invalid_controls(parameter: str, value) -> None:
    kwargs = {"max_candidates": 1, parameter: value}
    with fusion_native.FermatCursor(8051) as cursor:
        with pytest.raises(ValueError):
            cursor.advance(**kwargs)


@NATIVE
def test_cursor_closure_is_idempotent_and_closed_state_is_unusable() -> None:
    cursor = fusion_native.FermatCursor(8051)
    cursor.close()
    cursor.close()
    with pytest.raises(RuntimeError, match="closed"):
        cursor.snapshot()
    with pytest.raises(RuntimeError, match="closed"):
        cursor.advance(1)


@pytest.mark.parametrize("modulus", (8, 77, 143, 8051, 2**64 + 1, 2**80 + 7))
@pytest.mark.parametrize("degree", (1, 2, 7, 16, 31))
@pytest.mark.parametrize("shared", (False, True))
def test_shared_polynomial_matches_independent_mixed_step_products(
    modulus: int, degree: int, shared: bool,
) -> None:
    from vf_factor.shared_products import NormalizedProductEngine

    rng = random.Random(20261118 + modulus + degree)
    units = [s for s in (1, -1, 2, -2, 3, 5, 7, 30, modulus + 2)
             if gcd(s, modulus) == 1]
    starts = [0, modulus, -modulus, -1, 7, 11]
    starts += [rng.randrange(-2 * modulus, 2 * modulus + 1) for _ in range(9)]
    steps = [units[i % len(units)] for i in range(len(starts))]
    counters = {}
    engine = NormalizedProductEngine(modulus, degree, counters=counters, shared=shared)
    values = engine.products(starts, steps)
    assert len(values) == len(starts)
    for start, step, value in zip(starts, steps, values):
        expected = prod(start + step * i for i in range(degree)) % modulus
        if shared:
            assert int(value) * pow(step, degree, modulus) % modulus == expected
        else:
            assert int(value) == expected
        assert gcd(int(value), modulus) == gcd(expected, modulus)
    # Reusing an engine must share polynomial construction, not merely return
    # equivalent arithmetic once. Repeat in reverse order to check association.
    again = engine.products(list(reversed(starts)), list(reversed(steps)))
    assert list(map(int, again)) == list(reversed(list(map(int, values))))
    expected_builds = 1 if shared else len({s % modulus for s in steps})
    assert counters["polynomials_built"] == expected_builds


@pytest.mark.parametrize("shared", (False, True))
@pytest.mark.parametrize("step", (0, 7, -7, 11, 77, 154))
def test_shared_normalization_rejects_nonunit_steps(shared: bool, step: int) -> None:
    from vf_factor.shared_products import NormalizedProductEngine

    engine = NormalizedProductEngine(77, 4, shared=shared)
    with pytest.raises(ValueError):
        engine.products([1], [step])


@pytest.mark.parametrize("shared", (False, True))
def test_shared_product_empty_and_misaligned_batches(shared: bool) -> None:
    from vf_factor.shared_products import NormalizedProductEngine

    engine = NormalizedProductEngine(77, 4, shared=shared)
    assert engine.products([], []) == []
    with pytest.raises(ValueError):
        engine.products([1, 2], [5])
    with pytest.raises(ValueError):
        engine.products([1], [2, 5])


@pytest.mark.parametrize("shared", (False, True))
@pytest.mark.parametrize("modulus,degree,starts,steps,expected", (
    # A VF block containing only a trivial GCD=N must not hide a later factor.
    (77, 4, (77, 7), (5, 2), 7),
    (143, 4, (143, 11), (30, 2), 11),
    # Two copies of a prime make the block product vanish modulo a prime power.
    (49, 8, (7,), (5,), 7),
    # Signed unit steps must also recover the true scalar factors correctly.
    (77, 8, (7,), (-1,), 7),
))
def test_mixed_block_recovery_handles_gcd_n_and_trivial_hits(
    shared: bool, modulus: int, degree: int, starts: tuple[int, ...],
    steps: tuple[int, ...], expected: int,
) -> None:
    from vf_factor.shared_products import NormalizedProductEngine, recover_mixed_batch

    counters = {}
    engine = NormalizedProductEngine(modulus, degree, counters=counters, shared=shared)
    values = engine.products(list(starts), list(steps))
    requests = [dict(start=start, step=step, count=degree,
                     kind="vf" if index == 0 else "trial")
                for index, (start, step) in enumerate(zip(starts, steps))]
    assert gcd(prod(values), modulus) == modulus
    factor = recover_mixed_batch(modulus, requests, values, counters=counters)
    assert factor == expected
    assert counters["batch_gcd_calls"] == 1


@pytest.mark.parametrize("shared", (False, True))
def test_mixed_batch_with_only_trivial_hit_does_not_invent_factor(shared: bool) -> None:
    from vf_factor.shared_products import NormalizedProductEngine, recover_mixed_batch

    engine = NormalizedProductEngine(77, 4, shared=shared)
    values = engine.products([77], [5])
    result = recover_mixed_batch(77, [dict(start=77, step=5, count=4, kind="vf")], values)
    assert result is None


def test_all_coprime_mixed_blocks_use_one_aggregate_gcd() -> None:
    from vf_factor.shared_products import NormalizedProductEngine, recover_mixed_batch

    engine = NormalizedProductEngine(143, 3)
    starts, steps = [1, 2, 31], [2, 3, 5]
    values = engine.products(starts, steps)
    assert all(gcd(v, 143) == 1 for v in values)
    requests = [dict(start=a, step=s, count=3, kind="trial") for a, s in zip(starts, steps)]
    counters = {}
    assert recover_mixed_batch(143, requests, values, counters=counters) is None
    assert counters["gcd_calls"] == counters["batch_gcd_calls"] == 1


def test_shared_normalized_search_honors_expired_budget() -> None:
    from vf_factor.common import Budget, BudgetExceeded
    from vf_factor.shared_products import NormalizedProductEngine

    engine = NormalizedProductEngine(8051, 8, budget=Budget(0))
    with pytest.raises(BudgetExceeded):
        engine.products([7, 1001], [2, 5])


def test_reciprocal_bounds_against_every_small_integer_factor_pair() -> None:
    from vf_factor.fused import factor_upper_bound, fermat_midpoint_cap

    for n in range(3, 1002, 2):
        root = isqrt(n)
        pairs = [(d, (d + n // d) // 2) for d in range(1, root + 1)
                 if n % d == 0]
        selected_coverage = {0, ceil_root(n) - 1, ceil_root(n),
                             ceil_root(n) + 1, (n + 1) // 2}
        for _, midpoint in pairs:
            selected_coverage.update((midpoint - 1, midpoint, midpoint + 1))
        for covered in selected_coverage:
            upper = factor_upper_bound(n, covered)
            for divisor, midpoint in pairs:
                # Integrality of odd-factor midpoints is essential to the
                # sharper upper bound using the next integer R=A+1.
                assert (divisor <= upper) == (midpoint > covered), (n, covered, divisor)
        for lower in range(1, root + 1):
            cap = fermat_midpoint_cap(n, lower)
            for divisor, midpoint in pairs:
                assert (midpoint <= cap) == (divisor >= lower), (n, lower, divisor)


@pytest.mark.parametrize("p,q", (
    (83, 97), (101, 7907), (19, 531211069),
    (4294967279, 4294967291),
    (18446744073709551557, 18446744073709551629),
))
def test_reciprocal_bounds_include_the_exact_factor_at_integer_endpoints(p: int, q: int) -> None:
    from vf_factor.fused import factor_upper_bound, fermat_midpoint_cap

    n = p * q
    midpoint = (p + q) // 2
    assert factor_upper_bound(n, midpoint - 1) == p
    assert factor_upper_bound(n, midpoint) < p
    assert fermat_midpoint_cap(n, p) == midpoint
    for covered in (ceil_root(n) - 1, ceil_root(n), midpoint - 1, midpoint, midpoint + 1):
        assert (p <= factor_upper_bound(n, covered)) == (midpoint > covered)


def test_python_cursor_has_the_same_exact_coverage_and_resumption() -> None:
    from vf_factor.fused import PythonFermatCursor

    for n in range(3, 502, 2):
        cursor = PythonFermatCursor(n)
        factor, a, b = direct_fermat(n)
        expected_candidates = sum(eligible_midpoint(n, t)
                                  for t in range(ceil_root(n), a + 1))
        for iteration in range(1000):
            state = cursor.advance((1, 3, 7)[iteration % 3])
            validate_certificate(n, state)
            if state["factor"] or state["exhausted"]:
                break
        else:
            pytest.fail("Python cursor did not finish a small input")
        assert state["factor"] == factor
        assert state["last_midpoint"] == a
        assert state["candidates"] == expected_candidates
        if factor:
            assert state["half_gap"] == b
        assert cursor.advance(0)["delta_candidates"] == 0


@NATIVE
@pytest.mark.parametrize("n", (3, 5, 7, 101, 1009, 551, 8051, 101 * 5003))
def test_python_native_cursor_snapshot_parity_including_prime_endpoints(n: int) -> None:
    from vf_factor.fused import PythonFermatCursor

    python_cursor = PythonFermatCursor(n)
    fields = ("n", "factor", "candidates", "square_tests", "delta_candidates",
              "first_midpoint", "last_midpoint", "start_midpoint", "next_midpoint",
              "covered_through", "half_gap", "root_ceiling", "endpoint",
              "residue_count", "period", "status", "exhausted")
    with fusion_native.FermatCursor(n) as native_cursor:
        initial_python, initial_native = python_cursor.snapshot(), native_cursor.snapshot()
        assert {k: initial_python[k] for k in fields} == {k: initial_native[k] for k in fields}
        for iteration in range(1000):
            cap = (1, 3, 7, 19)[iteration % 4]
            actual_python = python_cursor.advance(cap)
            actual_native = native_cursor.advance(cap)
            assert {k: actual_python[k] for k in fields} == {k: actual_native[k] for k in fields}
            if actual_python["factor"] or actual_python["exhausted"]:
                break
        else:
            pytest.fail("Matched cursor scans did not finish a small input")
        if actual_python["exhausted"]:
            assert isprime(n)
            assert actual_python["covered_through"] == (n + 1) // 2
            assert actual_python["next_midpoint"] == (n + 1) // 2 + 1
        terminal_python, terminal_native = python_cursor.advance(0), native_cursor.advance(0)
        assert {k: terminal_python[k] for k in fields} == {k: terminal_native[k] for k in fields}


FUSED_OPTIONS = (
    {"native": False},
    pytest.param({"native": True}, marks=NATIVE),
    {"native": False, "feedback": False},
    {"native": False, "shared_polynomial": False},
    {"native": False, "vf_probes": False},
    {"native": False, "fermat_enabled": False},
)


def assert_fused_correct(result, n: int) -> None:
    assert result.n == n
    if isprime(n):
        assert result.status == "prime", result.to_dict()
        assert result.factor is None
    else:
        assert result.status == "factor", result.to_dict()
        assert 1 < result.factor < n
        assert n % result.factor == 0


@pytest.mark.parametrize("options", FUSED_OPTIONS)
def test_fused_algorithm_and_ablations_on_all_small_odd_inputs(options: dict) -> None:
    from vf_factor.fused import factor

    for n in range(3, 1002, 2):
        result = factor(n, timeout=None, max_steps=100_000,
                        initial_degree=4, max_degree=16, **options)
        assert_fused_correct(result, n)


@pytest.mark.parametrize("options", FUSED_OPTIONS)
def test_fused_semiprimes_beyond_the_shared_scalar_prefix(options: dict) -> None:
    from vf_factor.fused import factor

    rng = random.Random(20261119)
    for _ in range(25):
        p = int(nextprime(rng.randrange(121, 4000)))
        q = int(nextprime(rng.randrange(4000, 12000)))
        n = p * q
        # Constructed factors are solely an oracle; the implementation is
        # passed n alone, with no seed-factor or input-family metadata.
        result = factor(n, timeout=None, max_steps=1_000_000,
                        initial_degree=4, max_degree=32, initial_fermat_tests=0,
                        **options)
        assert_fused_correct(result, n)


@pytest.mark.parametrize("n", (7**3, 11**5, 127**3, 131**4, 131**5, 551, 798607))
@pytest.mark.parametrize("vf_probes", (False, True))
def test_fused_product_path_prime_powers_and_regressions(n: int, vf_probes: bool) -> None:
    from vf_factor.fused import factor

    result = factor(n, timeout=None, fermat_enabled=False, vf_probes=vf_probes,
                    initial_degree=4, max_degree=32, max_steps=1_000_000)
    assert_fused_correct(result, n)


@pytest.mark.parametrize("native", (False, True))
@pytest.mark.parametrize("n", (
    18446744073709551557 * 18446744073709551629,
    101 * int(nextprime(2**70)),
    ((1 << 90) + 7)**2,
))
def test_fused_arbitrary_precision_inputs_and_native_fallback(native: bool, n: int) -> None:
    from vf_factor.fused import factor

    result = factor(n, timeout=None, max_steps=100_000, native=native)
    assert result.status == "factor", result.to_dict()
    assert 1 < result.factor < n and n % result.factor == 0
    if not result.details.get("precheck"):
        assert result.details["fermat_backend"] == "python_bigint"


@pytest.mark.parametrize("cap", (0, 1, 3, 7, 31, 128, 1000))
@pytest.mark.parametrize("options", (
    {"native": False},
    pytest.param({"native": True}, marks=NATIVE),
    {"native": False, "fermat_enabled": False},
))
def test_fused_one_work_cap_covers_every_search_component(cap: int, options: dict) -> None:
    from vf_factor.fused import factor

    result = factor(1000003 * 2000003, timeout=None, max_steps=cap,
                    initial_degree=4, max_degree=16, **options)
    assert result.status == "budget", result.to_dict()
    assert result.factor is None
    counts = result.counters
    total = counts["scalar_tests"] + counts["represented_candidates"] + counts["fermat_tests"]
    assert total == counts["work_units"] <= cap


@pytest.mark.parametrize("native", (False, pytest.param(True, marks=NATIVE)))
def test_fused_zero_global_timeout_prevents_all_nonprecheck_work(native: bool) -> None:
    from vf_factor.fused import factor

    result = factor(1000003 * 2000003, timeout=0, native=native)
    assert result.status == "budget"
    assert result.factor is None
    assert result.counters["work_units"] == 0
    assert result.counters["rounds"] == 0


@pytest.mark.parametrize("n", (2, 3, 5, 10, 49, 2**64 - 1,
                                    18446744073709551557, ((1 << 90) + 7)**2))
def test_fused_shared_prechecks_precede_search_budgets(n: int) -> None:
    from vf_factor.fused import factor

    result = factor(n, timeout=0, max_steps=0)
    assert result.details["precheck"]
    assert result.counters["work_units"] == 0
    assert_fused_correct(result, n)


def test_fused_unsearched_large_prime_is_unresolved_not_certified() -> None:
    from vf_factor.fused import factor

    result = factor((1 << 89) - 1, timeout=None, max_steps=0)
    assert result.status == "budget"
    assert result.factor is None


@pytest.mark.parametrize("n", (-1, 0, 1, True, 77.0, "77"))
def test_fused_rejects_invalid_input(n) -> None:
    from vf_factor.fused import factor

    with pytest.raises(ValueError):
        factor(n)


@pytest.mark.parametrize("options", (
    {"timeout": -1}, {"timeout": True}, {"timeout": float("nan")},
    {"timeout": float("inf")}, {"max_steps": -1}, {"max_steps": True},
    {"max_steps": 1.5}, {"initial_degree": 0}, {"initial_degree": True},
    {"initial_degree": 16, "max_degree": 8}, {"max_degree": 0},
    {"initial_fermat_tests": -1}, {"fermat_ratio": 0},
    {"fermat_ratio": float("nan")}, {"native": "yes"}, {"feedback": 1},
    {"vf_probes": 1}, {"shared_polynomial": 1}, {"fermat_enabled": 1}, {"trace": 1},
))
def test_fused_rejects_invalid_search_controls(options: dict) -> None:
    from vf_factor.fused import factor

    with pytest.raises(ValueError):
        factor(8051, **options)


def test_fused_trace_certifies_both_frontiers_and_joint_polynomial_batches() -> None:
    from vf_factor.fused import factor, factor_upper_bound, fermat_midpoint_cap

    p, q = 1000003, 2000003
    n = p * q
    result = factor(n, timeout=None, max_steps=20_000, native=False,
                    initial_fermat_tests=0, initial_degree=4, max_degree=32,
                    trace=True)
    trace = result.details["trace"]
    assert trace
    assert any(event["event"] == "product_plan" and event["divisor_blocks"]
               and event["vf_blocks"] for event in trace)
    assert any(event["event"] == "product" for event in trace)
    assert any(event["event"] == "fermat" and event["chunk_tests"] for event in trace)
    assert result.counters["mixed_batches"] > 0
    assert result.counters["lower_bound_updates"] > 0
    assert result.counters["upper_bound_updates"] > 0
    previous_low, previous_high, previous_work = 7, isqrt(n), 0
    caps = []
    for event in trace:
        assert previous_low <= event["lo"] <= p
        assert p <= event["hi"] <= previous_high
        assert previous_work <= event["work_units"] <= 20_000
        if event["event"] != "fermat_factor":
            assert event["hi"] == min(isqrt(n), factor_upper_bound(n, event["fermat_covered"]))
        expected_cap = fermat_midpoint_cap(n, event["lo"])
        assert event["fermat_cap"] == expected_cap
        if event["event"] == "fermat":
            assert event["applied_cap"] == expected_cap
            caps.append(expected_cap)
        if event["event"] == "product_plan":
            assert event["strip_end"] <= event["hi"]
        previous_low, previous_high, previous_work = event["lo"], event["hi"], event["work_units"]
    assert caps and min(caps) < fermat_midpoint_cap(n, 7)
    assert result.counters["work_units"] <= 20_000

