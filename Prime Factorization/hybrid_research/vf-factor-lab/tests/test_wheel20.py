"""Independent checks for the complete Fermat wheel-20 extension.

The oracle enumerates integer midpoint/half-gap pairs directly. It does not
copy the native cyclic-gap construction. Native searches are optional; run
``python3 build_wheel20.py`` to enable them.
"""
from __future__ import annotations

from math import isqrt
import random

import pytest
from sympy import isprime, nextprime

from vf_factor import wheel20


METHODS = (
    wheel20.native_fermat_parity,
    wheel20.native_fermat_wheel20,
    wheel20.native_fermat_parity_qr64,
    wheel20.native_fermat_wheel20_qr64,
)
NATIVE = pytest.mark.skipif(
    not wheel20.is_available(),
    reason="optional wheel20 kernel not built; run python3 build_wheel20.py",
)

# These eight rows are the two separate parity lanes in the verified PDF.
# The digit-only table is their union, not a single four-pair row.
EXPECTED_PAIRS = {
    1: {(1, 0), (5, 2), (5, 8), (9, 0)},
    3: {(2, 1), (2, 9), (8, 1), (8, 9)},
    7: {(4, 3), (4, 7), (6, 3), (6, 7)},
    9: {(3, 0), (5, 4), (5, 6), (7, 0)},
    11: {(0, 3), (0, 7), (4, 5), (6, 5)},
    13: {(3, 4), (3, 6), (7, 4), (7, 6)},
    17: {(1, 2), (1, 8), (9, 2), (9, 8)},
    19: {(0, 1), (0, 9), (2, 5), (8, 5)},
}
QR64 = {r * r % 64 for r in range(64)}


def independent_pairs(n: int) -> set[tuple[int, int]]:
    return {(a, b) for a in range(10) for b in range(10)
            if (a * a - b * b - n) % 20 == 0}


@pytest.mark.parametrize("residue", range(1, 20, 2))
def test_all_odd_mod20_lanes_against_direct_enumeration(residue: int) -> None:
    actual = set(map(tuple, wheel20.admissible_pairs(residue)))
    assert actual == independent_pairs(residue)
    assert len(actual) == (9 if residue % 5 == 0 else 4)
    if residue in EXPECTED_PAIRS:
        assert actual == EXPECTED_PAIRS[residue]
    assert set(wheel20.midpoint_residues(residue)) == {a for a, _ in actual}


@pytest.mark.parametrize("digit", (1, 3, 5, 7, 9))
def test_digit_only_table_requires_both_parity_lanes(digit: int) -> None:
    first = set(map(tuple, wheel20.admissible_pairs(digit)))
    second = set(map(tuple, wheel20.admissible_pairs(digit + 10)))
    assert first.isdisjoint(second)
    expected = {(a, b) for a in range(10) for b in range(10)
                if (a * a - b * b - digit) % 10 == 0
                and (a - b) % 2 == 1}
    assert first | second == expected
    assert len(expected) == (18 if digit == 5 else 8)


def test_pair_classes_are_complete_for_lifted_integer_coordinates() -> None:
    tables = {r: set(map(tuple, wheel20.admissible_pairs(r)))
              for r in range(1, 20, 2)}
    for a in range(1, 151):
        for b in range(a):
            n = a * a - b * b
            if n % 2:
                assert (a % 10, b % 10) in tables[n % 20]
                assert (2 * a * b) % 4 == 0


@pytest.mark.parametrize("residue", EXPECTED_PAIRS)
def test_midpoint_projection_density_and_height_residues(residue: int) -> None:
    pairs = independent_pairs(residue)
    midpoint_count = len({a for a, _ in pairs})
    assert midpoint_count == (3 if residue in (1, 9, 11, 19) else 2)
    # Coupled height information is retained; it follows from these same
    # pairs rather than constituting an independent extra filter.
    heights = {4 * a * b % 40 for a, b in pairs}
    assert heights == ({0} if residue in (1, 9, 11, 19) else {8, 32})


def first_fermat_solution(n: int) -> tuple[int, int, int]:
    """Slow independent oracle over every integer midpoint, without sieving."""
    a = isqrt(n)
    if a * a < n:
        a += 1
    while a <= (n + 1) // 2:
        difference = a * a - n
        b = isqrt(difference)
        if b * b == difference and 1 < a - b < n:
            return a, b, a - b
        a += 1
    raise AssertionError(f"oracle found no proper factor for {n}")


def assert_search_matches_oracle(function, n: int) -> None:
    result = function(n, timeout=None)
    assert not result.details.get("precheck"), result.to_dict()
    a, b, factor = first_fermat_solution(n)
    assert result.status == "factor"
    assert result.factor == factor
    assert result.details["midpoint"] == a
    assert result.details["half_gap"] == b
    assert result.details["last_midpoint"] == a
    start = isqrt(n) + (isqrt(n) ** 2 != n)
    use_wheel = "wheel20" in function.__name__
    use_qr64 = "qr64" in function.__name__
    if use_wheel:
        residues = {x for x, _ in independent_pairs(n)}
        candidates = [r for r in range(start, a + 1) if r % 10 in residues]
    else:
        parity = ((n + 1) // 2) % 2
        candidates = [r for r in range(start, a + 1) if r % 2 == parity]
    square_tests = sum(not use_qr64 or (r * r - n) % 64 in QR64
                       for r in candidates)
    assert result.counters["candidates"] == len(candidates)
    assert result.counters["iterations"] == len(candidates)
    assert result.counters["square_tests"] == square_tests
    assert result.counters["residue_rejections"] == len(candidates) - square_tests
    assert result.details["first_midpoint"] == candidates[0]
    assert result.counters["scan_distance"] == a - start
    assert result.counters["midpoint_span"] == a - candidates[0] + 1


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
def test_551_regression_keeps_the_missing_lane(function) -> None:
    assert (4, 5) in set(map(tuple, wheel20.admissible_pairs(551 % 20)))
    result = function(551, timeout=None)
    assert result.factor == 19
    assert result.details["midpoint"] == 24
    assert result.details["half_gap"] == 5
    assert_search_matches_oracle(function, 551)


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
def test_all_small_odd_inputs(function) -> None:
    for n in range(3, 10_002, 2):
        result = function(n, timeout=None)
        assert result.n == n
        if isprime(n):
            assert result.status == "prime", result.to_dict()
            assert result.factor is None
        else:
            assert result.status == "factor", result.to_dict()
            assert 1 < result.factor < n
            assert n % result.factor == 0


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
def test_random_semiprime_order_and_exact_candidate_counts(function) -> None:
    rng = random.Random(20261117)
    for _ in range(60):
        p = int(nextprime(rng.randrange(101, 5000)))
        q = int(nextprime(rng.randrange(5001, 15000)))
        assert_search_matches_oracle(function, p * q)


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
@pytest.mark.parametrize("n", (
    4294967291 * 4294967279,
    (2**32 + 1) ** 2 - 1000014**2,
))
def test_exact_square_arithmetic_near_uint64_limit(function, n: int) -> None:
    assert 0 < n < 2**64
    assert_search_matches_oracle(function, n)
    if n == (2**32 + 1) ** 2 - 1000014**2:
        result = function(n, timeout=None)
        assert result.details["midpoint"] ** 2 > 2**64


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
@pytest.mark.parametrize("cap", (0, 1, 3, 10))
def test_hard_candidate_budget(function, cap: int) -> None:
    result = function(1000003 * 2000003, timeout=None, max_steps=cap)
    assert result.status == "budget"
    assert result.factor is None
    assert result.counters["candidates"] == cap
    assert result.counters["square_tests"] <= cap
    if cap == 0:
        assert result.details["last_midpoint"] is None


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
def test_zero_timeout_prevents_search(function) -> None:
    result = function(1000003 * 2000003, timeout=0)
    assert result.status == "budget"
    assert result.factor is None
    assert result.counters["candidates"] == 0


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
@pytest.mark.parametrize("n", (-1, 0, 1, True, 77.0, "77", 2**64, 2**64 + 1))
def test_native_range_and_input_validation(function, n) -> None:
    with pytest.raises(ValueError):
        function(n)


@NATIVE
@pytest.mark.parametrize("function", METHODS, ids=lambda f: f.__name__)
@pytest.mark.parametrize("n", (2, 3, 5, 10, 25, 49, 2**64 - 1,
                                18446744073709551557))
def test_prechecks_precede_search_budgets(function, n: int) -> None:
    result = function(n, timeout=0, max_steps=0)
    assert result.details.get("precheck")
    if isprime(n):
        assert result.status == "prime"
    else:
        assert result.status == "factor"
        assert 1 < result.factor < n and n % result.factor == 0
