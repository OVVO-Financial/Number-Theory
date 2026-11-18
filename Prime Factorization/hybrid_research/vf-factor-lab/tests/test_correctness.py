"""Independent arithmetic, regression, coverage, and algorithm tests.

No elapsed-time speed threshold is asserted. Timings belong to benchmarks.
The reference oracle is an independent sieve/SymPy and direct integer
products; expected factors are never passed to the search implementations.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt, prod
import random

import pytest
from sympy import isprime

from vf_factor import baselines
from vf_factor.common import Budget, BudgetExceeded, ceil_sqrt, is_prime64
from vf_factor.corpus import build_corpus, validate_corpus
from vf_factor.poly_engine import ap_fast, ap_scalar, progression_products
from vf_factor.vf_algorithm import (
    average_data, average_streams, level_order, product_trial,
    vf_fast, vf_fast_multi, vf_hybrid,
)


METHODS = {
    "trial_division": (baselines.trial_division, {}),
    "fermat": (baselines.fermat, {}),
    "pollard_brent": (baselines.pollard_brent, {"seed": 20261117}),
    "published_vf": (baselines.published_vf, {}),
    "ap_scalar": (ap_scalar, {"start": 7, "step": 5}),
    "ap_fast": (ap_fast, {"start": 7, "step": 5, "block_size": 4, "max_degree": 16}),
    # Disable scalar prefix to exercise the polynomial path even on small cases.
    "vf_fast": (vf_fast, {"scalar_prefix": 0, "initial_degree": 4, "max_degree": 16}),
    "vf_fast_multi": (vf_fast_multi, {"scalar_prefix": 0, "initial_degree": 4, "max_degree": 16}),
    "vf_hybrid": (vf_hybrid, {"scalar_prefix": 0, "initial_degree": 4, "max_degree": 16}),
    "product_trial": (product_trial, {"block_size": 4, "max_degree": 16}),
}


def assert_correct(result, n: int) -> None:
    assert result.n == n
    if isprime(n):
        assert result.status == "prime", result.to_dict()
        assert result.factor is None
    else:
        assert result.status == "factor", result.to_dict()
        assert isinstance(result.factor, int)
        assert 1 < result.factor < n
        assert n % result.factor == 0


@pytest.mark.parametrize("method", tuple(METHODS))
@pytest.mark.parametrize("n", range(3, 1002, 2))
def test_all_small_odd_inputs(method: str, n: int) -> None:
    function, options = METHODS[method]
    result = function(n, timeout=2.0, max_steps=250_000, **options)
    assert_correct(result, n)


@pytest.mark.parametrize("method", tuple(METHODS))
@pytest.mark.parametrize("n", (2, 4, 6, 10, 16, 26, 1024, 2 * 1000000007))
def test_even_inputs(method: str, n: int) -> None:
    function, options = METHODS[method]
    assert_correct(function(n, timeout=1.0, **options), n)


@pytest.mark.parametrize("method", tuple(METHODS))
@pytest.mark.parametrize("n", (0, 1, -77, True, 77.0, "77"))
def test_invalid_input(method: str, n) -> None:
    function, options = METHODS[method]
    with pytest.raises(ValueError):
        function(n, **options)


def test_prime64_against_independent_sieve() -> None:
    limit = 10_000
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            for multiple in range(p * p, limit + 1, p):
                sieve[multiple] = 0
    for n in range(limit + 1):
        assert is_prime64(n) == bool(sieve[n]), n


@pytest.mark.parametrize("n", (
    341, 561, 1105, 1729, 2465, 2821, 6601, 3215031751,
    341550071728321, 3825123056546413051, (1 << 64) - 1,
))
def test_prime64_rejects_pseudoprimes(n: int) -> None:
    assert not is_prime64(n)


def test_prime64_boundary() -> None:
    assert is_prime64(18446744073709551557)
    with pytest.raises(ValueError):
        is_prime64(1 << 64)


@pytest.mark.parametrize("method", tuple(METHODS))
def test_exact_large_square(method: str) -> None:
    # Beyond both binary64 exact integers and machine-word arithmetic.
    root = (1 << 100) + 12345
    n = root * root
    function, options = METHODS[method]
    result = function(n, timeout=1.0, **options)
    assert result.status == "factor"
    assert result.factor == root


@pytest.mark.parametrize("method", tuple(METHODS))
@pytest.mark.parametrize("n", (7**3, 11**3, 7**5, 101**3))
def test_odd_prime_powers(method: str, n: int) -> None:
    function, options = METHODS[method]
    assert_correct(function(n, timeout=2.0, **options), n)


@pytest.mark.parametrize("n", (8051, 798607))
@pytest.mark.parametrize("method", tuple(METHODS))
def test_published_small_examples(method: str, n: int) -> None:
    function, options = METHODS[method]
    assert_correct(function(n, timeout=2.0, **options), n)


@pytest.mark.parametrize("n,factor,level,distance,argument", (
    (798607, 101, 4, 35, 748713),
    (978508015703, 752867, 6, 50745, 963218792667),
))
def test_published_julia_output_reproduction(
    n: int, factor: int, level: int, distance: int, argument: int,
) -> None:
    result = baselines.published_vf(n, timeout=2.0)
    assert result.factor == factor
    assert result.details["selected_level"] == level
    assert result.details["realized_b"] == distance
    assert result.details["gcd_argument"] == argument


@pytest.mark.parametrize("modulus", (8, 15, 77, 143, 8051, (1 << 64) + 1, (1 << 80) + 7))
@pytest.mark.parametrize("degree", (1, 2, 3, 7, 16, 31, 64))
def test_fast_products_against_direct_integer_products(modulus: int, degree: int) -> None:
    rng = random.Random(20261117 + modulus + degree)
    starts = [0, -modulus, modulus, -1]
    starts += [rng.randrange(-2 * modulus, 2 * modulus + 1) for _ in range(7)]
    # Search APIs deliberately require a positive step. Negative starting
    # values and nonunit steps still exercise the composite-ring arithmetic.
    step = rng.choice((1, 2, 5, 7, 10, 15))
    actual = progression_products(modulus, starts, step, degree)
    expected = [prod(start + step * i for i in range(degree)) % modulus
                for start in starts]
    assert [int(value) for value in actual] == expected


def test_product_values_from_worked_example() -> None:
    starts = [997, 1017, 1037, 1057]
    products = progression_products(8051, starts, 5, 4)
    assert [int(p) for p in products] == [5123, 6424, 2870, 291]
    assert [gcd(int(p), 8051) for p in products] == [1, 1, 1, 97]


@pytest.mark.parametrize("degree", (1, 2, 3, 4, 7, 16, 31))
def test_gcd_equal_n_does_not_report_a_factor_or_stop_early(degree: int) -> None:
    result = ap_fast(77, start=77, step=5, count=16,
                     block_size=degree, timeout=2.0)
    assert result.factor in (7, 11)
    assert result.status == "factor"


def test_trivial_only_block_is_not_reported_as_prime_or_factor() -> None:
    for function in (ap_scalar, ap_fast):
        result = function(77, start=77, step=5, count=1, timeout=1.0)
        assert result.factor is None
        assert result.status not in ("factor", "prime")


def test_two_prime_factor_block_is_recovered() -> None:
    # The product contains both 7 and11, so a plain block GCD is77.
    result = ap_fast(77, start=7, step=1, count=8, block_size=8, timeout=1.0)
    assert result.factor in (7, 11)


def test_non_coprime_progression_step_exposes_factor() -> None:
    for function in (ap_scalar, ap_fast):
        result = function(77, start=1, step=7, count=10, timeout=1.0)
        assert result.factor == 7


@pytest.mark.parametrize("n", (511, 1351, 2291, 8251, 8051))
def test_average_formula_against_recursive_fraction_arithmetic(n: int) -> None:
    average = Fraction(isqrt(n) + (isqrt(n) ** 2 != n), 1)
    for k in range(1, 25):
        average = (average + n) / 2
        numerator, denominator = baselines.average_fraction(n, k)
        info = average_data(n, k)
        assert Fraction(numerator, denominator) == average
        assert info["numerator"] == numerator
        assert info["denominator"] == denominator
        assert info["floor"] == average.numerator // average.denominator


def test_average_formula_is_exact_for_huge_inputs() -> None:
    n = 10**80 + 15601
    average = Fraction(ceil_sqrt(n), 1)
    for k in range(1, 121):
        average = (average + n) / 2
        assert Fraction(*baselines.average_fraction(n, k)) == average


def test_original_table_residues() -> None:
    n = 8051
    for k, expected in ((3, [7, 15, 23]), (4, [7, 23, 55])):
        numerator, denominator = baselines.average_fraction(n, k)
        step = 5 * denominator
        surviving = [b for b in range(step)
                     if b % 10 in (3, 5, 7)
                     and (numerator - b) % denominator == 0]
        assert surviving == expected


@pytest.mark.parametrize("k,minus_count,plus_count,success", (
    (1, 402, 402, True), (2, 201, 201, True), (3, 101, 100, True),
    (4, 51, 50, True), (5, 25, 26, True), (6, 12, 13, False),
    (7, 6, 7, False), (8, 3, 4, False), (9, 2, 2, False),
    (10, 1, 1, False),
))
def test_original_finite_bound_coverage_diagnostic(
    k: int, minus_count: int, plus_count: int, success: bool,
) -> None:
    n = 8051
    numerator, denominator = baselines.average_fraction(n, k)
    B = (n - 9) // 6
    minus = [(numerator - b) // denominator for b in range(B + 1)
             if b % 10 in (3, 5, 7) and (numerator - b) % denominator == 0]
    plus = [(numerator + b) // denominator for b in range(B + 1)
            if b % 10 in (3, 5, 7) and (numerator + b) % denominator == 0]
    assert len(minus) == minus_count
    assert len(plus) == plus_count
    assert any(1 < gcd(value, n) < n for value in minus + plus) is success


@pytest.mark.parametrize("n", (511, 1351, 8251))
def test_literal_fermat_digit_transfer_fails_but_new_search_succeeds(n: int) -> None:
    # This records the old sieve's known failure; it is NOT a new search filter.
    assert n % 20 == 11
    for k in range(1, 12):
        numerator, denominator = baselines.average_fraction(n, k)
        assert numerator % 2 == 0
        for ending in (3, 5, 7):
            assert (numerator - ending) % 2 == 1
            assert (numerator + ending) % 2 == 1
    assert_correct(vf_fast(n, scalar_prefix=0, initial_degree=4,
                           max_degree=16, timeout=2.0), n)


@pytest.mark.parametrize("n", (8051, 1351, 2291, 798607, 978508015703))
def test_average_streams_represent_the_published_gcd_arguments(n: int) -> None:
    levels = level_order(n, "multilevel")
    streams = average_streams(n, levels, 32)
    for stream in streams:
        average = Fraction(*baselines.average_fraction(n, stream["level"]))
        floor = average.numerator // average.denominator
        ceiling = -((-average.numerator) // average.denominator)
        for j in range(min(32, stream["count"])):
            represented = stream["start"] + 5 * j
            if stream["direction"] == "descending":
                original = floor - 5 * j - stream["offset"]
            else:
                original = ceiling + 5 * j + stream["offset"]
            assert gcd(represented, n) == gcd(original, n)


@pytest.mark.parametrize("method", tuple(METHODS))
def test_zero_timeout_does_not_claim_prime_on_unsearched_composite(method: str) -> None:
    n = 1000000007 * 1000000009
    function, options = METHODS[method]
    result = function(n, timeout=0.0, **options)
    assert result.factor is None
    assert result.status == "budget"


@pytest.mark.parametrize("method", tuple(METHODS))
def test_zero_work_budget_prevents_nonprecheck_search(method: str) -> None:
    function, options = METHODS[method]
    result = function(77, timeout=1.0, max_steps=0, **options)
    assert result.factor is None
    assert result.status == "budget"


def test_budget_rejects_extra_work() -> None:
    budget = Budget(max_steps=2)
    budget.check(2)
    with pytest.raises(BudgetExceeded):
        budget.check(1)
    with pytest.raises(BudgetExceeded):
        Budget(seconds=0).check()


def test_corpus_is_reproducible_and_reference_validated() -> None:
    heldout = build_corpus()
    assert heldout == build_corpus()
    validate_corpus(heldout)
    assert len(heldout) == 83
    assert sum(case["expected_status"] == "factor" for case in heldout) == 75
    tuning = build_corpus(split="tuning")
    validate_corpus(tuning)
    random_families = {"balanced", "close", "unbalanced"}
    h = {c["n"] for c in heldout if c["family"] in random_families}
    t = {c["n"] for c in tuning if c["family"] in random_families}
    assert h.isdisjoint(t)
    for family in random_families:
        assert sum(c["family"] == family for c in heldout) >= 8


def test_corpus_rejects_unknown_split() -> None:
    with pytest.raises(ValueError):
        build_corpus(split="optimized-after-seeing-results")
