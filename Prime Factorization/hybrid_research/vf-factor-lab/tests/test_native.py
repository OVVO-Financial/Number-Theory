"""Reproducible correctness checks for the optional native64 baselines.

Build with ``python3 build_native.py`` before running these tests. They skip
cleanly when the native kernel has not been built. No speed threshold is
asserted; this file verifies exactness, source fidelity and guard behavior.
"""
from __future__ import annotations

from math import isqrt
import random

import pytest
from sympy import isprime, nextprime

from vf_factor import baselines, native


pytestmark = pytest.mark.skipif(
    not native.is_available(),
    reason="optional native64 kernel not built; run python3 build_native.py",
)

METHODS = (
    native.published_vf,
    native.pollard_brent,
    native.trial_division,
    native.fermat,
)


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


@pytest.mark.parametrize("function", METHODS, ids=lambda function: function.__name__)
def test_native_all_small_odd_inputs(function) -> None:
    # 5,000 inputs per method, matching the 20,000-check validation record.
    for n in range(3, 10_002, 2):
        assert_correct(function(n, timeout=None), n)


@pytest.mark.parametrize("n,factor,level,distance,argument,calls", (
    (8051, 83, 3, 0, 7055, 1),
    (798607, 101, 4, 35, 748713, 30),
    (978508015703, 752867, 6, 50745, 963218792667, 40598),
    (511, 7, 2, 10, 378, 14),
    (1351, 7, 3, 10, 1176, 13),
    (2291, 79, 4, 15, 2133, 21),
))
def test_native_published_vf_exact_source_order(
    n: int, factor: int, level: int, distance: int, argument: int, calls: int,
) -> None:
    actual = native.published_vf(n, timeout=None)
    reference = baselines.published_vf(n, timeout=None)
    assert actual.factor == reference.factor == factor
    assert actual.details["selected_level"] == level
    assert actual.details["realized_b"] == distance
    assert actual.details["gcd_argument"] == argument
    assert actual.counters["gcd_calls"] == reference.counters["gcd_calls"] == calls
    for key in ("selected_level", "realized_b", "gcd_argument", "direction", "offset"):
        assert actual.details[key] == reference.details[key]


def random_cases_and_roots():
    """Use exactly the seed and generation sequence in the validation record."""
    rng = random.Random(57)
    numbers = []
    for _ in range(100):
        p = int(nextprime(rng.randrange(100, 50_000)))
        q = int(nextprime(rng.randrange(100, 50_000)))
        numbers.append(p * q)
    roots = [0, 1, 2, 3, 4, 5, 2**64 - 1, 2**64 - 2,
             2**63, 2**32, 2**32 - 1]
    roots.extend(rng.getrandbits(64) for _ in range(500))
    for r in (2**32 - 1, 2**31 + 17, 2**16 - 1):
        roots.extend((r * r - 1, r * r, r * r + 1))
    return numbers, roots


@pytest.mark.parametrize("function", METHODS, ids=lambda function: function.__name__)
def test_native_seeded_random_semiprimes(function) -> None:
    numbers, _ = random_cases_and_roots()
    for n in numbers:
        # The generated factors are used only to construct independent test
        # inputs. The algorithm receives n alone, without factor metadata.
        assert_correct(function(n, timeout=None), n)


def test_native_random_vf_search_order_with_hard_cap() -> None:
    numbers, _ = random_cases_and_roots()
    for n in numbers:
        reference = baselines.published_vf(n, timeout=None, max_steps=10_000)
        actual = native.published_vf(n, timeout=None, max_steps=10_000)
        assert actual.status == reference.status
        assert actual.factor == reference.factor
        if reference.details.get("precheck"):
            assert actual.details.get("precheck")
            continue
        assert actual.counters["gcd_calls"] == reference.counters["gcd_calls"]
        assert actual.details["selected_level"] == reference.details["selected_level"]
        if actual.factor:
            for key in ("realized_b", "gcd_argument", "direction", "offset"):
                assert actual.details[key] == reference.details[key]


def test_native_integer_square_root_520_cases() -> None:
    _, roots = random_cases_and_roots()
    assert len(roots) == 520
    kernel = native._load()
    for n in roots:
        assert kernel.vf_native_isqrt(n) == isqrt(n), n


@pytest.mark.parametrize("function", METHODS, ids=lambda function: function.__name__)
@pytest.mark.parametrize("cap", (0, 1, 3, 10))
def test_native_hard_step_limits(function, cap: int) -> None:
    result = function(1000003 * 2000003, timeout=None, max_steps=cap)
    assert result.status == "budget"
    assert result.factor is None
    assert result.counters["candidates"] <= cap
    assert result.details["budget_reason"] == "step budget"


@pytest.mark.parametrize("function", METHODS, ids=lambda function: function.__name__)
def test_native_zero_timeout(function) -> None:
    result = function(1000003 * 2000003, timeout=0)
    assert result.status == "budget"
    assert result.factor is None
    assert result.details["budget_reason"] == "time budget"


@pytest.mark.parametrize("function", METHODS, ids=lambda function: function.__name__)
@pytest.mark.parametrize("n", (-1, 0, 1, True, 77.0, "77", 2**64, 2**64 + 1))
def test_native_invalid_and_out_of_range_inputs(function, n) -> None:
    with pytest.raises(ValueError):
        function(n)


@pytest.mark.parametrize("function", METHODS, ids=lambda function: function.__name__)
@pytest.mark.parametrize("n", (2, 4, 49, 2**64 - 1,
                                (2**32 - 5)**2, 18446744073709551557))
def test_native_shared_prechecks_at_machine_boundaries(function, n: int) -> None:
    # Prechecks intentionally precede the search budget for every backend.
    result = function(n, timeout=0, max_steps=0)
    assert_correct(result, n)
    assert result.details.get("precheck")


def test_native_backend_and_timing_metadata() -> None:
    info = native.backend_info()
    assert info["abi_version"] == 1
    assert info["rho_rng"] == "std::mt19937_64"
    result = native.published_vf(8051, timeout=None)
    assert result.details["backend"] == "native64"
    assert result.algorithm == "published_vf_native64"
    # This asserts accounting, not a performance threshold.
    assert result.elapsed_seconds >= result.details["kernel_seconds"] >= 0
