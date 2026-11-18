"""Exact baseline searches and the published iterated-averages heuristic.

``published_vf`` reproduces the search order in the project's later
``Iterated Averages.md`` Julia code, not an inferred replacement for it.
It uses exact rational remainders, rejects gcd == n, and has explicit budgets.
The common prechecks are applied equally to every algorithm.

``max_steps`` means divisor tests (trial division), a-candidates (Fermat),
polynomial-map evaluations (Brent rho), or GCD candidates (published VF).
Timeouts are checked at coarse checkpoints to avoid a clock call per candidate.
"""
from __future__ import annotations

from math import gcd, isqrt
from random import Random

from .common import (
    Budget,
    BudgetExceeded,
    SearchResult,
    ceil_sqrt,
    checked_result,
    precheck,
)


def _begin(n: int, algorithm: str, timeout: float | None):
    budget = Budget(timeout)
    result = SearchResult(algorithm=algorithm, n=n)
    found = precheck(n)
    if found is not None:
        result.factor, result.status = found
        result.details["precheck"] = True
        result.elapsed_seconds = budget.elapsed
        return budget, result, True
    return budget, result, False


def _end(result: SearchResult, budget: Budget) -> SearchResult:
    result.elapsed_seconds = budget.elapsed
    return checked_result(result)


def _step_limit(count: int, max_steps: int | None) -> None:
    if max_steps is not None and count >= max_steps:
        raise BudgetExceeded("step budget")


def average_fraction(n: int, level: int) -> tuple[int, int]:
    """Return numerator and denominator of A_level, without reducing it."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer >= 2")
    if not isinstance(level, int) or level < 1:
        raise ValueError("level must be a positive integer")
    denominator = 1 << level
    return (denominator - 1) * n + ceil_sqrt(n), denominator


def published_vf_parameters(n: int) -> dict:
    """Return the later Julia VF routine's exact level-selection parameters.

    Julia initializes A_1, appends floor(log2((n-9)/6)) further averages,
    and retains ceil(log2(number_of_averages)) of them. Its highest-remainder
    rule selects the first tie, then raises a selected level 1 to level 2.
    Inputs below 49 are padded to two levels; every composite in that range
    is already handled by the shared prechecks.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer >= 2")
    seed = ceil_sqrt(n)
    max_imaginary_floor = max(0, (n - 9) // 6)
    original_count = max(1, max_imaginary_floor.bit_length())
    original_cut = (original_count - 1).bit_length()
    cut = max(2, original_cut)
    scaled_remainders = [
        ((seed - n) % (1 << k)) << (cut - k)
        for k in range(1, cut + 1)
    ]
    selected_level = max(2, 1 + max(range(cut), key=scaled_remainders.__getitem__))
    numerator, denominator = average_fraction(n, selected_level)
    floor_average, remainder = divmod(numerator, denominator)
    return {
        "seed": seed,
        "max_imaginary_floor": max_imaginary_floor,
        "average_count": original_count,
        "cut": cut,
        "selected_level": selected_level,
        "numerator": numerator,
        "denominator": denominator,
        "floor_average": floor_average,
        "ceil_average": floor_average + (remainder != 0),
        "remainder_numerator": remainder,
        "offsets": (0, 1, 2) if n % 10 in (1, 9) else (0, 1),
        "fractional_remainders": [
            {"level": k, "numerator": (seed - n) % (1 << k), "denominator": 1 << k}
            for k in range(1, cut + 1)
        ],
        "small_input_level_padding": original_cut < 2,
    }


def trial_division(
    n: int, *, timeout: float | None = 2.0, max_steps: int | None = None, **_options
) -> SearchResult:
    """Trial divide by the eight reduced residue classes modulo 30."""
    budget, result, done = _begin(n, "trial_division", timeout)
    counts = result.counters
    counts.update(divisor_tests=0, candidates=0)
    result.details["wheel"] = 30
    if done:
        return result
    limit = isqrt(n)
    divisor, index = 7, 0
    increments = (4, 2, 4, 2, 4, 6, 2, 6)
    tests = 0
    try:
        budget.check()
        while divisor <= limit:
            _step_limit(tests, max_steps)
            if tests & 1023 == 0:
                budget.check()
            tests += 1
            if n % divisor == 0:
                result.factor = divisor
                break
            divisor += increments[index]
            index = (index + 1) & 7
        else:
            result.status = "prime"
            result.details["coverage"] = "all possible prime divisors <= floor(sqrt(n))"
    except BudgetExceeded as exc:
        result.details["budget_reason"] = str(exc)
    counts.update(divisor_tests=tests, candidates=tests)
    return _end(result, budget)


_SQUARE_RESIDUES_64 = frozenset(i * i % 64 for i in range(64))


def fermat(
    n: int, *, timeout: float | None = 2.0, max_steps: int | None = None, **_options
) -> SearchResult:
    """Parity-sieved Fermat search, with a square-residue filter modulo 64."""
    budget, result, done = _begin(n, "fermat", timeout)
    counts = result.counters
    counts.update(candidates=0, square_tests=0, residue_rejections=0)
    if done:
        return result
    a = ceil_sqrt(n)
    required_parity = 1 if n % 4 == 1 else 0
    if a % 2 != required_parity:
        a += 1
    b_squared = a * a - n
    last_a = (n + 1) // 2
    candidates = square_tests = rejected = 0
    try:
        budget.check()
        while a <= last_a:
            _step_limit(candidates, max_steps)
            if candidates & 1023 == 0:
                budget.check()
            candidates += 1
            if (b_squared & 63) in _SQUARE_RESIDUES_64:
                square_tests += 1
                b = isqrt(b_squared)
                if b * b == b_squared:
                    factor = a - b
                    if 1 < factor < n:
                        result.factor = factor
                        result.details.update(midpoint=a, half_gap=b)
                        break
                    if factor == 1:
                        result.status = "prime"
                        result.details["coverage"] = "all Fermat midpoints through the trivial pair"
                        break
            else:
                rejected += 1
            b_squared += 4 * a + 4
            a += 2
        else:
            result.status = "prime"
            result.details["coverage"] = "all Fermat midpoints through the trivial pair"
    except BudgetExceeded as exc:
        result.details["budget_reason"] = str(exc)
    counts.update(candidates=candidates, square_tests=square_tests, residue_rejections=rejected)
    return _end(result, budget)


def pollard_brent(
    n: int, *, timeout: float | None = 2.0, max_steps: int | None = None,
    seed: int = 0, batch_size: int = 128, **_options
) -> SearchResult:
    """Reproducible Brent rho with batched GCDs and degenerate-block recovery."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    budget, result, done = _begin(n, "pollard_brent", timeout)
    counts = result.counters
    counts.update(map_evaluations=0, modular_multiplications=0, gcd_calls=0,
                  restarts=0, degenerate_blocks=0, recovery_steps=0)
    result.details.update(seed=seed, batch_size=batch_size)
    if done:
        return result
    rng = Random(seed)
    maps = products = gcd_calls = restarts = degenerate = recovery = 0

    def advance(value: int, c: int) -> int:
        nonlocal maps
        _step_limit(maps, max_steps)
        if maps & 1023 == 0:
            budget.check()
        maps += 1
        return (value * value + c) % n

    try:
        budget.check()
        while True:
            y = rng.randrange(1, n - 1)
            c = rng.randrange(1, n - 1)
            cycle_length = 1
            g = 1
            while g == 1:
                x = y
                for _ in range(cycle_length):
                    y = advance(y, c)
                consumed = 0
                while consumed < cycle_length and g == 1:
                    saved_y = y
                    product = 1
                    take = min(batch_size, cycle_length - consumed)
                    for _ in range(take):
                        y = advance(y, c)
                        product = product * abs(x - y) % n
                        products += 1
                    g = gcd(product, n)
                    gcd_calls += 1
                    consumed += take
                cycle_length *= 2
            if g == n:
                degenerate += 1
                g = 1
                while g == 1:
                    saved_y = advance(saved_y, c)
                    g = gcd(abs(x - saved_y), n)
                    gcd_calls += 1
                    recovery += 1
            if 1 < g < n:
                result.factor = g
                break
            restarts += 1
    except BudgetExceeded as exc:
        result.details["budget_reason"] = str(exc)
    counts.update(map_evaluations=maps, candidates=maps,
                  modular_multiplications=maps + products, gcd_calls=gcd_calls,
                  restarts=restarts, degenerate_blocks=degenerate,
                  recovery_steps=recovery)
    return _end(result, budget)


def published_vf(
    n: int, *, timeout: float | None = 2.0, max_steps: int | None = None, **_options
) -> SearchResult:
    """Exact-integer reproduction of Julia VF(n, 5), with safe termination.

    The selected level, order of branch tests, offsets, and ascending boundary
    follow the source. gcd == n is rejected rather than returned as a factor.
    The descending offset-zero branch is unbounded in the source; searching
    2*floor(sqrt(n))+2 rounds suffices for completeness after handling 2 and 5.
    The paper's Fermat half-gap bound is NOT imposed on this later routine.
    """
    budget, result, done = _begin(n, "published_vf", timeout)
    counts = result.counters
    counts.update(candidates=0, gcd_calls=0, rounds=0,
                  descending_gcds=0, ascending_gcds=0, trivial_gcd_hits=0)
    result.details.update(source="Iterated Averages.md: VF(n, 5)",
                          exact_arithmetic=True, fallback_used=False)
    if done:
        return result
    params = published_vf_parameters(n)
    result.details.update(params)
    floor_average = params["floor_average"]
    ceil_average = params["ceil_average"]
    offsets = params["offsets"]
    round_limit = 2 * isqrt(n) + 2
    calls = descending = ascending = trivial = rounds = 0
    distance = 0
    try:
        budget.check()
        for round_index in range(round_limit):
            if round_index & 127 == 0:
                budget.check()
            rounds += 1
            distance = 5 * round_index
            base = floor_average - distance
            for offset in offsets:
                _step_limit(calls, max_steps)
                value = base - offset
                g = gcd(value, n)
                calls += 1
                descending += 1
                if 1 < g < n:
                    result.factor = g
                    result.details.update(direction="descending", offset=offset,
                                          realized_b=distance, gcd_argument=value)
                    break
                trivial += g == n
            if result.factor:
                break
            base = ceil_average + distance
            if base < n - 2:
                for offset in offsets:
                    _step_limit(calls, max_steps)
                    value = base + offset
                    g = gcd(value, n)
                    calls += 1
                    ascending += 1
                    if 1 < g < n:
                        result.factor = g
                        result.details.update(direction="ascending", offset=offset,
                                              realized_b=distance, gcd_argument=value)
                        break
                    trivial += g == n
            if result.factor:
                break
        else:
            result.status = "prime"
            result.details["coverage"] = "2*floor(sqrt(n))+2 unrestricted descending step-5 terms"
    except BudgetExceeded as exc:
        result.details["budget_reason"] = str(exc)
    counts.update(candidates=calls, gcd_calls=calls, rounds=rounds,
                  descending_gcds=descending, ascending_gcds=ascending,
                  trivial_gcd_hits=trivial, scan_distance=distance)
    return _end(result, budget)
