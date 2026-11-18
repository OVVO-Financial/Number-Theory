"""Exact iterated-average schedules searched with polynomial block products.

No actual factors enter level selection.  The published decimal-neighbour
schedule is treated as a choice of GCD progressions, not a theorem restricting
the true Fermat half-gap.  An unrestricted descending progression gives the
complete 2*sqrt(n) index bound; an exhausted time budget is never called prime.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt

from .common import Budget, BudgetExceeded, SearchResult, ceil_sqrt, checked_result, precheck


def average_data(n: int, level: int) -> dict:
    if level < 1:
        raise ValueError("level must be >= 1")
    denominator = 1 << level
    numerator = (denominator - 1) * n + ceil_sqrt(n)
    floor, remainder = divmod(numerator, denominator)
    return dict(level=level, numerator=numerator, denominator=denominator,
                floor=floor, ceil=floor + bool(remainder), remainder=remainder)


def level_order(n: int, policy: str = "published") -> list[int]:
    from .baselines import published_vf_parameters
    params = published_vf_parameters(n)
    selected = params["selected_level"]
    if policy == "published":
        return [selected]
    if policy == "multilevel":
        candidates = list(range(1, params["cut"] + 1))
        candidates.sort(key=lambda k: (-Fraction(average_data(n, k)["remainder"], 1 << k), k))
        return [selected] + [k for k in candidates if k != selected]
    if policy == "first":
        return [1]
    raise ValueError("unknown level policy")


def average_streams(n: int, levels: list[int], total_indices: int) -> list[dict]:
    offsets = (0, 1, 2) if n % 10 in (1, 9) else (0, 1)
    streams = []
    for level in levels:
        avg = average_data(n, level)
        # Descending x becomes increasing y=n-x with the same gcd.
        for offset in offsets:
            streams.append(dict(start=n - avg["floor"] + offset, level=level,
                                direction="descending", offset=offset, count=total_indices))
        # Preserve the published condition on the ascending BASE before offsets.
        asc_count = max(0, (n - 3 - avg["ceil"]) // 5 + 1)
        for offset in offsets:
            streams.append(dict(start=avg["ceil"] + offset, level=level,
                                direction="ascending", offset=offset,
                                count=min(total_indices, asc_count)))
    return streams


def vf_fast(n: int, *, timeout: float | None = 2.0, max_steps: int | None = None,
            policy: str = "published", scalar_prefix: int = 32,
            initial_degree: int = 32, max_degree: int = 4096,
            hybrid_fermat_steps: int = 0) -> SearchResult:
    """Find one proper factor; use all levels only with policy='multilevel'.

    The optional short Fermat pass is reported explicitly and is disabled in
    the principal ablation, so gains are not falsely attributed to averages.
    """
    from .poly_engine import ProgressionProductEngine, find_factor_in_block
    name = "vf_fast_multi" if policy == "multilevel" else "vf_fast"
    if hybrid_fermat_steps:
        name = "vf_hybrid"
    budget = Budget(timeout)
    counters = dict(gcd_calls=0, scalar_candidates=0, represented_candidates=0,
                    polynomial_builds=0, evaluations=0, fermat_tests=0)
    result = SearchResult(name, n, counters=counters)

    def reserve(terms):
        # A block is evaluated only if its entire represented range fits.
        # Recovery GCDs revisit terms already represented by that block.
        used = counters["represented_candidates"] + counters["fermat_tests"]
        if max_steps is not None and used + terms > max_steps:
            raise BudgetExceeded("represented-candidate budget")

    def finish(factor=None, status="factor", **details):
        result.factor, result.status = factor, status
        result.details.update(details)
        result.elapsed_seconds = budget.elapsed
        return checked_result(result)

    early = precheck(n)
    if early is not None:
        factor, status = early
        return finish(factor, status, precheck=True)
    if initial_degree < 1 or max_degree < initial_degree or scalar_prefix < 0:
        raise ValueError("invalid block configuration")
    try:
        budget.check()
        if hybrid_fermat_steps:
            a = ceil_sqrt(n)
            parity = 0 if n % 4 == 3 else 1
            a += (a - parity) % 2
            for test_index in range(hybrid_fermat_steps):
                if test_index % 256 == 0:
                    budget.check()
                reserve(1)
                counters["fermat_tests"] += 1
                d = a*a-n
                b = isqrt(d)
                if b*b == d and 1 < a-b < n:
                    return finish(a-b, phase="fermat")
                a += 2
            budget.check()

        levels = level_order(n, policy)
        total_indices = 2 * isqrt(n) + 2
        streams = average_streams(n, levels, total_indices)
        result.details.update(levels=levels, policy=policy, streams=len(streams),
                              complete_index_bound=total_indices,
                              scalar_prefix=scalar_prefix, max_degree=max_degree)
        prefix = min(scalar_prefix, total_indices)
        for j in range(prefix):
            for stream in streams:
                if j >= stream["count"]:
                    continue
                reserve(1)
                counters["gcd_calls"] += 1
                counters["scalar_candidates"] += 1
                counters["represented_candidates"] += 1
                g = gcd(stream["start"] + 5*j, n)
                if 1 < g < n:
                    return finish(g, phase="scalar_prefix", index=j, stream=stream)
            budget.check()

        cursor, degree = prefix, initial_degree
        engines = {}
        while cursor < total_indices:
            budget.check()
            width = min(degree * degree, total_indices - cursor)
            full_blocks = width // degree
            requests, starts = [], []
            # Interleave by block number, then by stream, to retain the same
            # broad search order as the scalar routine.
            for block in range(full_blocks):
                j = cursor + block * degree
                for stream in streams:
                    if j + degree <= stream["count"]:
                        starts.append(stream["start"] + 5*j)
                        requests.append((stream, j))
            if requests:
                reserve(len(starts) * degree)
                if degree not in engines:
                    engines[degree] = ProgressionProductEngine(
                        n, 5, degree, budget=budget, counters=counters)
                    counters["polynomial_builds"] += 1
                values = engines[degree].products(starts)
                for check_index, (value, (stream, j), start) in enumerate(zip(values, requests, starts)):
                    if check_index % 256 == 0:
                        budget.check()
                    counters["gcd_calls"] += 1
                    g = gcd(int(value), n)
                    if 1 < g < n:
                        return finish(g, phase="polynomial", degree=degree,
                                      index_range=[j, j+degree-1], stream=stream)
                    if g == n:
                        g = find_factor_in_block(n, start, 5, degree,
                                                 counters=counters, budget=budget)
                        if g is not None:
                            return finish(g, phase="block_recovery", degree=degree,
                                          index_range=[j, j+degree-1], stream=stream)
                budget.check()
            # Handle tails exactly: no hidden widening of the requested range.
            for stream in streams:
                stop = min(cursor + width, stream["count"])
                available = max(0, stop - cursor)
                tail_start = cursor + (available // degree) * degree
                for j in range(tail_start, stop):
                    reserve(1)
                    counters["gcd_calls"] += 1
                    counters["scalar_candidates"] += 1
                    counters["represented_candidates"] += 1
                    g = gcd(stream["start"] + 5*j, n)
                    if 1 < g < n:
                        return finish(g, phase="tail", index=j, stream=stream)
                budget.check()
            cursor += width
            degree = min(2 * degree, max_degree)
        return finish(None, "prime", phase="complete_coverage")
    except BudgetExceeded as exc:
        return finish(None, "budget", reason=str(exc))


def vf_fast_multi(n: int, **kwargs) -> SearchResult:
    return vf_fast(n, policy="multilevel", **kwargs)


def vf_first(n: int, **kwargs) -> SearchResult:
    """Ablation: identical engine and streams, but always average level one."""
    result = vf_fast(n, policy="first", **kwargs)
    result.algorithm = "vf_first"
    return result


def vf_hybrid(n: int, **kwargs) -> SearchResult:
    return vf_fast(n, hybrid_fermat_steps=64, **kwargs)


def product_trial(n: int, *, timeout: float | None = 2.0, max_steps=None,
                  max_degree: int = 4096, **kwargs) -> SearchResult:
    """Fast product baseline over ordinary odd trial divisors, no averages."""
    from .poly_engine import ap_fast
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer >= 2")
    count = max(0, (isqrt(n) - 7) // 2 + 1)
    result = ap_fast(n, start=7, step=2, count=count, timeout=timeout,
                     max_steps=max_steps, max_degree=max_degree, **kwargs)
    result.algorithm = "product_trial"
    result.details["coverage"] = "all odd trial divisors 7..floor(sqrt(n)); common prechecks handle 2,3,5"
    if result.status == "exhausted":
        result.status = "prime"
    return result
