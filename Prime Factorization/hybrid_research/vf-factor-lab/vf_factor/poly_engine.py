"""Fast exact searches of arithmetic progressions modulo a composite integer.

The block polynomial F(X)=prod(X + step*i, i=0..m-1) is built by a
balanced product tree. Evaluating F at start + step*m*j searches m
candidate terms per evaluation. For machine-word moduli we use FLINT's
nmod polynomial multiplication and a monic product/remainder tree; for
larger moduli we use FLINT's native fmpz_mod multipoint evaluation.

No inversion of an arbitrary nonzero residue is required. All polynomial
divisors in the remainder tree are monic, so composite moduli are valid.
This is Strassen-style product search, not a new complexity exponent.

Timeouts are cooperative: a native polynomial operation already in flight
cannot be interrupted, so a timeout may overshoot by that operation's cost.
"""
from __future__ import annotations

from ctypes import c_ulong, sizeof
from math import gcd, isqrt
from time import perf_counter

from flint import fmpz_mod_poly_ctx, nmod_poly

from .common import (
    Budget,
    BudgetExceeded,
    SearchResult,
    checked_result,
    precheck,
)


WORD_BITS = 8 * sizeof(c_ulong)


def _inc(counters: dict, key: str, amount=1) -> None:
    counters[key] = counters.get(key, 0) + amount


def _check(budget: Budget | None) -> None:
    if budget is not None:
        budget.check()


def _validate_progression(n: int, step: int, count: int | None = None) -> None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer >= 2")
    if not isinstance(step, int) or isinstance(step, bool) or step <= 0:
        raise ValueError("step must be a positive integer")
    if count is not None and (
        not isinstance(count, int) or isinstance(count, bool) or count < 0
    ):
        raise ValueError("count must be a nonnegative integer")


class ProgressionProductEngine:
    """A reusable block polynomial for one (N, step, block_size) triple.

    ``products(starts)`` returns the products, modulo N, of the progressions
    ``start + step*i`` for ``i=0..block_size-1``. Reuse this object across
    several averaging levels within one factorization run. No global cache
    exists, so measured setup cost is never silently reused across inputs.
    """

    def __init__(
        self,
        n: int,
        step: int,
        block_size: int,
        *,
        budget: Budget | None = None,
        counters: dict | None = None,
        backend: str = "auto",
    ):
        _validate_progression(n, step)
        if not isinstance(block_size, int) or isinstance(block_size, bool) or block_size < 1:
            raise ValueError("block_size must be a positive integer")
        if backend not in ("auto", "nmod", "fmpz"):
            raise ValueError("backend must be auto, nmod, or fmpz")
        if backend == "nmod" and n >= 1 << WORD_BITS:
            raise ValueError("nmod requires a machine-word modulus")
        self.n = n
        self.step = step
        self.block_size = block_size
        self.budget = budget
        self.counters = counters if counters is not None else {}
        self.backend = (
            "nmod_tree" if backend == "nmod" or (backend == "auto" and n < 1 << WORD_BITS)
            else "fmpz_native"
        )
        self.context = fmpz_mod_poly_ctx(n) if self.backend == "fmpz_native" else None
        self.counters["polynomial_backend"] = self.backend
        self.counters["polynomial_degree"] = max(
            self.counters.get("polynomial_degree", 0), block_size
        )
        _inc(self.counters, "polynomials_built")
        started = perf_counter()
        try:
            _check(self.budget)
            level = []
            for i in range(block_size):
                level.append(self._poly([(step * i) % n, 1]))
                if i % 1024 == 0:
                    _check(self.budget)
            while len(level) > 1:
                _check(self.budget)
                next_level = []
                for i in range(0, len(level), 2):
                    if i + 1 < len(level):
                        next_level.append(level[i] * level[i + 1])
                        _inc(self.counters, "polynomial_multiplications")
                    else:
                        next_level.append(level[i])
                level = next_level
            self.polynomial = level[0]
            _check(self.budget)
        finally:
            _inc(self.counters, "polynomial_build_seconds", perf_counter() - started)

    def _poly(self, coeffs):
        if self.context is not None:
            return self.context(coeffs)
        return nmod_poly(coeffs, self.n)

    def _nmod_evaluate(self, points: list[int]) -> list[int]:
        if len(points) == 1:
            return [int(self.polynomial(points[0]))]
        tree = [[self._poly([(-point) % self.n, 1]) for point in points]]
        while len(tree[-1]) > 1:
            _check(self.budget)
            level = tree[-1]
            following = []
            for i in range(0, len(level), 2):
                if i + 1 < len(level):
                    following.append(level[i] * level[i + 1])
                    _inc(self.counters, "polynomial_multiplications")
                else:
                    following.append(level[i])
            tree.append(following)
        # Every divisor below is monic. In particular, no field assumption
        # or primality of N is needed for these polynomial remainders.
        remainders = [self.polynomial % tree[-1][0]]
        _inc(self.counters, "polynomial_remainders")
        for level in range(len(tree) - 2, -1, -1):
            _check(self.budget)
            remainders = [
                remainders[i // 2] % divisor for i, divisor in enumerate(tree[level])
            ]
            _inc(self.counters, "polynomial_remainders", len(remainders))
        return [int(remainder[0]) for remainder in remainders]

    def products(self, starts: list[int]) -> list[int]:
        """Compute all block products with fast multipoint evaluation."""
        if not starts:
            return []
        _check(self.budget)
        points = [int(start) % self.n for start in starts]
        started = perf_counter()
        try:
            if self.backend == "nmod_tree":
                values = self._nmod_evaluate(points)
            else:
                values = [int(value) for value in self.polynomial.multipoint_evaluate(points)]
            _inc(self.counters, "evaluations", len(points))
            _inc(self.counters, "represented_candidates", len(points) * self.block_size)
        finally:
            _inc(self.counters, "polynomial_evaluation_seconds", perf_counter() - started)
        _check(self.budget)
        return values


def progression_products(
    n: int,
    starts: list[int],
    step: int,
    block_size: int,
    *,
    budget: Budget | None = None,
    counters: dict | None = None,
) -> list[int]:
    """Convenience entry point; reuse ProgressionProductEngine for repeats."""
    if not starts:
        _validate_progression(n, step)
        if block_size < 1:
            raise ValueError("block_size must be positive")
        return []
    engine = ProgressionProductEngine(
        n, step, block_size, budget=budget, counters=counters
    )
    return engine.products(starts)


def find_factor_in_block(
    n: int,
    start: int,
    step: int,
    count: int,
    *,
    counters: dict | None = None,
    budget: Budget | None = None,
) -> int | None:
    """Inspect one exceptional block; trivial GCD=N hits are skipped."""
    _validate_progression(n, step, count)
    counters = counters if counters is not None else {}
    value = start
    for i in range(count):
        if i % 256 == 0:
            _check(budget)
        factor = gcd(value, n)
        _inc(counters, "gcd_calls")
        _inc(counters, "rescue_candidates")
        if 1 < factor < n:
            return factor
        value += step
    _check(budget)
    return None


def factor_from_products(
    n: int,
    starts: list[int],
    products: list[int],
    step: int,
    block_size: int,
    *,
    counters: dict | None = None,
    budget: Budget | None = None,
) -> tuple[int | None, dict]:
    """Inspect block GCDs; return a proper factor and its block metadata."""
    if len(starts) != len(products):
        raise ValueError("starts and products must have equal lengths")
    counters = counters if counters is not None else {}
    for index, (start, product) in enumerate(zip(starts, products)):
        if index % 256 == 0:
            _check(budget)
        factor = gcd(product, n)
        _inc(counters, "gcd_calls")
        _inc(counters, "blocks_checked")
        if 1 < factor < n:
            return factor, {"block_index": index, "block_start": start, "source": "block_gcd"}
        if factor == n:
            _inc(counters, "gcd_n_blocks")
            factor = find_factor_in_block(
                n, start, step, block_size, counters=counters, budget=budget
            )
            if factor is not None:
                return factor, {"block_index": index, "block_start": start, "source": "block_rescue"}
    _check(budget)
    return None, {}


def _begin(algorithm, n, step, timeout, max_steps):
    _validate_progression(n, step)
    if max_steps is not None and (
        not isinstance(max_steps, int) or isinstance(max_steps, bool) or max_steps < 0
    ):
        raise ValueError("max_steps must be a nonnegative integer")
    budget = Budget(timeout, max_steps)
    result = SearchResult(algorithm, n)
    early = precheck(n)
    if early is not None:
        result.factor, result.status = early
        result.details["precheck"] = True
        result.elapsed_seconds = budget.elapsed
        return result, budget, True
    step_factor = gcd(step, n)
    _inc(result.counters, "gcd_calls")
    if 1 < step_factor < n:
        result.factor, result.status = step_factor, "factor"
        result.details["precheck"] = "step_gcd"
        result.elapsed_seconds = budget.elapsed
        return result, budget, True
    if step_factor == n:
        raise ValueError("step must not be a multiple of n")
    return result, budget, False


def ap_scalar(
    n: int,
    *,
    start: int = 1,
    step: int = 5,
    count: int | None = None,
    timeout: float | None = 2.0,
    max_steps: int | None = None,
) -> SearchResult:
    """Individual GCD scan; default range guarantees a split if composite.

    ``count`` and ``max_steps`` measure candidate terms, not CPU operations.
    A deliberately shorter finite count returns ``exhausted``, never prime.
    """
    _validate_progression(n, step, count)
    result, budget, done = _begin("ap_scalar", n, step, timeout, max_steps)
    if done:
        return checked_result(result)
    complete_count = 2 * isqrt(n) + 2
    requested = complete_count if count is None else count
    limit = min(requested, max_steps) if max_steps is not None else requested
    result.details.update(start=start, step=step, requested_candidates=requested)
    value = start
    processed = 0
    try:
        budget.check()
        for index in range(limit):
            if index % 256 == 0:
                budget.check()
            factor = gcd(value, n)
            _inc(result.counters, "gcd_calls")
            processed += 1
            if 1 < factor < n:
                result.factor, result.status = factor, "factor"
                result.details["candidate_index"] = index
                break
            value += step
        else:
            result.status = (
                "budget" if limit < requested
                else "prime" if requested >= complete_count
                else "exhausted"
            )
        budget.steps = processed
    except BudgetExceeded as error:
        result.status = "budget"
        result.details["budget_reason"] = str(error)
    finally:
        result.counters["represented_candidates"] = processed
        result.elapsed_seconds = budget.elapsed
    return checked_result(result)


def ap_fast(
    n: int,
    *,
    start: int = 1,
    step: int = 5,
    count: int | None = None,
    timeout: float | None = 2.0,
    max_steps: int | None = None,
    block_size: int | None = None,
    max_degree: int = 8192,
) -> SearchResult:
    """Fast product-tree/multipoint search of the same progression.

    The default block size is ceil(sqrt(candidate_count)), capped by
    max_degree to bound memory. The cap can increase the number of batches;
    it does not remove any requested candidates. Each run includes its own
    polynomial setup. No previous input's polynomial is reused.
    """
    _validate_progression(n, step, count)
    if (
        not isinstance(max_degree, int) or isinstance(max_degree, bool) or max_degree < 1
        or block_size is not None and (
            not isinstance(block_size, int) or isinstance(block_size, bool) or block_size < 1
        )
    ):
        raise ValueError("polynomial degree limits must be positive")
    result, budget, done = _begin("ap_fast", n, step, timeout, max_steps)
    if done:
        return checked_result(result)
    complete_count = 2 * isqrt(n) + 2
    requested = complete_count if count is None else count
    limit = min(requested, max_steps) if max_steps is not None else requested
    result.details.update(start=start, step=step, requested_candidates=requested)
    if limit == 0:
        result.status = "budget" if limit < requested else "exhausted"
        result.elapsed_seconds = budget.elapsed
        return result
    if block_size is None:
        m = isqrt(limit)
        m += m * m < limit
        m = min(max_degree, m)
    else:
        m = min(block_size, max_degree, limit)
    result.details["block_size"] = m
    result.details["max_degree"] = max_degree
    completed = 0
    try:
        engine = ProgressionProductEngine(
            n, step, m, budget=budget, counters=result.counters
        )
        while completed + m <= limit:
            budget.check()
            num_blocks = min(m, (limit - completed) // m)
            starts = [start + step * (completed + m * j) for j in range(num_blocks)]
            products = engine.products(starts)
            budget.check(num_blocks * m)
            factor, info = factor_from_products(
                n, starts, products, step, m, counters=result.counters, budget=budget
            )
            if factor is not None:
                result.factor, result.status = factor, "factor"
                info["batch_candidate_offset"] = completed
                result.details.update(info)
                break
            completed += num_blocks * m
        else:
            # An incomplete final block is scanned exactly; we never include
            # extra terms outside a caller's explicitly requested interval.
            tail = limit - completed
            if tail:
                previous_rescue_count = result.counters.get("rescue_candidates", 0)
                factor = find_factor_in_block(
                    n, start + step * completed, step, tail,
                    counters=result.counters, budget=budget,
                )
                tail_scanned = result.counters.get("rescue_candidates", 0) - previous_rescue_count
                _inc(result.counters, "represented_candidates", tail_scanned)
                budget.check(tail_scanned)
                if factor is not None:
                    result.factor, result.status = factor, "factor"
                    result.details["source"] = "partial_tail"
            if result.factor is None:
                result.status = (
                    "budget" if limit < requested
                    else "prime" if requested >= complete_count
                    else "exhausted"
                )
    except BudgetExceeded as error:
        result.status = "budget"
        result.details["budget_reason"] = str(error)
    finally:
        result.elapsed_seconds = budget.elapsed
    return checked_result(result)
