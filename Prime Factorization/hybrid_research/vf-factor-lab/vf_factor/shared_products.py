"""One polynomial engine for exact GCD searches of different progressions.

For a unit step s modulo N, a block beginning at alpha satisfies

    product(alpha + s*j) = s**m * P_m(alpha / s)  (mod N),
    P_m(X) = product(X + j, j=0..m-1).

The factor s**m is a unit, so it can be omitted when extracting a GCD.
Divisor-wheel and VF witness blocks can therefore share the same polynomial
and one multipoint-evaluation batch. The unshared mode is a matched control
that constructs one original-step polynomial for each distinct step.

These helpers make no search-coverage claims. Only the caller knows which
requests form a complete, contiguous prefix of possible divisors.
"""
from __future__ import annotations

from math import gcd

from .common import Budget
from .poly_engine import ProgressionProductEngine


def _inc(counters: dict, key: str, amount: int = 1) -> None:
    counters[key] = counters.get(key, 0) + amount


def _check(budget: Budget | None) -> None:
    if budget is not None:
        budget.check()


def _integer(value, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    return value


def _modulus(n: int) -> int:
    _integer(n, "n")
    if n < 2:
        raise ValueError("n must be an integer >= 2")
    return n


class NormalizedProductEngine:
    """Reusable polynomial batches for unit-step arithmetic progressions.

    Construction is lazy: an empty request does not build any polynomial.
    The existing polynomial engine owns the evaluations, represented-candidate,
    build-count and polynomial-time counters. This wrapper additionally counts
    nonempty product batches, mixed-step batches, unit checks, and inversions.

    In shared mode the returned residues are GCD-equivalent to the original
    products, not necessarily equal to them. In unshared mode they are the
    original products modulo N. Signed unit steps are supported in both modes.
    """

    def __init__(self, n: int, degree: int, budget: Budget | None = None,
                 counters: dict | None = None, shared: bool = True):
        self.n = _modulus(n)
        _integer(degree, "degree")
        if degree < 1:
            raise ValueError("degree must be a positive integer")
        if not isinstance(shared, bool):
            raise ValueError("shared must be a bool")
        self.degree = degree
        self.budget = budget
        self.counters = counters if counters is not None else {}
        self.shared = shared
        self._engines: dict[int, ProgressionProductEngine] = {}
        self._inverses: dict[int, int] = {1: 1}
        self._unit_steps: set[int] = {1}

    def _step(self, step: int) -> int:
        _integer(step, "step")
        residue = step % self.n
        if residue not in self._unit_steps:
            _check(self.budget)
            factor = gcd(residue, self.n)
            _inc(self.counters, "gcd_calls")
            _inc(self.counters, "step_gcd_calls")
            if factor != 1:
                raise ValueError("every progression step must be coprime to n")
            self._unit_steps.add(residue)
        return residue

    def _engine(self, step: int) -> ProgressionProductEngine:
        engine = self._engines.get(step)
        if engine is None:
            engine = ProgressionProductEngine(
                self.n, step, self.degree,
                budget=self.budget, counters=self.counters,
            )
            self._engines[step] = engine
        return engine

    def products(self, starts, steps) -> list[int]:
        """Return one GCD-equivalent block product per (start, step) pair."""
        starts, steps = list(starts), list(steps)
        if len(starts) != len(steps):
            raise ValueError("starts and steps must have equal lengths")
        if not starts:
            return []
        _check(self.budget)
        # Validate all requests before starting a polynomial construction.
        for start in starts:
            _integer(start, "start")
        canonical_steps = [self._step(step) for step in steps]
        _inc(self.counters, "product_batches")
        if len(set(canonical_steps)) > 1:
            _inc(self.counters, "mixed_batches")

        if self.shared:
            points = []
            for index, (start, step) in enumerate(zip(starts, canonical_steps)):
                if index % 256 == 0:
                    _check(self.budget)
                inverse = self._inverses.get(step)
                if inverse is None:
                    inverse = pow(step, -1, self.n)
                    self._inverses[step] = inverse
                    _inc(self.counters, "step_inversions")
                points.append((start * inverse) % self.n)
            _inc(self.counters, "normalized_blocks", len(points))
            return self._engine(1).products(points)

        groups: dict[int, list[int]] = {}
        for index, step in enumerate(canonical_steps):
            groups.setdefault(step, []).append(index)
        values = [0] * len(starts)
        for step, indices in groups.items():
            _check(self.budget)
            group_values = self._engine(step).products([starts[i] for i in indices])
            for index, value in zip(indices, group_values):
                values[index] = value
        _check(self.budget)
        return values


def recover_mixed_batch(n: int, requests, values, *, budget: Budget | None = None,
                        counters: dict | None = None) -> int | None:
    """Extract a proper factor from a mixed batch, or finish without one.

    Each request supplies its original integer ``start``, ``step``, ``count``
    and a caller-owned ``kind`` label. Values can be actual products or the
    unit-normalized products returned above. The caller is responsible for
    their correspondence. Empty batches perform no GCD.

    The usual path performs exactly one aggregate GCD. If it equals N, block
    GCDs isolate exceptional blocks and scalar recovery skips trivial hits.
    In particular, a VF block containing only multiples of all of N cannot
    conceal a useful factor in a later block. A budget exception leaves the
    batch unfinished; the caller must not certify its divisor coverage.
    """
    n = _modulus(n)
    requests, values = list(requests), list(values)
    if len(requests) != len(values):
        raise ValueError("requests and values must have equal lengths")
    if not requests:
        return None
    counters = counters if counters is not None else {}
    _check(budget)
    for request in requests:
        if not isinstance(request, dict):
            raise ValueError("each request must be a dictionary")
        try:
            start, step, count = (request[key] for key in ("start", "step", "count"))
        except KeyError as exc:
            raise ValueError("requests require start, step, and count") from exc
        _integer(start, "start")
        _integer(step, "step")
        _integer(count, "count")
        if count < 0:
            raise ValueError("count must be nonnegative")
    reduced_values = []
    combined = 1
    for index, value in enumerate(values):
        if index % 256 == 0:
            _check(budget)
        residue = _integer(value, "block value") % n
        reduced_values.append(residue)
        combined = (combined * residue) % n
    _check(budget)
    factor = gcd(combined, n)
    _inc(counters, "gcd_calls")
    _inc(counters, "batch_gcd_calls")
    if 1 < factor < n:
        return factor
    if factor == 1:
        _check(budget)
        return None

    _inc(counters, "gcd_n_batches")
    for request, value in zip(requests, reduced_values):
        _check(budget)
        factor = gcd(value, n)
        _inc(counters, "gcd_calls")
        _inc(counters, "recovery_block_gcd_calls")
        _inc(counters, "blocks_checked")
        if 1 < factor < n:
            return factor
        if factor == n:
            _inc(counters, "gcd_n_blocks")
            candidate = request["start"]
            for index in range(request["count"]):
                if index % 256 == 0:
                    _check(budget)
                factor = gcd(candidate, n)
                _inc(counters, "gcd_calls")
                _inc(counters, "rescue_candidates")
                if 1 < factor < n:
                    return factor
                candidate += request["step"]
    _check(budget)
    return None
