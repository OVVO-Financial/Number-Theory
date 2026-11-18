"""Shared result, budget, and exact 64-bit primality handling.

Factorization algorithms use no known factors or reference factorizations.
All methods receive identical small-prime/square/prime prechecks.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from math import isqrt
from time import perf_counter


@dataclass
class SearchResult:
    algorithm: str
    n: int
    factor: int | None = None
    status: str = "budget"
    counters: dict = field(default_factory=dict)
    details: dict = field(default_factory=dict)
    elapsed_seconds: float = 0.0

    def to_dict(self) -> dict:
        result = asdict(self)
        result["cofactor"] = self.n // self.factor if self.factor else None
        return result


class BudgetExceeded(Exception):
    pass


class Budget:
    def __init__(self, seconds: float | None = None, max_steps: int | None = None):
        self.started = perf_counter()
        self.deadline = self.started + seconds if seconds is not None else None
        self.max_steps = max_steps
        self.steps = 0

    def check(self, steps: int = 0) -> None:
        self.steps += steps
        if self.max_steps is not None and self.steps > self.max_steps:
            raise BudgetExceeded("step budget")
        if self.deadline is not None and perf_counter() >= self.deadline:
            raise BudgetExceeded("time budget")

    @property
    def elapsed(self) -> float:
        return perf_counter() - self.started


def ceil_sqrt(n: int) -> int:
    r = isqrt(n)
    return r + (r * r != n)


def is_prime64(n: int) -> bool:
    """Deterministic Miller-Rabin for 0 <= n < 2**64; reject larger input.

    The seven-base 64-bit test is used only as a shared precheck, never as
    a source of factors. Large integers are left to the search algorithms.
    """
    if n >= 1 << 64:
        raise ValueError("is_prime64 is restricted to integers below 2**64")
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        a %= n
        if a == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def precheck(n: int) -> tuple[int | None, str] | None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer >= 2")
    for p in (2, 3, 5):
        if n == p:
            return None, "prime"
        if n % p == 0:
            return p, "factor"
    r = isqrt(n)
    if r * r == n:
        return r, "factor"
    if n < 1 << 64 and is_prime64(n):
        return None, "prime"
    return None


def checked_result(result: SearchResult) -> SearchResult:
    if result.factor is not None:
        if not (1 < result.factor < result.n and result.n % result.factor == 0):
            raise ArithmeticError(f"invalid factor produced: {result}")
        result.status = "factor"
    return result
