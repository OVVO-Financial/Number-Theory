"""Deterministic benchmark inputs, generated independently of timed algorithms.

This module is a test-data generator, not part of any factorization algorithm.
SymPy is used only to select/validate primes and produce reference metadata.
Timed algorithms must receive only ``case['n']``, never ``p`` or ``q``.

The held-out seed and family counts are fixed before benchmarking. The tuning
split uses a different seed. Published and regression inputs intentionally
overlap splits because they are explicit examples, not performance holdouts.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from sympy import isprime, nextprime, prevprime

DEFAULT_SEED = 20261117
TARGET_BITS = (24, 32, 40, 48, 56, 64)


def _prime_between(rng: random.Random, lower: int, upper: int) -> int:
    """Choose a reproducible prime in [lower, upper), without timing feedback."""
    if lower >= upper:
        raise ValueError("empty prime interval")
    candidate = int(nextprime(rng.randrange(max(2, lower), upper)))
    if candidate >= upper:
        candidate = int(prevprime(upper))
    if candidate < lower:
        raise ValueError("prime interval contains no prime")
    return candidate


def _random_pair(rng: random.Random, bits: int, pbits: int) -> tuple[int, int]:
    """Choose p, q with an exactly bits-bit product and specified p size."""
    lower_n, upper_n = 1 << (bits - 1), 1 << bits
    for _ in range(1000):
        p = _prime_between(rng, 1 << (pbits - 1), 1 << pbits)
        qlo = max(p + 1, (lower_n + p - 1) // p)
        qhi = (upper_n - 1) // p + 1
        if qlo >= qhi:
            continue
        q = _prime_between(rng, qlo, qhi)
        if p != q and (p * q).bit_length() == bits:
            return min(p, q), max(p, q)
    raise RuntimeError("could not generate semiprime with requested dimensions")


def _balanced_pair(rng: random.Random, bits: int) -> tuple[int, int]:
    # Require a visibly nonclose ratio; no algorithm run is used to select it.
    for _ in range(1000):
        p, q = _random_pair(rng, bits, bits // 2)
        if 6 * p <= 5 * q <= 16 * p:
            return p, q
    raise RuntimeError("could not generate nonclose balanced pair")


def _close_pair(rng: random.Random, bits: int) -> tuple[int, int]:
    half = bits // 2
    # For even bits, p >= .75*2**half ensures an exactly bits-bit product.
    lower = 3 * (1 << half) // 4
    upper = 7 * (1 << half) // 8
    p = _prime_between(rng, lower, upper)
    gap_start = rng.randrange(2, 48)
    q = int(nextprime(p + gap_start))
    if (p * q).bit_length() != bits:
        raise ArithmeticError("incorrect close-factor target bit count")
    return p, q


def build_corpus(seed: int = DEFAULT_SEED, split: str = "heldout") -> list[dict]:
    """Return 83 heldout cases (75 composites, 8 primes), or a smaller tune set.

    The 75 composites comprise 3 published examples, 4 anchor regressions,
    24 balanced nonclose semiprimes, 18 close semiprimes, 18 unbalanced
    semiprimes, 4 squares, and 4 odd prime powers. Family metadata enables
    separate comparisons rather than averaging incompatible input regimes.
    """
    if split not in ("heldout", "tuning"):
        raise ValueError("split must be 'heldout' or 'tuning'")
    effective_seed = seed if split == "heldout" else seed + 1_000_003
    rng = random.Random(effective_seed)
    cases: list[dict] = []
    used: set[int] = set()

    def add(case_id: str, family: str, p: int | None, q: int | None,
            *, n: int | None = None, target_bits: int | None = None,
            provenance: str = "fixed-seed prime selection", **extra) -> None:
        if n is None:
            assert p is not None and q is not None
            n = p * q
        if n in used:
            raise ArithmeticError(f"duplicate benchmark integer: {n}")
        used.add(n)
        cases.append({
            "case_id": case_id, "family": family, "n": n,
            "p": p, "q": q, "bits": n.bit_length(),
            "target_bits": target_bits or n.bit_length(),
            "split": split, "seed": effective_seed,
            "provenance": provenance,
            "expected_status": "prime" if p is None else "factor",
            **extra,
        })

    for n, p, q in ((8051, 83, 97), (798607, 101, 7907),
                    (978508015703, 752867, 1299709)):
        assert p * q == n and isprime(p) and isprime(q)
        add(f"published-{n}", "published", p, q,
            provenance="Complex Space Factorization / Iterated Averages examples")
    for n, p, q in ((511, 7, 73), (1351, 7, 193),
                    (2291, 29, 79), (8251, 37, 223)):
        add(f"anchor-{n}", "anchor_regression", p, q,
            provenance="fixed seed / true midpoint regression")

    balanced_count = 4 if split == "heldout" else 1
    close_count = 3 if split == "heldout" else 1
    unbalanced_count = 3 if split == "heldout" else 1
    small_factor_bits = {
        24: (8, 9, 10), 32: (8, 10, 12), 40: (8, 12, 16),
        48: (8, 14, 20), 56: (8, 16, 20), 64: (8, 16, 20),
    }
    for bits in TARGET_BITS:
        for i in range(balanced_count):
            p, q = _balanced_pair(rng, bits)
            add(f"balanced-{bits}-{i}", "balanced", p, q,
                target_bits=bits, factor_ratio=q / p)
        for i in range(close_count):
            p, q = _close_pair(rng, bits)
            add(f"close-{bits}-{i}", "close", p, q,
                target_bits=bits, factor_gap=q - p)
        for i, pbits in enumerate(small_factor_bits[bits][:unbalanced_count]):
            p, q = _random_pair(rng, bits, pbits)
            add(f"unbalanced-{bits}-{pbits}", "unbalanced", p, q,
                target_bits=bits, small_factor_bits=p.bit_length())

    for bits in (16, 32, 48, 64):
        half = bits // 2
        p = _prime_between(rng, 3 * (1 << half) // 4, 7 * (1 << half) // 8)
        add(f"square-{bits}", "square", p, p,
            target_bits=bits, exponent=2)
    for p in (7, 101, 1009, 65537):
        add(f"cube-{p}", "prime_power", p, p * p,
            provenance="fixed odd prime cube", exponent=3)
    for bits in (16, 24, 32, 40, 48, 56, 63, 64):
        prime = _prime_between(rng, 1 << (bits - 1), 1 << bits)
        add(f"prime-{bits}", "prime", None, None, n=prime,
            target_bits=bits)
    return cases


def validate_corpus(cases: list[dict]) -> None:
    """Independent metadata validation; never run in a timed factor call."""
    ids, values = set(), set()
    for case in cases:
        n, p, q = case["n"], case["p"], case["q"]
        assert case["case_id"] not in ids
        assert n not in values
        ids.add(case["case_id"])
        values.add(n)
        assert n.bit_length() == case["bits"] == case["target_bits"]
        if p is None:
            assert q is None and isprime(n)
        else:
            assert 1 < p <= q < n and p * q == n and isprime(p)
            if case["family"] == "prime_power":
                assert q == p * p
            else:
                assert isprime(q)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", choices=("heldout", "tuning"), default="heldout")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = build_corpus(args.seed, args.split)
    validate_corpus(cases)
    payload = json.dumps(cases, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
