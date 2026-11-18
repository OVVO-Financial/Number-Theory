from __future__ import annotations
import argparse
import json

from .registry import get_algorithms


def main():
    algorithms = get_algorithms()
    parser = argparse.ArgumentParser(description="Exact iterated-average factorization laboratory")
    parser.add_argument("n", type=int, help="integer >= 2")
    parser.add_argument("--method", choices=algorithms, default="vf_hybrid")
    parser.add_argument("--timeout", type=float, default=10.0,
                        help="seconds per factor search; timeout is cooperative")
    parser.add_argument("--max-degree", type=int, default=4096)
    parser.add_argument("--full", action="store_true", help="recursively split composite cofactors")
    args = parser.parse_args()
    function = algorithms[args.method]
    if not args.full:
        result = function(args.n, timeout=args.timeout, max_degree=args.max_degree)
        print(json.dumps(result.to_dict(), indent=2))
        return
    pending, primes, unresolved, runs = [args.n], [], [], []
    while pending:
        n = pending.pop()
        result = function(n, timeout=args.timeout, max_degree=args.max_degree)
        runs.append(result.to_dict())
        if result.factor:
            pending.extend((result.factor, n // result.factor))
        elif result.status == "prime":
            primes.append(n)
        else:
            unresolved.append(n)
    print(json.dumps(dict(input=args.n, prime_factors=sorted(primes),
                          unresolved=sorted(unresolved), runs=runs), indent=2))


if __name__ == "__main__":
    main()
