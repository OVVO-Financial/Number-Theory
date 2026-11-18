"""Command-line entry point for the complete Fermat wheel and its controls."""
import argparse
import json

from vf_factor import wheel20


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("--timeout", type=float, default=10.0)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--qr64", action="store_true", help="add the modulo-64 square filter")
    parser.add_argument("--parity-only", action="store_true", help="use the matched control without wheel 20")
    args = parser.parse_args()
    name = "native_fermat_" + ("parity" if args.parity_only else "wheel20")
    if args.qr64:
        name += "_qr64"
    result = getattr(wheel20, name)(args.n, timeout=args.timeout, max_steps=args.max_steps)
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    main()
