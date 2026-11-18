"""Reproduce the complete Fermat wheel-20 comparison using the original corpus.

The four Fermat variants share one C++ kernel. Two independent switches select
the complete wheel and the modulo-64 square filter. Three supplemental cases
check the omitted-lane regression and a pair of inputs with identical local
Fermat lane data but very different factor locations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from vf_factor import benchmark, native, wheel20
from vf_factor.vf_algorithm import product_trial

ROOT = Path(__file__).resolve().parent


def methods():
    return {
        "native_fermat_parity": wheel20.native_fermat_parity,
        "native_fermat_wheel20": wheel20.native_fermat_wheel20,
        "native_fermat_parity_qr64": wheel20.native_fermat_parity_qr64,
        "native_fermat_wheel20_qr64": wheel20.native_fermat_wheel20_qr64,
        "product_trial": product_trial,
        "native_pollard_brent": native.pollard_brent,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=1.0)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--outdir", default=str(ROOT / "results/wheel20"))
    args = parser.parse_args()
    cases = json.loads((ROOT / "results/main/benchmark-corpus.json").read_text())
    for case_id, p, q in (
        ("wheel20-regression-551", 19, 29),
        ("same-lane-near", 100459, 100469),
        ("same-lane-far", 19, 531211069),
    ):
        n = p * q
        cases.append(dict(case_id=case_id, family="wheel20_supplement", n=n,
                          bits=n.bit_length(), target_bits=n.bit_length(),
                          p=p, q=q, factorization=[[p, 1], [q, 1]]))
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = ("vf_factor/wheel20.py", "vf_factor/wheel20_kernel.cpp", "build_wheel20.py",
             "wheel20_benchmark.py", "vf_factor/benchmark.py", "vf_factor/poly_engine.py",
             "vf_factor/vf_algorithm.py")
    protocol = {
        "source_url": "https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Complete_Fermat_Sieve_Verified.pdf",
        "baseline": "Original product_trial unchanged; degree cap 4096; per-input setup timed",
        "comparison": "Matched native parity/wheel20, each with and without QR64",
        "corpus": "Original 83 fixed cases, followed by 3 preselected regressions/structural examples",
        "main_analysis": "Report original corpus separately from supplemental examples",
        "timeout_seconds": args.timeout,
        "completed_repeats": args.repeats,
        "algorithm_order_seed": 20261118,
        "native": wheel20.backend_info(),
        "source_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in files},
    }
    (out / "protocol.json").write_text(json.dumps(protocol, indent=2))
    registry = methods()
    rows = benchmark.run(cases, list(registry), timeout=args.timeout, repeats=args.repeats,
                         seed=20261118, outdir=str(out), algorithm_registry=registry)
    assert not any(r["status"] == "error" for r in rows), "Benchmark contains errors; inspect raw log"


if __name__ == "__main__":
    main()
