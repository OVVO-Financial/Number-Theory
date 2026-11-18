"""Exact complete Fermat wheel-20 search and matched native controls.

The verified terminal-digit table is indexed by N modulo 20. Its projection
onto the midpoint R has period ten, with two or three retained residues when
gcd(N,10)=1. A cyclic gap table visits exactly those midpoints in order.

The plain parity and wheel methods use the same integer square-root code.
The separately named QR64 variants add the same square-residue prefilter to
both, isolating the wheel's incremental benefit. No known factors are input.
Build the optional kernel with ``python3 build_wheel20.py``.
"""
from __future__ import annotations

import ctypes
import math
from pathlib import Path
from time import perf_counter

from .common import SearchResult, checked_result, precheck


def admissible_pairs(n_mod20: int) -> tuple[tuple[int, int], ...]:
    """Return the complete table of (R mod 10, i mod 10) for one odd class."""
    if not isinstance(n_mod20, int) or isinstance(n_mod20, bool) or n_mod20 % 2 == 0:
        raise ValueError("the complete odd-N table requires an odd integer class")
    residue = n_mod20 % 20
    return tuple((r, b) for r in range(10) for b in range(10)
                 if (r * r - b * b) % 20 == residue)


def midpoint_residues(n: int) -> tuple[int, ...]:
    """Project the exact pair table onto R modulo ten."""
    return tuple(sorted({r for r, _ in admissible_pairs(n)}))


class _WheelResult(ctypes.Structure):
    _fields_ = [
        ("factor", ctypes.c_uint64),
        ("candidates", ctypes.c_uint64),
        ("square_tests", ctypes.c_uint64),
        ("first_midpoint", ctypes.c_uint64),
        ("last_midpoint", ctypes.c_uint64),
        ("start_midpoint", ctypes.c_uint64),
        ("root_ceiling", ctypes.c_uint64),
        ("half_gap", ctypes.c_uint64),
        ("residue_mask", ctypes.c_uint32),
        ("residue_count", ctypes.c_uint32),
        ("status", ctypes.c_uint32),
        ("error", ctypes.c_uint32),
        ("elapsed_seconds", ctypes.c_double),
    ]


_LIBRARY = None


def _load():
    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    directory = Path(__file__).resolve().parent
    for suffix in (".so", ".dylib", ".dll"):
        path = directory / ("wheel20_kernel" + suffix)
        if path.exists():
            lib = ctypes.CDLL(str(path))
            break
    else:
        raise RuntimeError("Wheel-20 kernel is not built; run python3 build_wheel20.py")
    lib.wheel20_abi_version.restype = ctypes.c_uint
    lib.wheel20_result_size.restype = ctypes.c_uint
    lib.wheel20_compiler.restype = ctypes.c_char_p
    lib.wheel20_isqrt.argtypes = [ctypes.c_uint64]
    lib.wheel20_isqrt.restype = ctypes.c_uint64
    lib.wheel20_run.argtypes = [ctypes.c_uint64, ctypes.c_uint, ctypes.c_double,
                                ctypes.c_uint64, ctypes.POINTER(_WheelResult)]
    lib.wheel20_run.restype = ctypes.c_int
    if lib.wheel20_abi_version() != 1 or lib.wheel20_result_size() != ctypes.sizeof(_WheelResult):
        raise RuntimeError("Wheel-20 ABI mismatch; rebuild with python3 build_wheel20.py")
    _LIBRARY = lib
    return lib


def is_available() -> bool:
    try:
        _load()
        return True
    except (RuntimeError, OSError):
        return False


def backend_info() -> dict:
    lib = _load()
    return {
        "backend": "C++17 native64 via ctypes",
        "compiler": lib.wheel20_compiler().decode(),
        "abi_version": lib.wheel20_abi_version(),
        "integer_range": "2 <= n < 2**64",
        "square_arithmetic": "unsigned __int128 and exact integer Newton square roots",
        "wheel_table_index": "N modulo 20",
        "midpoint_period": 10,
        "timing": "includes common prechecks, per-input gap-table setup, FFI and result construction",
    }


_NAMES = ("native_fermat_parity", "native_fermat_wheel20",
          "native_fermat_parity_qr64", "native_fermat_wheel20_qr64")


def _run(n, mode, timeout, max_steps):
    started = perf_counter()
    if not isinstance(n, int) or isinstance(n, bool) or not 2 <= n < 1 << 64:
        raise ValueError("native64 requires an integer 2 <= n < 2**64")
    if max_steps is not None and (not isinstance(max_steps, int) or isinstance(max_steps, bool)
                                  or not 0 <= max_steps < 1 << 64):
        raise ValueError("max_steps must be None or an integer from 0 through 2**64-1")
    if timeout is not None and (not isinstance(timeout, (int, float))
                                or not math.isfinite(timeout)):
        raise ValueError("timeout must be None or a finite number of seconds")
    wheel, qr64 = bool(mode & 1), bool(mode & 2)
    result = SearchResult(algorithm=_NAMES[mode], n=n,
                          details={"backend": "native64_wheel20", "fallback_used": False,
                                   "wheel_modulus": 20 if wheel else 4,
                                   "midpoint_period": 10 if wheel else 2, "qr64": qr64})
    found = precheck(n)
    if found is not None:
        result.factor, result.status = found
        result.details["precheck"] = True
        result.elapsed_seconds = perf_counter() - started
        return checked_result(result)
    lib = _load()
    remaining = -1.0 if timeout is None else max(0.0, timeout - (perf_counter() - started))
    out = _WheelResult()
    error = lib.wheel20_run(n, mode, remaining,
                           (1 << 64) - 1 if max_steps is None else max_steps,
                           ctypes.byref(out))
    if error or out.status == 3:
        raise RuntimeError(f"Wheel-20 search failed with code {out.error}")
    result.factor = int(out.factor) if out.factor else None
    result.status = ("budget", "factor", "prime")[out.status]
    count = int(out.candidates)
    result.counters.update(
        candidates=count, iterations=count, square_tests=int(out.square_tests),
        residue_rejections=int(out.candidates - out.square_tests),
        scan_distance=int(out.last_midpoint - out.root_ceiling) if count else 0,
        midpoint_span=int(out.last_midpoint - out.first_midpoint + 1) if count else 0,
    )
    result.details.update(
        kernel_seconds=out.elapsed_seconds,
        first_midpoint=int(out.first_midpoint) if count else None,
        last_midpoint=int(out.last_midpoint) if count else None,
        start_midpoint=int(out.start_midpoint), ceil_sqrt=int(out.root_ceiling),
        midpoint_residues=tuple(r for r in range(10 if wheel else 2)
                               if out.residue_mask & (1 << r)),
        source="Complete_Fermat_Sieve_Verified.pdf",
    )
    if result.status == "budget":
        result.details["budget_reason"] = "time budget" if out.error == 1 else "step budget"
    if result.factor:
        result.details.update(midpoint=int(out.last_midpoint), half_gap=int(out.half_gap))
    if result.status == "prime":
        result.details["coverage"] = "complete deterministic Fermat search exhausted"
    result.elapsed_seconds = perf_counter() - started
    return checked_result(result)


def native_fermat_parity(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 0, timeout, max_steps)


def native_fermat_wheel20(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 1, timeout, max_steps)


def native_fermat_parity_qr64(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 2, timeout, max_steps)


def native_fermat_wheel20_qr64(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 3, timeout, max_steps)
