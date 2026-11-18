"""Optional compiled 64-bit baselines for fair comparison with native FLINT.

Build with ``python3 build_native.py``. These wrappers include shared Python
prechecks, FFI overhead, and result construction in ``elapsed_seconds``.
``kernel_seconds`` separately reports time spent in the compiled search.
Native Brent uses C++ mt19937_64, not Python's Random generator.
"""
from __future__ import annotations

import ctypes
from pathlib import Path
from time import perf_counter

from .common import SearchResult, checked_result, precheck


class _NativeResult(ctypes.Structure):
    _fields_ = [
        ("factor", ctypes.c_uint64),
        ("steps", ctypes.c_uint64),
        ("gcd_calls", ctypes.c_uint64),
        ("modular_multiplications", ctypes.c_uint64),
        ("square_tests", ctypes.c_uint64),
        ("rounds", ctypes.c_uint64),
        ("restarts", ctypes.c_uint64),
        ("trivial_gcd_hits", ctypes.c_uint64),
        ("auxiliary", ctypes.c_uint64),
        ("argument_abs", ctypes.c_uint64),
        ("status", ctypes.c_uint32),
        ("selected_level", ctypes.c_uint32),
        ("direction", ctypes.c_int32),
        ("offset", ctypes.c_uint32),
        ("argument_negative", ctypes.c_uint32),
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
        path = directory / ("native_kernel" + suffix)
        if path.exists():
            lib = ctypes.CDLL(str(path))
            break
    else:
        raise RuntimeError("Native kernel is not built; run python3 build_native.py")
    lib.vf_native_abi_version.restype = ctypes.c_uint
    lib.vf_native_result_size.restype = ctypes.c_uint
    lib.vf_native_compiler.restype = ctypes.c_char_p
    lib.vf_native_isqrt.argtypes = [ctypes.c_uint64]
    lib.vf_native_isqrt.restype = ctypes.c_uint64
    lib.vf_native_run.argtypes = [
        ctypes.c_uint64, ctypes.c_uint, ctypes.c_double, ctypes.c_uint64,
        ctypes.c_uint64, ctypes.c_uint, ctypes.POINTER(_NativeResult),
    ]
    lib.vf_native_run.restype = ctypes.c_int
    if lib.vf_native_abi_version() != 1 or lib.vf_native_result_size() != ctypes.sizeof(_NativeResult):
        raise RuntimeError("Native ABI mismatch; rebuild with python3 build_native.py")
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
        "compiler": lib.vf_native_compiler().decode(),
        "abi_version": lib.vf_native_abi_version(),
        "integer_range": "2 <= n < 2**64",
        "modular_products": "unsigned __int128",
        "rho_rng": "std::mt19937_64",
        "timing": "wrapper wall time includes common prechecks and FFI; kernel time reported separately",
    }


_NAMES = ("published_vf", "pollard_brent", "trial_division", "fermat")


def _run(n, algorithm, timeout, max_steps, seed=0, batch_size=128):
    started = perf_counter()
    if not isinstance(n, int) or isinstance(n, bool) or not 2 <= n < 1 << 64:
        raise ValueError("native64 requires an integer 2 <= n < 2**64")
    if max_steps is not None and (not isinstance(max_steps, int) or not 0 <= max_steps < 1 << 64):
        raise ValueError("max_steps must be None or an integer from 0 through 2**64-1")
    if not isinstance(batch_size, int) or not 1 <= batch_size < 1 << 32:
        raise ValueError("batch_size must be a positive uint32")
    result = SearchResult(algorithm=_NAMES[algorithm] + "_native64", n=n,
                          details={"backend": "native64", "fallback_used": False})
    found = precheck(n)
    if found is not None:
        result.factor, result.status = found
        result.details["precheck"] = True
        result.elapsed_seconds = perf_counter() - started
        return checked_result(result)
    lib = _load()
    # Include prechecks/loading in the wall-time allowance. A negative value
    # is reserved by the C ABI for the explicitly unlimited None case.
    remaining = -1.0 if timeout is None else max(0.0, timeout - (perf_counter() - started))
    out = _NativeResult()
    error = lib.vf_native_run(
        n, algorithm, remaining, (1 << 64) - 1 if max_steps is None else max_steps,
        seed % (1 << 64), batch_size, ctypes.byref(out),
    )
    if error or out.status == 3:
        raise RuntimeError(f"Native search failed with code {out.error}")
    result.factor = int(out.factor) if out.factor else None
    result.status = ("budget", "factor", "prime")[out.status]
    result.counters.update(candidates=int(out.steps), gcd_calls=int(out.gcd_calls))
    result.details["kernel_seconds"] = out.elapsed_seconds
    if result.status == "budget":
        result.details["budget_reason"] = "time budget" if out.error == 1 else "step budget"
    if algorithm == 0:
        result.counters.update(rounds=int(out.rounds), trivial_gcd_hits=int(out.trivial_gcd_hits),
                               scan_distance=int(out.auxiliary))
        result.details.update(selected_level=int(out.selected_level),
                              source="Iterated Averages.md: VF(n, 5)", exact_arithmetic=True)
        if result.factor:
            result.details.update(
                realized_b=int(out.auxiliary),
                gcd_argument=int(out.argument_abs) * (-1 if out.argument_negative else 1),
                direction="descending" if out.direction < 0 else "ascending",
                offset=int(out.offset),
            )
    elif algorithm == 1:
        result.counters.update(map_evaluations=int(out.steps),
                               modular_multiplications=int(out.modular_multiplications),
                               restarts=int(out.restarts),
                               degenerate_blocks=int(out.trivial_gcd_hits),
                               recovery_steps=int(out.auxiliary))
        result.details.update(seed=seed, batch_size=batch_size, rng="std::mt19937_64")
    elif algorithm == 2:
        result.counters["divisor_tests"] = int(out.steps)
        result.details["wheel"] = 30
    elif algorithm == 3:
        result.counters.update(square_tests=int(out.square_tests),
                               residue_rejections=int(out.steps - out.square_tests))
        if result.factor:
            result.details.update(midpoint=int(out.argument_abs), half_gap=int(out.auxiliary))
    if result.status == "prime":
        result.details["coverage"] = "complete deterministic search exhausted"
    result.elapsed_seconds = perf_counter() - started
    return checked_result(result)


def published_vf(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 0, timeout, max_steps)


def pollard_brent(n, *, timeout=2.0, max_steps=None, seed=0, batch_size=128, **_options):
    return _run(n, 1, timeout, max_steps, seed, batch_size)


def trial_division(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 2, timeout, max_steps)


def fermat(n, *, timeout=2.0, max_steps=None, **_options):
    return _run(n, 3, timeout, max_steps)
