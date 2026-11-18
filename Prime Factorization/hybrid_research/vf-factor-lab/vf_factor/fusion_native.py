"""Resumable native Fermat cursor with a complete combined period-160 wheel.

There are no primality or small-factor prechecks here: the unified coordinator
performs those once. The cursor preserves its wheel, next midpoint, difference
of squares and coverage certificate across calls.

For 64-bit square arguments, a hardware estimate is corrected with exact
integer comparisons. Larger square arguments use an integer Newton fallback.

``covered_through`` certifies that every midpoint from ceil(sqrt(n)) through
that value has been ruled out as a proper Fermat factorization. Congruence-only
gaps are included without square tests. A successful midpoint is never included
in excluded coverage. A temporary ``max_midpoint`` does not destroy the state.
"""
from __future__ import annotations

import ctypes
import math
from pathlib import Path
from time import perf_counter

_UINT64_MAX = (1 << 64) - 1
_STATUS = ("ready", "factor", "exhausted", "candidate_budget", "time_budget", "range_limit")
_LIBRARY = None


class _Snapshot(ctypes.Structure):
    _fields_ = [
        ("n", ctypes.c_uint64), ("factor", ctypes.c_uint64),
        ("candidates", ctypes.c_uint64), ("square_tests", ctypes.c_uint64),
        ("delta_candidates", ctypes.c_uint64),
        ("first_midpoint", ctypes.c_uint64), ("last_midpoint", ctypes.c_uint64),
        ("start_midpoint", ctypes.c_uint64), ("next_midpoint", ctypes.c_uint64),
        ("covered_through", ctypes.c_uint64), ("half_gap", ctypes.c_uint64),
        ("root_ceiling", ctypes.c_uint64), ("endpoint", ctypes.c_uint64),
        ("residue_count", ctypes.c_uint32), ("period", ctypes.c_uint32),
        ("status", ctypes.c_uint32), ("exhausted", ctypes.c_uint32),
        ("elapsed_seconds", ctypes.c_double), ("total_kernel_seconds", ctypes.c_double),
    ]


def _load():
    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    directory = Path(__file__).resolve().parent
    for suffix in (".so", ".dylib", ".dll"):
        path = directory / ("fusion_kernel" + suffix)
        if path.exists():
            lib = ctypes.CDLL(str(path))
            break
    else:
        raise RuntimeError("Fusion kernel is not built; run python3 build_fusion.py")
    lib.fusion_abi_version.restype = ctypes.c_uint
    lib.fusion_snapshot_size.restype = ctypes.c_uint
    lib.fusion_compiler.restype = ctypes.c_char_p
    lib.fusion_isqrt.argtypes = [ctypes.c_uint64]
    lib.fusion_isqrt.restype = ctypes.c_uint64
    lib.fusion_fermat_create.argtypes = [ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint)]
    lib.fusion_fermat_create.restype = ctypes.c_void_p
    lib.fusion_fermat_destroy.argtypes = [ctypes.c_void_p]
    lib.fusion_fermat_destroy.restype = None
    lib.fusion_fermat_snapshot.argtypes = [ctypes.c_void_p, ctypes.POINTER(_Snapshot)]
    lib.fusion_fermat_snapshot.restype = ctypes.c_int
    lib.fusion_fermat_advance.argtypes = [ctypes.c_void_p, ctypes.c_uint64, ctypes.c_double,
                                         ctypes.c_uint64, ctypes.POINTER(_Snapshot)]
    lib.fusion_fermat_advance.restype = ctypes.c_int
    if lib.fusion_abi_version() != 1 or lib.fusion_snapshot_size() != ctypes.sizeof(_Snapshot):
        raise RuntimeError("Fusion ABI mismatch; rebuild with python3 build_fusion.py")
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
        "backend": "resumable C++17 native64 via ctypes",
        "compiler": lib.fusion_compiler().decode(),
        "abi_version": lib.fusion_abi_version(),
        "midpoint_period": 160,
        "filters": "complete N mod20 paired sieve intersected with QR64",
        "integer_range": "odd integers 3 <= n < 2**64",
        "square_arithmetic": "hardware estimate with exact integer correction; integer Newton fallback above 64 bits; unsigned __int128 products",
        "coverage": "failed square tests and necessary-congruence gaps; excludes any successful midpoint",
    }


def _integer(value, name, *, minimum=0):
    if not isinstance(value, int) or isinstance(value, bool) or not minimum <= value <= _UINT64_MAX:
        raise ValueError(f"{name} must be an integer from {minimum} through 2**64-1")
    return value


def _decode(out: _Snapshot) -> dict:
    result = {name: int(getattr(out, name)) for name, kind in out._fields_
              if kind not in (ctypes.c_double,)}
    result.update(status=_STATUS[out.status], factor=int(out.factor) or None,
                  exhausted=bool(out.exhausted), elapsed_seconds=float(out.elapsed_seconds),
                  total_kernel_seconds=float(out.total_kernel_seconds),
                  first_midpoint=int(out.first_midpoint) if out.candidates else None,
                  last_midpoint=int(out.last_midpoint) if out.candidates else None,
                  half_gap=int(out.half_gap) if out.factor else None)
    return result


class FermatCursor:
    """Persistent complete-sieve Fermat state; close it or use a context manager.

    Candidate limits are per call. Timeouts are relative, cooperative budgets
    for one advance call; the coordinator supplies its remaining global budget.
    A midpoint limit only pauses the cursor. Congruence-only certified gaps may
    extend past that temporary limit, but no square test is performed past it.
    Zero work or time budgets perform no square tests and preserve coverage.
    """

    def __init__(self, n: int):
        started = perf_counter()
        _integer(n, "n", minimum=3)
        if n % 2 == 0:
            raise ValueError("FermatCursor requires odd n")
        self.n = n
        self._library = _load()
        error = ctypes.c_uint()
        self._pointer = self._library.fusion_fermat_create(n, ctypes.byref(error))
        if not self._pointer:
            raise RuntimeError(f"Fusion cursor initialization failed with code {error.value}")
        self._output = _Snapshot()
        self._library.fusion_fermat_snapshot(self._pointer, ctypes.byref(self._output))
        self._last = _decode(self._output)
        self.setup_seconds = perf_counter() - started

    @property
    def covered_through(self):
        return self._last["covered_through"]

    @property
    def next_midpoint(self):
        return self._last["next_midpoint"]

    @property
    def candidates(self):
        return self._last["candidates"]

    def snapshot(self) -> dict:
        if not self._pointer:
            raise RuntimeError("Fermat cursor is closed")
        return dict(self._last)

    def advance(self, max_candidates: int, timeout_seconds: float | None = None,
                max_midpoint: int | None = None) -> dict:
        if not self._pointer:
            raise RuntimeError("Fermat cursor is closed")
        maximum = _integer(max_candidates, "max_candidates")
        limit = _UINT64_MAX if max_midpoint is None else _integer(max_midpoint, "max_midpoint")
        if timeout_seconds is None:
            seconds = -1.0
        elif (not isinstance(timeout_seconds, (int, float)) or isinstance(timeout_seconds, bool)
              or not math.isfinite(timeout_seconds) or timeout_seconds < 0):
            raise ValueError("timeout_seconds must be None or a finite nonnegative number")
        else:
            seconds = float(timeout_seconds)
        error = self._library.fusion_fermat_advance(
            self._pointer, maximum, seconds, limit, ctypes.byref(self._output))
        if error:
            raise RuntimeError(f"Fusion cursor advance failed with code {error}")
        self._last = _decode(self._output)
        return dict(self._last)

    def close(self):
        pointer = getattr(self, "_pointer", None)
        if pointer:
            self._library.fusion_fermat_destroy(pointer)
            self._pointer = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def __del__(self):
        self.close()

