"""One factorizer with reciprocal Fermat/divisor bounds and shared products.

The complete wheel-20 and QR64 conditions are intersected into a period-160
Fermat cursor. Failed midpoint tests shrink the upper bound on the smaller
factor; a certified divisor prefix shrinks the Fermat endpoint. A single
unit-step block polynomial serves divisor-wheel and iterated-average requests.
All components consume one deadline and one optional represented-work cap.
"""
from __future__ import annotations

import argparse
import json
import math
from math import gcd, isqrt
from time import perf_counter

from .common import Budget, BudgetExceeded, SearchResult, ceil_sqrt, checked_result, precheck
from .shared_products import NormalizedProductEngine, recover_mixed_batch
from .vf_algorithm import average_streams, level_order

DIVISOR_RESIDUES = (1, 7, 11, 13, 17, 19, 23, 29)
QR20 = frozenset((i*i) % 20 for i in range(20))
QR64 = frozenset((i*i) % 64 for i in range(64))


def factor_upper_bound(n: int, covered_through: int) -> int:
    """Largest possible smaller divisor after all midpoints <= A are excluded."""
    b = max(covered_through + 1, ceil_sqrt(n))
    return b - ceil_sqrt(b*b - n)


def fermat_midpoint_cap(n: int, lo: int) -> int:
    """Largest possible midpoint if every proper divisor below lo is excluded."""
    if lo < 1:
        raise ValueError("lo must be positive")
    return (lo*lo + n) // (2*lo)


class PythonFermatCursor:
    """Exact arbitrary-precision implementation of the same resumable wheel."""

    def __init__(self, n):
        if not isinstance(n, int) or isinstance(n, bool) or n < 3 or n % 2 == 0:
            raise ValueError("Fermat cursor requires an odd integer >= 3")
        self.n = n
        self.root_ceiling = ceil_sqrt(n)
        self.residues = tuple(r for r in range(160)
                              if (r*r-n) % 20 in QR20 and (r*r-n) % 64 in QR64)
        self.gaps = tuple((self.residues[(i+1) % len(self.residues)] - r) % 160 or 160
                          for i, r in enumerate(self.residues))
        a = self.root_ceiling
        rem = a % 160
        self.index = next((i for i,r in enumerate(self.residues) if r >= rem), 0)
        a += (self.residues[self.index] - rem) % 160
        self.next_midpoint = a
        self.start_midpoint = a
        self.endpoint = (n+1)//2
        self.covered_through = min(a-1,self.endpoint)
        self.difference = a*a-n
        self.candidates = 0
        self.first_midpoint = None
        self.last_midpoint = None
        self.factor = None
        self.half_gap = None
        self.total_kernel_seconds = 0.0
        self.exhausted = False
        self.status = "ready"

    def snapshot(self):
        return dict(n=self.n, factor=self.factor, candidates=self.candidates,
                    square_tests=self.candidates, delta_candidates=0,
                    first_midpoint=self.first_midpoint, last_midpoint=self.last_midpoint,
                    start_midpoint=self.start_midpoint, next_midpoint=self.next_midpoint,
                    covered_through=self.covered_through, half_gap=self.half_gap,
                    root_ceiling=self.root_ceiling, endpoint=self.endpoint,
                    residue_count=len(self.residues), period=160, status=self.status,
                    exhausted=self.exhausted, elapsed_seconds=0.0,
                    total_kernel_seconds=self.total_kernel_seconds)

    def advance(self, max_candidates, timeout_seconds=None, max_midpoint=None):
        if not isinstance(max_candidates, int) or isinstance(max_candidates, bool) or max_candidates < 0:
            raise ValueError("max_candidates must be a nonnegative integer")
        if timeout_seconds is not None and (not isinstance(timeout_seconds, (int,float))
                or isinstance(timeout_seconds,bool) or not math.isfinite(timeout_seconds) or timeout_seconds < 0):
            raise ValueError("timeout_seconds must be finite and nonnegative")
        if max_midpoint is not None and (not isinstance(max_midpoint,int)
                or isinstance(max_midpoint,bool) or max_midpoint < 0):
            raise ValueError("max_midpoint must be a nonnegative integer")
        started, before = perf_counter(), self.candidates
        deadline = None if timeout_seconds is None else started + timeout_seconds
        limit = self.endpoint if max_midpoint is None else min(max_midpoint, self.endpoint)
        self.status = "factor" if self.factor else "exhausted" if self.exhausted else "ready"
        while not self.factor and not self.exhausted:
            used = self.candidates-before
            if used >= max_candidates:
                self.status = "candidate_budget"
                break
            if used % 64 == 0 and deadline is not None and perf_counter() >= deadline:
                self.status = "time_budget"
                break
            if self.next_midpoint > limit:
                self.status = "range_limit"
                break
            a = self.next_midpoint
            self.candidates += 1
            self.last_midpoint = a
            if self.first_midpoint is None:
                self.first_midpoint = a
            b = isqrt(self.difference)
            if b*b == self.difference:
                d = a-b
                if 1 < d < self.n:
                    self.factor, self.half_gap, self.status = d, b, "factor"
                    break
                if d == 1:
                    self.covered_through = self.endpoint
                    self.next_midpoint = self.endpoint+1
                    self.exhausted, self.status = True, "exhausted"
                    break
            gap = self.gaps[self.index]
            self.difference += 2*a*gap + gap*gap
            self.next_midpoint += gap
            self.index = (self.index+1) % len(self.residues)
            self.covered_through = min(self.next_midpoint-1,self.endpoint)
            if self.next_midpoint > self.endpoint:
                self.exhausted, self.status = True, "exhausted"
        elapsed = perf_counter()-started
        self.total_kernel_seconds += elapsed
        out = self.snapshot()
        out.update(delta_candidates=self.candidates-before, elapsed_seconds=elapsed)
        return out

    def close(self):
        pass


def _positive_integer(value, name, minimum=1):
    if not isinstance(value,int) or isinstance(value,bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def factor(n: int, *, timeout=2.0, max_steps=None, native=True, vf_probes=True,
           feedback=True, fermat_enabled=True, shared_polynomial=True,
           initial_degree=16, max_degree=512, fermat_ratio=2.0,
           initial_fermat_tests=64, trace=False) -> SearchResult:
    """Find a proper factor with one shared search state and global budget.

    max_steps counts initial scalar divisor tests, represented AP positions,
    and exact Fermat square tests. Recovery revisits already represented terms.
    Time limits are cooperative; an in-flight polynomial operation can finish
    after the deadline. Common prechecks precede the search work cap.

    vf_probes, feedback, fermat_enabled, and shared_polynomial expose controlled
    ablations. Defaults enable the full synthesis. No reference factors or
    corpus labels enter the algorithm.
    """
    started = perf_counter()
    if timeout is not None and (not isinstance(timeout,(int,float)) or isinstance(timeout,bool)
            or not math.isfinite(timeout) or timeout < 0):
        raise ValueError("timeout must be None or a finite nonnegative number")
    if max_steps is not None:
        _positive_integer(max_steps,"max_steps",0)
    _positive_integer(initial_degree,"initial_degree")
    _positive_integer(max_degree,"max_degree")
    _positive_integer(initial_fermat_tests,"initial_fermat_tests",0)
    if max_degree < initial_degree:
        raise ValueError("max_degree must be >= initial_degree")
    if not isinstance(fermat_ratio,(int,float)) or isinstance(fermat_ratio,bool) or not math.isfinite(fermat_ratio) or fermat_ratio <= 0:
        raise ValueError("fermat_ratio must be finite and positive")
    for name,value in (("native",native),("vf_probes",vf_probes),("feedback",feedback),
                       ("fermat_enabled",fermat_enabled),("shared_polynomial",shared_polynomial),("trace",trace)):
        if not isinstance(value,bool):
            raise ValueError(f"{name} must be a bool")
    budget = Budget(timeout)
    budget.started = started
    budget.deadline = None if timeout is None else started + timeout
    counters = dict(scalar_tests=0, fermat_tests=0, represented_candidates=0,
                    divisor_candidates=0, vf_candidates=0, gcd_calls=0, rounds=0,
                    fermat_chunks=0, fermat_seconds=0.0, product_seconds=0.0,
                    upper_bound_updates=0, lower_bound_updates=0,
                    clipped_strip_positions=0, fermat_bound_reduction=0)
    result = SearchResult("fused",n,counters=counters,
                          details=dict(shared_bounds=feedback, shared_polynomial=shared_polynomial,
                                       vf_probes=vf_probes, fermat_enabled=fermat_enabled,
                                       initial_degree=initial_degree,max_degree=max_degree,
                                       fermat_ratio=fermat_ratio,fermat_period=160))
    history = []
    cursor = None
    lo, hi, square_root = 7, 0, 0
    engines = {}
    streams = None

    def work_used():
        return counters["scalar_tests"] + counters["represented_candidates"] + counters["fermat_tests"]

    def reserve(amount):
        budget.check()
        if max_steps is not None and work_used()+amount > max_steps:
            raise BudgetExceeded("represented-work budget")

    def remaining_seconds():
        return None if budget.deadline is None else max(0.0,budget.deadline-perf_counter())

    def record(event, **extra):
        if trace:
            history.append(dict(event=event,lo=lo,hi=hi,
                                fermat_next=cursor.next_midpoint if cursor else None,
                                fermat_covered=cursor.covered_through if cursor else None,
                                fermat_cap=fermat_midpoint_cap(n,lo) if lo <= square_root else None,
                                work_units=work_used(),**extra))

    def finish(divisor=None,status="factor",**details):
        result.factor, result.status = divisor,status
        result.details.update(details,lo=lo,hi=hi)
        counters["work_units"] = work_used()
        if cursor:
            result.details["fermat_next"] = cursor.next_midpoint
            result.details["fermat_covered"] = cursor.covered_through
        if trace:
            result.details["trace"] = history
        result.elapsed_seconds = perf_counter()-started
        return checked_result(result)

    def update_upper():
        nonlocal hi
        if feedback and cursor:
            previous = hi
            hi = min(hi,factor_upper_bound(n,cursor.covered_through))
            if hi < previous:
                counters["upper_bound_updates"] += 1
                counters["fermat_bound_reduction"] = square_root-hi

    def advance_fermat(maximum, seconds):
        reserve(0)
        if max_steps is not None:
            maximum = min(maximum,max_steps-work_used())
        if maximum <= 0:
            raise BudgetExceeded("represented-work budget")
        remaining = remaining_seconds()
        if remaining is not None:
            seconds = remaining if seconds is None else min(seconds,remaining)
        cap = fermat_midpoint_cap(n,lo if feedback else 7)
        before = perf_counter()
        state = cursor.advance(maximum,timeout_seconds=seconds,max_midpoint=cap)
        counters["fermat_seconds"] += perf_counter()-before
        counters["fermat_tests"] += state["delta_candidates"]
        counters["fermat_chunks"] += 1
        if state["factor"]:
            record("fermat_factor",factor=state["factor"])
            return state["factor"]
        update_upper()
        record("fermat",chunk_tests=state["delta_candidates"],applied_cap=cap)
        budget.check()
        return None

    early = precheck(n)
    if early is not None:
        return finish(early[0],early[1],phase="precheck",precheck=True)
    square_root = isqrt(n)
    hi = square_root
    try:
        reserve(1)
        if fermat_enabled:
            from . import fusion_native
            use_native = native and n < 1 << 64 and fusion_native.is_available()
            cursor = fusion_native.FermatCursor(n) if use_native else PythonFermatCursor(n)
            result.details["fermat_backend"] = "native64" if use_native else "python_bigint"
            update_upper()
        record("init")
        if cursor and initial_fermat_tests:
            found = advance_fermat(initial_fermat_tests,remaining_seconds())
            if found:
                return finish(found,phase="fermat")
        if lo > hi:
            return finish(None,"prime",phase="frontiers_met")

        # One short shared divisor prefix leaves the subsequent strips aligned
        # to the wheel modulo 30. No factor metadata enters this choice.
        prefix_end = min(hi,119)
        for d in range(lo,prefix_end+1):
            if d % 30 not in DIVISOR_RESIDUES:
                continue
            reserve(1)
            counters["scalar_tests"] += 1
            counters["divisor_candidates"] += 1
            counters["gcd_calls"] += 1
            g = gcd(d,n)
            if 1 < g < n:
                record("prefix_factor",factor=g)
                return finish(g,phase="divisor_prefix")
        lo = prefix_end+1
        counters["lower_bound_updates"] += 1
        record("divisor_prefix")
        degree = initial_degree

        while lo <= hi:
            reserve(1)
            counters["rounds"] += 1
            tick = perf_counter()
            # About m block evaluations represent m^2 candidates. Eight
            # step-30 streams cover a contiguous strip of possible divisors.
            available = max(1,(hi-lo+1)*8//30)
            m = min(degree,max(1,isqrt(available)))
            blocks_per_lane = max(1,m//8)
            base = (lo//30)*30
            nominal_end = base + 30*m*blocks_per_lane-1
            end = min(hi,nominal_end)
            requests, tails = [], []
            for residue in DIVISOR_RESIDUES:
                first = base+residue
                if first < lo:
                    first += 30
                count = max(0,(end-first)//30+1)
                full,tail = divmod(count,m)
                for j in range(full):
                    requests.append(dict(start=first+30*m*j,step=30,count=m,kind="divisor"))
                for j in range(tail):
                    tails.append(first+30*m*full+30*j)
            divisor_blocks = len(requests)
            if vf_probes:
                if streams is None:
                    levels = level_order(n,"published")
                    streams = average_streams(n,levels,2*square_root+2)
                    for stream in streams:
                        stream["cursor"] = 0
                    result.details["vf_levels"] = levels
                for stream in streams:
                    if stream["cursor"]+m <= stream["count"]:
                        requests.append(dict(start=stream["start"]+5*stream["cursor"],
                                             step=5,count=m,kind="vf"))
                        stream["cursor"] += m
            vf_blocks = len(requests)-divisor_blocks
            reserve(len(requests)*m + len(tails))
            record("product_plan",degree=m,strip_end=end,nominal_strip_end=nominal_end,
                   divisor_blocks=divisor_blocks,vf_blocks=vf_blocks,tail_count=len(tails))
            if requests:
                engine = engines.get(m)
                if engine is None:
                    engine = NormalizedProductEngine(n,m,budget=budget,counters=counters,
                                                     shared=shared_polynomial)
                    engines[m] = engine
                values = engine.products([r["start"] for r in requests],[r["step"] for r in requests])
                counters["divisor_candidates"] += divisor_blocks*m
                counters["vf_candidates"] += vf_blocks*m
                found = recover_mixed_batch(n,requests,values,budget=budget,counters=counters)
                if found:
                    counters["product_seconds"] += perf_counter()-tick
                    record("product_factor",factor=found,degree=m,
                           divisor_blocks=divisor_blocks,vf_blocks=vf_blocks)
                    return finish(found,phase="shared_product" if vf_blocks else "divisor_product")
            for d in tails:
                reserve(1)
                counters["scalar_tests"] += 1
                counters["divisor_candidates"] += 1
                counters["gcd_calls"] += 1
                found = gcd(d,n)
                if 1 < found < n:
                    counters["product_seconds"] += perf_counter()-tick
                    record("tail_factor",factor=found)
                    return finish(found,phase="divisor_tail")
            # This update is atomic: all wheel residues and final tails in the
            # numeric strip have now been certified, including exception recovery.
            lo = end+1
            counters["lower_bound_updates"] += 1
            if end < nominal_end:
                counters["clipped_strip_positions"] += max(0,(nominal_end-end)*8//30)
            elapsed = perf_counter()-tick
            counters["product_seconds"] += elapsed
            record("product",degree=m,strip_end=end,divisor_blocks=divisor_blocks,vf_blocks=vf_blocks)
            budget.check()
            if lo > hi:
                break
            if cursor:
                allocation = min(.02,max(.00005,fermat_ratio*elapsed))
                found = advance_fermat((1 << 63)-1,allocation)
                if found:
                    return finish(found,phase="fermat")
            degree = min(2*degree,max_degree)
        return finish(None,"prime",phase="frontiers_met")
    except BudgetExceeded as exc:
        return finish(None,"budget",reason=str(exc),phase="budget")
    finally:
        if cursor:
            cursor.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n",type=int)
    parser.add_argument("--timeout",type=float,default=10.0)
    parser.add_argument("--max-steps",type=int)
    parser.add_argument("--python",action="store_true",help="use arbitrary-precision Python Fermat cursor")
    parser.add_argument("--trace",action="store_true")
    parser.add_argument("--full",action="store_true",help="recursively factor with one overall deadline")
    args=parser.parse_args()
    options=dict(max_steps=args.max_steps,native=not args.python,trace=args.trace)
    if not args.full:
        print(json.dumps(factor(args.n,timeout=args.timeout,**options).to_dict(),indent=2))
        return
    _positive_integer(args.n,"n",2)
    started=perf_counter()
    pending,primes,unresolved,runs=[args.n],[],[],[]
    work_remaining=args.max_steps
    while pending:
        n=pending.pop()
        remaining=max(0.0,args.timeout-(perf_counter()-started))
        run=factor(n,timeout=remaining,**{**options,"max_steps":work_remaining})
        runs.append(run.to_dict())
        if work_remaining is not None:
            work_remaining=max(0,work_remaining-run.counters.get("work_units",0))
        if run.factor:
            pending.extend((run.factor,n//run.factor))
        elif run.status=="prime":
            primes.append(n)
        else:
            unresolved.append(n)
    print(json.dumps(dict(input=args.n,prime_factors=sorted(primes),unresolved=sorted(unresolved),
                          elapsed_seconds=perf_counter()-started,runs=runs),indent=2))


if __name__=="__main__":
    main()
