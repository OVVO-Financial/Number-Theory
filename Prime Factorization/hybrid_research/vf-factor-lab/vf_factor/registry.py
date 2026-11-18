"""Named algorithm entry points used by both CLI and benchmark runner."""
from .baselines import fermat, pollard_brent, published_vf, trial_division
from .common import Budget, SearchResult, checked_result, precheck
from .vf_algorithm import product_trial, vf_fast, vf_fast_multi, vf_first, vf_hybrid


def flint_native(n, *, timeout=2.0, **kwargs):
    """Reference FLINT full factorization; native C call is not preemptible.

    The benchmark uses this only on <=64-bit inputs. It reports full
    factorization time, which is at least as much work as our one-factor goal.
    """
    from flint import fmpz
    budget = Budget(timeout)
    result = SearchResult("flint_native", n)
    early = precheck(n)
    if early is not None:
        result.factor, result.status = early
        result.details["precheck"] = True
    else:
        factors = [(int(p), int(e)) for p, e in fmpz(n).factor()]
        result.details["full_factorization"] = factors
        result.details["native_call_timeout_not_preemptible"] = True
        for p, _ in factors:
            if 1 < p < n:
                result.factor = p
                break
        if result.factor is None:
            result.status = "prime"
    result.elapsed_seconds = budget.elapsed
    return checked_result(result)


ALGORITHMS = dict(published_vf=published_vf, vf_fast=vf_fast,
                  vf_fast_multi=vf_fast_multi, vf_hybrid=vf_hybrid,
                  vf_first=vf_first,
                  product_trial=product_trial, trial_division=trial_division,
                  fermat=fermat, pollard_brent=pollard_brent,
                  flint_native=flint_native)


def get_algorithms():
    methods = dict(ALGORITHMS)
    try:
        from . import native
        if native.is_available():
            for name in ("published_vf", "pollard_brent", "trial_division", "fermat"):
                if hasattr(native, name):
                    methods["native_" + name] = getattr(native, name)
    except (ImportError, OSError):
        pass
    return methods
