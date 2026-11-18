#pragma once

#include "gmp_factor.hpp"
#include "gmp_product_worker.hpp"
#include "gmp_regions.hpp"
#include <atomic>
#include <condition_variable>
#include <ctime>
#include <deque>
#include <exception>
#include <mutex>
#include <thread>

namespace vf {

struct ParallelOptions {
    Options algorithm;
    std::size_t threads = 4;
    std::size_t degree = 128;
    std::uint64_t fermat_span = 655360;
    double job_seconds = .025;
    std::optional<Z> factor_lo, factor_hi;
    void validate() const {
        algorithm.validate();
        if (!threads || threads > 256) throw std::invalid_argument("threads must be between 1 and 256");
        if (!degree || degree > 65536) throw std::invalid_argument("degree must be between 1 and 65536");
        if (!fermat_span) throw std::invalid_argument("Fermat span must be positive");
        if (!std::isfinite(job_seconds) || job_seconds <= 0)
            throw std::invalid_argument("job time quantum must be finite and positive");
        if ((factor_lo && *factor_lo < 2) || (factor_hi && *factor_hi < 2) ||
            (factor_lo && factor_hi && *factor_lo > *factor_hi))
            throw std::invalid_argument("factor band must have 2 <= lower <= upper");
    }
};

struct ParallelTrace {
    std::string event, kind, status;
    std::uint64_t job_id = 0, worker_id = 0, job_work = 0, work_units = 0;
    Z requested_lo = 0, requested_hi = -1;
    std::optional<Z> certified_lo, certified_hi, lo, hi;
};

struct ParallelResult {
    Z n, factor = 0, search_lo = 7, search_hi = 0;
    std::optional<Z> lo, hi;
    std::string status = "budget", reason, phase;
    Counters counters;
    std::vector<ParallelTrace> trace;
    std::vector<ClosedInterval> certificates;
    std::size_t threads = 0;
    double elapsed_seconds = 0, process_cpu_seconds = 0;
    std::uint64_t outstanding_reservations = 0, seed_work_units = 0;
};

namespace regional_detail {

struct Job {
    std::uint64_t id = 0, quota = 0;
    bool fermat = false;
    Z lo, hi;
    std::vector<Request> auxiliary;
};

struct JobResult {
    Z factor = 0;
    std::optional<ClosedInterval> certificate;
    Z covered_coordinate = -1;
    std::string status = "cancelled";
    Counters counters;
    std::uint64_t work = 0;
    double elapsed_seconds = 0;
};

inline void merge_counters(Counters& total, const Counters& part) {
    for (const auto& [key,value] : part.counts) {
        if (key == "work_units") continue;
        if (key == "polynomial_degree") total.counts[key] = std::max(total.get(key),value);
        else total.inc(key,value);
    }
    for (const auto& [key,value] : part.seconds) total.add_seconds(key,value);
}

} // namespace regional_detail

// Every mutable GMP object used for arithmetic belongs to one worker. The
// coordinator alone owns the ledger, reservations, result and trace, under a
// mutex. Immutable N may be read concurrently, as permitted by GMP.
inline ParallelResult parallel_find_factor(const Z& n, const ParallelOptions& options,
                                           Budget& global_time,
                                           std::optional<std::uint64_t> max_work) {
    options.validate();
    if (n < 2) throw std::invalid_argument("N must be an integer >= 2");
    if (global_time.remaining() < 0)
        throw std::invalid_argument("parallel_find_factor requires a finite global deadline");
    const auto started = Clock::now();
    const auto cpu_started = std::clock();
    ParallelResult result; result.n = n; result.threads = options.threads;
    auto& counts = result.counters;
    for (const char* key : {"work_units","fermat_tests","scalar_tests","represented_candidates",
         "jobs_started","jobs_finished","jobs_completed","jobs_cancelled","jobs_budget",
         "factor_jobs","reservation_high_water","peak_inflight_work","active_workers_high_water",
         "product_jobs","fermat_jobs","certificates_committed"}) counts.counts[key] = 0;
    auto finalize_time = [&]() {
        result.elapsed_seconds = seconds_since(started);
        result.process_cpu_seconds = double(std::clock()-cpu_started)/CLOCKS_PER_SEC;
    };

    // Invoke only the established common prechecks, with zero search work.
    Options pre = options.algorithm;
    pre.fermat_enabled = false; pre.vf_probes = false; pre.trace = false;
    pre.max_steps = 0; pre.timeout = std::max(0.0,global_time.remaining());
    const Result initial = factor(n,pre);
    if (initial.status != "budget") {
        result.factor = initial.factor; result.status = initial.status;
        result.phase = "precheck"; finalize_time(); return result;
    }
    const Z root = floor_sqrt(n);
    const Z lower = std::max(Z(7), options.factor_lo.value_or(Z(7)));
    const Z upper = std::min(root,options.factor_hi.value_or(root));
    if (lower > upper) throw std::invalid_argument("factor band does not intersect 7..floor(sqrt(N))");
    result.search_lo = lower; result.search_hi = upper;
    const bool full_domain = lower == 7 && upper == root;
    CoverageLedger ledger(lower,upper);
    const auto midpoint_band = midpoint_interval_for_factors(n,lower,upper);
    Z next_product = lower;
    Z next_fermat = midpoint_band ? midpoint_band->lo : Z(1);
    Z fermat_end = midpoint_band ? midpoint_band->hi : Z(0);
    std::uint64_t spent = 0, reserved = 0;
    const std::uint64_t allowance = max_work.value_or(std::numeric_limits<std::uint64_t>::max());
    auto record = [&](const std::string& event, const std::string& kind,
                      const std::string& status, std::uint64_t job_id,
                      std::uint64_t worker_id, const Z& requested_lo, const Z& requested_hi,
                      const std::optional<ClosedInterval>& certificate, std::uint64_t work) {
        if (!options.algorithm.trace) return;
        ParallelTrace trace;
        trace.event=event; trace.kind=kind; trace.status=status;
        trace.job_id=job_id; trace.worker_id=worker_id;
        trace.requested_lo=requested_lo; trace.requested_hi=requested_hi;
        trace.job_work=work; trace.work_units=spent;
        if (certificate) { trace.certified_lo=certificate->lo; trace.certified_hi=certificate->hi; }
        trace.lo=ledger.first_uncovered(); trace.hi=ledger.last_uncovered();
        result.trace.push_back(std::move(trace));
    };
    auto commit = [&](std::optional<ClosedInterval> certificate) -> std::optional<ClosedInterval> {
        if (!certificate) return std::nullopt;
        certificate->lo = std::max(lower,certificate->lo);
        certificate->hi = std::min(upper,certificate->hi);
        if (certificate->empty()) return std::nullopt;
        if (ledger.add(*certificate)) counts.inc("certificates_committed");
        return certificate;
    };
    auto snapshot_finish = [&]() {
        result.lo=ledger.first_uncovered(); result.hi=ledger.last_uncovered();
        result.outstanding_reservations=reserved;
        counts.counts["work_units"]=spent;
        for (const auto& interval : ledger.intervals()) result.certificates.push_back({interval.first,interval.second});
        finalize_time();
    };

    // Preserve the cheap close-factor path before paying thread startup costs.
    if (options.algorithm.fermat_enabled && midpoint_band && allowance &&
        options.algorithm.initial_fermat_tests && global_time.remaining() > 0) {
        RegionFermatCursor seed(n,next_fermat,fermat_end,options.algorithm.native64);
        const auto quota=std::min(allowance,options.algorithm.initial_fermat_tests);
        seed.advance(quota,global_time.remaining());
        spent=seed.candidates;
        result.seed_work_units=spent;
        counts.counts["reservation_high_water"]=spent;
        counts.inc("fermat_tests",spent);
        const auto certificate=commit(seed.excluded_factors(lower));
        record("seed","fermat",seed.factor != 0 ? "factor" : "completed",0,0,
               seed.start,seed.covered_through,certificate,spent);
        next_fermat=seed.covered_through+1;
        if (seed.factor != 0) {
            result.factor=seed.factor; result.status="factor"; result.phase="seed_fermat";
            snapshot_finish(); return result;
        }
    }
    if (ledger.complete()) {
        result.status=full_domain ? "prime" : "region_exhausted";
        result.phase="coverage_complete"; snapshot_finish(); return result;
    }
    if (global_time.remaining() == 0 || spent == allowance) {
        result.reason=global_time.remaining()==0 ? "global time budget" : "global represented-work budget";
        result.phase="budget"; snapshot_finish(); return result;
    }

    using regional_detail::Job;
    using regional_detail::JobResult;
    std::mutex mutex;
    std::condition_variable changed;
    std::atomic<bool> stop{false};
    std::exception_ptr failure;
    std::deque<ClosedInterval> retry_product, retry_fermat;
    std::uint64_t next_id=1;
    std::size_t active=0;
    std::size_t active_fermat=0, active_product=0;
    std::uint64_t finished_fermat=0, finished_product=0;
    double elapsed_fermat=0, elapsed_product=0;
    const std::size_t per_lane=std::max<std::size_t>(1,options.degree/8);
    const Z product_width=Z(30)*static_cast<unsigned long>(options.degree)*static_cast<unsigned long>(per_lane);
    std::vector<AverageStream> streams;
    if (options.algorithm.vf_probes) streams=make_average_streams(n).second;

    auto next_product_region = [&]() -> std::optional<ClosedInterval> {
        for (;;) {
            const bool retry=!retry_product.empty();
            ClosedInterval band;
            if (retry) { band=retry_product.front(); retry_product.pop_front(); }
            else band={next_product,upper};
            const auto begin=ledger.first_uncovered_in(band.lo,band.hi);
            if (!begin) { if (retry) continue; return std::nullopt; }
            Z end=std::min(band.hi,Z(*begin+product_width-1));
            const auto covered=ledger.intervals().upper_bound(*begin);
            if (covered!=ledger.intervals().end()) end=std::min(end,Z(covered->first-1));
            if (retry) {
                if (end<band.hi) retry_product.push_back({Z(end+1),band.hi});
            } else next_product=end+1;
            return ClosedInterval{*begin,end};
        }
    };
    auto next_fermat_region = [&]() -> std::optional<ClosedInterval> {
        if (!options.algorithm.fermat_enabled || !midpoint_band) return std::nullopt;
        const auto lo=ledger.first_uncovered(), hi=ledger.last_uncovered();
        if (!lo || !hi) return std::nullopt;
        const auto current=midpoint_interval_for_factors(n,*lo,*hi);
        if (!current) return std::nullopt;
        for (;;) {
            const bool retry=!retry_fermat.empty();
            ClosedInterval band;
            if (retry) { band=retry_fermat.front(); retry_fermat.pop_front(); }
            else band={next_fermat,fermat_end};
            band.lo=std::max(band.lo,current->lo); band.hi=std::min(band.hi,current->hi);
            if (band.empty()) { if (retry) continue; return std::nullopt; }
            const Z end=std::min(band.hi,Z(band.lo+z_from_u64(options.fermat_span)-1));
            if (retry) { if (end<band.hi) retry_fermat.push_back({Z(end+1),band.hi}); }
            else next_fermat=end+1;
            return ClosedInterval{band.lo,end};
        }
    };
    auto claim = [&](std::unique_lock<std::mutex>& lock) -> std::optional<Job> {
        for (;;) {
            if (stop.load(std::memory_order_relaxed)) return std::nullopt;
            if (global_time.remaining()==0) {
                result.reason="global time budget"; stop.store(true); changed.notify_all(); return std::nullopt;
            }
            if (ledger.complete()) {
                result.status=full_domain ? "prime" : "region_exhausted";
                result.phase="coverage_complete"; stop.store(true); changed.notify_all(); return std::nullopt;
            }
            if (spent > allowance || reserved > allowance-spent) throw std::logic_error("global lease overflow");
            const auto available=allowance-spent-reserved;
            if (!available) {
                if (!active) {
                    result.reason="global represented-work budget"; stop.store(true); changed.notify_all(); return std::nullopt;
                }
                changed.wait_for(lock,std::chrono::milliseconds(1)); continue;
            }
            // Match time allocation rather than job counts: GMP Fermat and
            // polynomial jobs have different costs at different input sizes.
            // Include estimated in-flight cost to avoid assigning every free
            // worker to the same branch while its first result is pending.
            const double average_fermat=finished_fermat ? elapsed_fermat/finished_fermat : .001;
            const double average_product=finished_product ? elapsed_product/finished_product : .001;
            const double assigned_fermat=elapsed_fermat+active_fermat*average_fermat;
            const double assigned_product=elapsed_product+active_product*average_product;
            bool fermat=options.algorithm.fermat_enabled &&
                assigned_fermat <= options.algorithm.fermat_ratio*assigned_product;
            auto band=fermat ? next_fermat_region() : next_product_region();
            if (!band) { fermat=!fermat; band=fermat ? next_fermat_region() : next_product_region(); }
            if (!band) {
                if (active) { changed.wait_for(lock,std::chrono::milliseconds(1)); continue; }
                // Missing coverage can never be promoted to a prime claim.
                result.reason="no runnable region with uncertified coverage";
                stop.store(true); changed.notify_all(); return std::nullopt;
            }
            Job job; job.id=next_id++; job.fermat=fermat; job.lo=band->lo; job.hi=band->hi;
            std::uint64_t needed;
            if (fermat) needed=z_to_u64(Z(job.hi-job.lo+1));
            else {
                const Z term_bound=(job.hi-job.lo+1)*8/30+8;
                needed=z_to_u64(term_bound);
                for (auto& stream : streams) {
                    if (stream.cursor+static_cast<unsigned long>(options.degree)<=stream.count) {
                        job.auxiliary.push_back({Z(stream.start+5*stream.cursor),5,options.degree,"vf"});
                        stream.cursor+=static_cast<unsigned long>(options.degree);
                    }
                }
                needed+=static_cast<std::uint64_t>(job.auxiliary.size())*options.degree;
            }
            job.quota=std::min(needed,available);
            reserved+=job.quota; ++active;
            if (fermat) ++active_fermat; else ++active_product;
            counts.inc("jobs_started"); counts.inc(fermat ? "fermat_jobs" : "product_jobs");
            counts.counts["reservation_high_water"]=std::max(counts.get("reservation_high_water"),spent+reserved);
            counts.counts["peak_inflight_work"]=std::max(counts.get("peak_inflight_work"),reserved);
            counts.counts["active_workers_high_water"]=std::max<std::uint64_t>(counts.get("active_workers_high_water"),active);
            return job;
        }
    };

    auto worker = [&](std::size_t worker_id) {
        try {
            ProductRegionWorker products(n,options.algorithm.shared_polynomial);
            for (;;) {
                std::optional<Job> current;
                { std::unique_lock<std::mutex> lock(mutex); current=claim(lock); }
                if (!current) return;
                const Job& job=*current;
                JobResult outcome; outcome.covered_coordinate=job.lo-1;
                const auto job_started=Clock::now();
                try {
                    if (!stop.load(std::memory_order_relaxed)) {
                        const double seconds=std::min(options.job_seconds,std::max(0.0,global_time.remaining()));
                        if (job.fermat) {
                            RegionFermatCursor cursor(n,job.lo,job.hi,options.algorithm.native64);
                            const auto began=Clock::now();
                            while (!cursor.completed && cursor.factor==0 && cursor.candidates<job.quota &&
                                   !stop.load(std::memory_order_relaxed)) {
                                const double left=std::min(global_time.remaining(),seconds-seconds_since(began));
                                if (left<=0) break;
                                const auto chunk=std::min<std::uint64_t>(16384,job.quota-cursor.candidates);
                                const auto used=cursor.advance(chunk,left);
                                if (!used && !cursor.completed && cursor.factor==0) break;
                            }
                            outcome.factor=cursor.factor;
                            outcome.work=cursor.candidates;
                            outcome.counters.inc("fermat_tests",outcome.work);
                            outcome.certificate=cursor.excluded_factors(lower);
                            outcome.covered_coordinate=cursor.covered_through;
                            outcome.status=cursor.factor!=0 ? "factor" : cursor.completed ? "completed" :
                                stop.load(std::memory_order_relaxed) ? "cancelled" : "budget";
                        } else {
                            auto r=products.run(job.lo,job.hi,options.degree,job.auxiliary,seconds,job.quota,&stop);
                            outcome.factor=r.factor; outcome.work=r.work_units;
                            outcome.counters=std::move(r.counters);
                            outcome.covered_coordinate=r.covered_hi;
                            if (r.covered_hi>=r.covered_lo) outcome.certificate=ClosedInterval{r.covered_lo,r.covered_hi};
                            outcome.status=r.status;
                        }
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(mutex);
                    // An unexpected error invalidates this run. Reserve its
                    // full possible cost until all workers have been joined.
                    reserved-=job.quota; spent+=job.quota; --active;
                    if (job.fermat) --active_fermat; else --active_product;
                    if (!failure) failure=std::current_exception();
                    stop.store(true); changed.notify_all(); return;
                }
                outcome.elapsed_seconds=seconds_since(job_started);
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (outcome.work>job.quota) throw std::logic_error("worker exceeded reserved quota");
                    reserved-=job.quota; spent+=outcome.work; --active;
                    if (job.fermat) {
                        --active_fermat; ++finished_fermat; elapsed_fermat+=outcome.elapsed_seconds;
                        counts.add_seconds("fermat_job_seconds",outcome.elapsed_seconds);
                    } else {
                        --active_product; ++finished_product; elapsed_product+=outcome.elapsed_seconds;
                        counts.add_seconds("product_job_seconds",outcome.elapsed_seconds);
                    }
                    regional_detail::merge_counters(counts,outcome.counters);
                    counts.inc("jobs_finished");
                    if (outcome.status=="completed") counts.inc("jobs_completed");
                    else if (outcome.status=="cancelled") counts.inc("jobs_cancelled");
                    else if (outcome.status=="factor") counts.inc("factor_jobs");
                    else counts.inc("jobs_budget");
                    const auto certificate=commit(outcome.certificate);
                    if (outcome.factor!=0) {
                        if (outcome.factor<=1 || outcome.factor>=n || !mpz_divisible_p(n.get_mpz_t(),outcome.factor.get_mpz_t()))
                            throw std::logic_error("invalid factor from regional worker");
                        if (result.factor==0) {
                            result.factor=outcome.factor; result.status="factor";
                            result.reason.clear();
                            result.phase=job.fermat ? "parallel_fermat" : "parallel_product";
                        }
                        stop.store(true);
                    }
                    record("job_result",job.fermat ? "fermat" : "product",outcome.status,
                           job.id,worker_id,job.lo,job.hi,certificate,outcome.work);
                    if (!stop.load(std::memory_order_relaxed) && outcome.covered_coordinate<job.hi) {
                        auto& retry=job.fermat ? retry_fermat : retry_product;
                        retry.push_back({Z(std::max(job.lo,Z(outcome.covered_coordinate+1))),job.hi});
                    }
                    changed.notify_all();
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(mutex);
            if (!failure) failure=std::current_exception();
            stop.store(true); changed.notify_all();
        }
    };
    std::vector<std::thread> pool;
    try {
        pool.reserve(options.threads);
        for (std::size_t i=0;i<options.threads;++i) pool.emplace_back(worker,i);
    } catch (...) {
        stop.store(true); changed.notify_all();
        for (auto& thread : pool) if (thread.joinable()) thread.join();
        throw;
    }
    for (auto& thread : pool) thread.join();
    if (failure) std::rethrow_exception(failure);
    if (reserved || active) throw std::logic_error("unfinished work reservations after joining workers");
    if (result.phase.empty()) result.phase=result.status=="budget" ? "budget" : "coverage_complete";
    snapshot_finish();
    return result;
}

inline ParallelResult factor_parallel(const Z& n, const ParallelOptions& options = {}) {
    options.validate(); Budget global(options.algorithm.timeout);
    return parallel_find_factor(n,options,global,options.algorithm.max_steps);
}

struct ParallelFullResult {
    Z input;
    std::vector<Z> prime_factors, probable_prime_factors, unresolved;
    std::vector<ParallelResult> runs;
    std::uint64_t work_units=0;
    double elapsed_seconds=0, process_cpu_seconds=0;
};

inline ParallelFullResult factor_parallel_full(const Z& n, const ParallelOptions& options = {}) {
    options.validate();
    if (n<2) throw std::invalid_argument("N must be >= 2");
    if (options.factor_lo || options.factor_hi) throw std::invalid_argument("--full cannot be combined with a restricted factor band");
    Budget global(options.algorithm.timeout);
    const auto cpu_started=std::clock();
    ParallelFullResult result; result.input=n;
    auto work_remaining=options.algorithm.max_steps;
    std::vector<Z> pending{n};
    while (!pending.empty()) {
        if (global.remaining()==0) {
            result.unresolved.insert(result.unresolved.end(),pending.begin(),pending.end()); break;
        }
        Z current=std::move(pending.back()); pending.pop_back();
        auto run=parallel_find_factor(current,options,global,work_remaining);
        const auto work=run.counters.get("work_units");
        if (work_remaining) *work_remaining-=work;
        result.work_units+=work;
        if (run.factor!=0) { pending.push_back(run.factor); pending.push_back(current/run.factor); }
        else if (run.status=="prime") result.prime_factors.push_back(current);
        else if (run.status=="probable_prime") result.probable_prime_factors.push_back(current);
        else result.unresolved.push_back(current);
        result.runs.push_back(std::move(run));
    }
    for (auto* values : {&result.prime_factors,&result.probable_prime_factors,&result.unresolved}) std::sort(values->begin(),values->end());
    result.elapsed_seconds=global.elapsed();
    result.process_cpu_seconds=double(std::clock()-cpu_started)/CLOCKS_PER_SEC;
    return result;
}

} // namespace vf
