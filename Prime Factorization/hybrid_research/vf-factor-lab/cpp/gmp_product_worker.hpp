#pragma once

// Thread-confined divisor-region searches with reusable normalized polynomials.
// The coordinator assigns independent regions and reserves a global work quota.
// This worker never shares mutable GMP objects, budgets, or engine state.

#include "gmp_products.hpp"
#include <array>
#include <atomic>
#include <memory>

namespace vf {

struct ProductRegionResult {
    Z factor = 0;
    Z covered_lo = 0, covered_hi = -1;
    std::string status;
    Counters counters;
    double elapsed_seconds = 0;
    std::uint64_t work_units = 0;
    std::string reason;
};

class ProductRegionWorker {
    Z n_, square_root_;
    bool shared_;
    // Assigning fresh values to these members retains their addresses. Cached
    // engines can therefore safely refer to them across sequential jobs.
    Budget budget_;
    Counters counters_;
    std::map<std::size_t, std::unique_ptr<NormalizedProductEngine>> engines_;

    inline static constexpr std::array<unsigned,8> residues_ = {{1,7,11,13,17,19,23,29}};

    struct Strip {
        Z end;
        std::vector<Request> requests;
        std::vector<Z> tails;
        std::uint64_t divisor_terms = 0;
    };

    static bool stopped(const std::atomic<bool>* stop) {
        return stop && stop->load(std::memory_order_relaxed);
    }

    static Z next_wheel_candidate(const Z& start) {
        const unsigned remainder = mpz_fdiv_ui(start.get_mpz_t(),30);
        const auto it = std::lower_bound(residues_.begin(),residues_.end(),remainder);
        if (it == residues_.end()) return Z(start + 30 - remainder + residues_[0]);
        return Z(start + *it - remainder);
    }

    static std::size_t affordable_degree(std::size_t degree, std::uint64_t available) {
        // degree <= 65536, so this square is safe in uint64_t.
        if (available >= static_cast<std::uint64_t>(degree)*degree) return degree;
        if (!available) return 1;
        Z allowance(std::to_string(available)), root;
        mpz_sqrt(root.get_mpz_t(),allowance.get_mpz_t());
        return std::max<std::size_t>(1,root.get_ui());
    }

    static Strip plan_strip(const Z& lo, const Z& hi, std::size_t degree) {
        const std::size_t blocks_per_lane = std::max<std::size_t>(1,degree/8);
        const Z base = (lo/30)*30;
        const Z nominal_end = base+Z(30)*static_cast<unsigned long>(degree)*
                                      static_cast<unsigned long>(blocks_per_lane)-1;
        Strip strip;
        strip.end = std::min(hi,nominal_end);
        strip.requests.reserve(8*blocks_per_lane+6);
        for (unsigned residue : residues_) {
            Z first = base+residue;
            if (first < lo) first += 30;
            std::size_t count = 0;
            if (first <= strip.end) {
                const Z count_z = (strip.end-first)/30+1;
                // Every count is bounded by the configured strip width.
                count = count_z.get_ui();
            }
            strip.divisor_terms += count;
            const std::size_t full = count/degree, tail = count%degree;
            for (std::size_t j=0; j<full; ++j)
                strip.requests.push_back({Z(first+Z(30)*static_cast<unsigned long>(degree)*
                                                          static_cast<unsigned long>(j)),
                                          30,degree,"divisor"});
            for (std::size_t j=0; j<tail; ++j)
                strip.tails.push_back(first+Z(30)*static_cast<unsigned long>(degree)*
                                                static_cast<unsigned long>(full)+30*static_cast<unsigned long>(j));
        }
        return strip;
    }

public:
    explicit ProductRegionWorker(const Z& n, bool shared=true)
        : n_(n), shared_(shared), budget_(-1) {
        if (n_ < 2 || mpz_gcd_ui(nullptr,n_.get_mpz_t(),30) != 1)
            throw std::invalid_argument("product worker requires N >= 2 coprime to 30");
        mpz_sqrt(square_root_.get_mpz_t(),n_.get_mpz_t());
    }

    // Engine references depend on this object's stable address. Each worker
    // is owned by one thread and runs at most one job at a time.
    ProductRegionWorker(const ProductRegionWorker&) = delete;
    ProductRegionWorker& operator=(const ProductRegionWorker&) = delete;
    ProductRegionWorker(ProductRegionWorker&&) = delete;
    ProductRegionWorker& operator=(ProductRegionWorker&&) = delete;

    ProductRegionResult run(const Z& lo, const Z& hi, std::size_t degree,
                            const std::vector<Request>& auxiliary, double seconds,
                            std::uint64_t max_work,
                            const std::atomic<bool>* stop=nullptr) {
        if (lo < 7 || hi < lo || hi > square_root_)
            throw std::invalid_argument("product region must satisfy 7 <= lo <= hi <= floor(sqrt(N))");
        if (!degree || degree > 65536)
            throw std::invalid_argument("product degree must be between 1 and 65536");
        if (!std::isfinite(seconds) || (seconds < 0 && seconds != -1))
            throw std::invalid_argument("product timeout must be nonnegative, or -1 for unlimited");
        if (auxiliary.size() > std::numeric_limits<std::uint64_t>::max()/degree)
            throw std::overflow_error("auxiliary work count overflow");
        for (const auto& request : auxiliary) {
            if (request.count != degree)
                throw std::invalid_argument("auxiliary blocks must have the requested degree");
            Z gcd, step = request.step;
            mpz_gcd(gcd.get_mpz_t(),step.get_mpz_t(),n_.get_mpz_t());
            if (gcd != 1) throw std::invalid_argument("auxiliary steps must be coprime to N");
        }
        budget_ = Budget(seconds,max_work);
        counters_ = Counters{};
        for (const char* key : {"scalar_tests","represented_candidates","divisor_candidates",
             "vf_candidates","gcd_calls","completed_strips","scalar_fallbacks",
             "auxiliary_skipped_for_budget","engine_cache_hits"})
            counters_.counts[key] = 0;
        Z current = lo, covered = lo-1;
        bool first_strip = true;

        auto finish = [&](const std::string& status, const Z& factor=Z(0),
                          const std::string& reason="") {
            if (factor != 0 && (factor <= 1 || factor >= n_ ||
                !mpz_divisible_p(n_.get_mpz_t(),factor.get_mpz_t())))
                throw std::logic_error("product worker produced an invalid factor");
            ProductRegionResult result;
            result.factor = factor; result.covered_lo = lo; result.covered_hi = covered;
            result.status = status; result.reason = reason;
            result.work_units = budget_.used();
            counters_.counts["work_units"] = result.work_units;
            counters_.counts["cached_degree_engines"] = engines_.size();
            result.counters = counters_;
            result.elapsed_seconds = budget_.elapsed();
            return result;
        };

        try {
            while (current <= hi) {
                if (stopped(stop)) return finish("cancelled");
                budget_.check();
                const std::uint64_t allowance = budget_.available_work();
                const std::size_t m = affordable_degree(degree,allowance);
                Strip strip = plan_strip(current,hi,m);

                // Very small quotas cannot fit all eight wheel lanes of even
                // a degree-one strip. Test a contiguous scalar prefix instead
                // of returning an unchanged zero-progress job indefinitely.
                if (strip.divisor_terms > allowance) {
                    counters_.inc("scalar_fallbacks");
                    if (first_strip && !auxiliary.empty())
                        counters_.inc("auxiliary_skipped_for_budget",auxiliary.size());
                    first_strip = false;
                    while (current <= hi) {
                        if (stopped(stop)) return finish("cancelled");
                        budget_.check();
                        const Z candidate = next_wheel_candidate(current);
                        if (candidate > hi) { covered = hi; current = hi+1; break; }
                        // The gap consists only of multiples of 2, 3, or 5.
                        covered = candidate-1;
                        current = candidate;
                        budget_.reserve(1); budget_.charge(1);
                        counters_.inc("scalar_tests"); counters_.inc("divisor_candidates");
                        counters_.inc("gcd_calls");
                        Z factor;
                        mpz_gcd(factor.get_mpz_t(),candidate.get_mpz_t(),n_.get_mpz_t());
                        if (factor > 1 && factor < n_) return finish("factor",factor);
                        const Z next = next_wheel_candidate(Z(candidate+1));
                        covered = std::min(hi,Z(next-1));
                        current = covered+1;
                    }
                    continue;
                }

                const std::size_t divisor_blocks = strip.requests.size();
                if (first_strip && !auxiliary.empty()) {
                    const std::uint64_t auxiliary_work = static_cast<std::uint64_t>(auxiliary.size())*degree;
                    if (m == degree && auxiliary_work <= allowance-strip.divisor_terms)
                        strip.requests.insert(strip.requests.end(),auxiliary.begin(),auxiliary.end());
                    else counters_.inc("auxiliary_skipped_for_budget",auxiliary.size());
                }
                first_strip = false;
                const std::size_t auxiliary_blocks = strip.requests.size()-divisor_blocks;
                const std::uint64_t request_work = static_cast<std::uint64_t>(strip.requests.size())*m;
                budget_.reserve(request_work+strip.tails.size());
                if (!strip.requests.empty()) {
                    auto found_engine = engines_.find(m);
                    if (found_engine == engines_.end()) {
                        auto engine = std::make_unique<NormalizedProductEngine>(n_,m,&budget_,&counters_,shared_);
                        found_engine = engines_.emplace(m,std::move(engine)).first;
                    } else counters_.inc("engine_cache_hits");
                    counters_.counts["polynomial_degree"] =
                        std::max<std::uint64_t>(counters_.get("polynomial_degree"),m);
                    std::vector<Z> starts; std::vector<long> steps;
                    starts.reserve(strip.requests.size()); steps.reserve(strip.requests.size());
                    for (const auto& request : strip.requests) {
                        starts.push_back(request.start); steps.push_back(request.step);
                    }
                    const auto represented_before = counters_.get("represented_candidates");
                    std::vector<Z> values;
                    try { values = found_engine->second->products(starts,steps); }
                    catch (...) {
                        budget_.charge(counters_.get("represented_candidates")-represented_before);
                        throw;
                    }
                    budget_.charge(counters_.get("represented_candidates")-represented_before);
                    counters_.inc("divisor_candidates",static_cast<std::uint64_t>(divisor_blocks)*m);
                    counters_.inc("vf_candidates",static_cast<std::uint64_t>(auxiliary_blocks)*m);
                    if (stopped(stop)) return finish("cancelled");
                    const Z factor = recover_mixed_batch(n_,strip.requests,values,&budget_,&counters_);
                    if (factor != 0) return finish("factor",factor);
                }
                for (const Z& candidate : strip.tails) {
                    if (stopped(stop)) return finish("cancelled");
                    budget_.reserve(1); budget_.charge(1);
                    counters_.inc("scalar_tests"); counters_.inc("divisor_candidates");
                    counters_.inc("gcd_calls");
                    Z factor;
                    mpz_gcd(factor.get_mpz_t(),candidate.get_mpz_t(),n_.get_mpz_t());
                    if (factor > 1 && factor < n_) return finish("factor",factor);
                }
                // Only this assignment certifies a whole polynomial strip.
                // Earlier exceptions or cancellation leave that strip pending.
                covered = strip.end;
                current = strip.end+1;
                counters_.inc("completed_strips");
            }
            return finish("completed");
        } catch (const BudgetExceeded& error) {
            return finish("budget",Z(0),error.what());
        }
    }
};

} // namespace vf
