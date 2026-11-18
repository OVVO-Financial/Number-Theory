#pragma once

#include "gmp_arithmetic.hpp"
#include "gmp_products.hpp"
#include <array>
#include <memory>
#include <utility>

namespace vf {

struct Options {
    double timeout = 10.0;
    std::optional<std::uint64_t> max_steps;
    bool vf_probes = true, feedback = true, fermat_enabled = true;
    bool shared_polynomial = true, native64 = true, trace = false;
    bool probable_prime_check = true;
    std::size_t initial_degree = 16, max_degree = 512;
    double fermat_ratio = 2.0;
    std::uint64_t initial_fermat_tests = 64;
    void validate() const {
        if (!std::isfinite(timeout) || timeout < 0)
            throw std::invalid_argument("timeout must be finite and nonnegative");
        if (!initial_degree || max_degree < initial_degree || max_degree > 65536)
            throw std::invalid_argument("degrees must satisfy 1 <= initial <= maximum <= 65536");
        if (!std::isfinite(fermat_ratio) || fermat_ratio <= 0)
            throw std::invalid_argument("Fermat ratio must be finite and positive");
    }
};

struct TraceEvent {
    std::string event;
    std::map<std::string, std::optional<Z>> fields;
};

struct Result {
    Z n, factor = 0, lo = 7, hi = 0;
    std::optional<Z> fermat_next, fermat_covered;
    std::string status = "budget", phase, reason, fermat_backend = "disabled";
    double elapsed_seconds = 0;
    Counters counters;
    std::vector<TraceEvent> trace;
    std::vector<std::size_t> vf_levels;
    Options options;
};

struct AverageStream { Z start, count, cursor = 0; };

// Exact port of published_vf_parameters / average_streams. Selection depends
// only on N and the seed ceil(sqrt(N)); no factor information is supplied.
inline std::pair<std::size_t, std::vector<AverageStream>> make_average_streams(const Z& n) {
    const Z seed = ceil_sqrt(n);
    Z max_imaginary = floor_div(Z(n - 9), Z(6));
    if (max_imaginary < 0) max_imaginary = 0;
    const std::size_t average_count = max_imaginary == 0 ? 1 : mpz_sizeinbase(max_imaginary.get_mpz_t(), 2);
    std::size_t cut = 0;
    for (std::size_t value = average_count - 1; value; value >>= 1) ++cut;
    cut = std::max<std::size_t>(2, cut);
    std::size_t best_level = 1;
    Z best = -1;
    const Z difference = seed - n;
    for (std::size_t k = 1; k <= cut; ++k) {
        Z remainder;
        mpz_fdiv_r_2exp(remainder.get_mpz_t(), difference.get_mpz_t(), k);
        mpz_mul_2exp(remainder.get_mpz_t(), remainder.get_mpz_t(), cut-k);
        if (remainder > best) { best = remainder; best_level = k; }
    }
    const std::size_t level = std::max<std::size_t>(2, best_level);
    Z denominator = 1;
    mpz_mul_2exp(denominator.get_mpz_t(), denominator.get_mpz_t(), level);
    const Z numerator = (denominator - 1) * n + seed;
    Z low, remainder;
    mpz_fdiv_qr(low.get_mpz_t(), remainder.get_mpz_t(), numerator.get_mpz_t(), denominator.get_mpz_t());
    const Z high = low + (remainder != 0 ? 1 : 0);
    const Z total = 2 * floor_sqrt(n) + 2;
    Z ascending = floor_div(Z(n - 3 - high), Z(5)) + 1;
    if (ascending < 0) ascending = 0;
    if (ascending > total) ascending = total;
    const auto digit = mpz_fdiv_ui(n.get_mpz_t(), 10);
    const unsigned offsets = digit == 1 || digit == 9 ? 3 : 2;
    std::vector<AverageStream> streams;
    for (unsigned offset = 0; offset < offsets; ++offset)
        streams.push_back({Z(n - low + offset), total, Z(0)});
    for (unsigned offset = 0; offset < offsets; ++offset)
        streams.push_back({Z(high + offset), ascending, Z(0)});
    return {level, std::move(streams)};
}

inline Result find_factor(const Z& n, const Options& options, Budget& budget) {
    options.validate();
    if (n < 2) throw std::invalid_argument("N must be an integer >= 2");
    const auto started = Clock::now();
    const auto work_before = budget.used();
    Result result;
    result.n = n;
    result.options = options;
    auto& counts = result.counters;
    for (const char* key : {"scalar_tests", "fermat_tests", "represented_candidates",
         "divisor_candidates", "vf_candidates", "gcd_calls", "rounds", "fermat_chunks",
         "upper_bound_updates", "lower_bound_updates", "clipped_strip_positions"})
        counts.counts[key] = 0;
    counts.seconds["fermat_seconds"] = 0;
    counts.seconds["product_seconds"] = 0;
    counts.integers["fermat_bound_reduction"] = 0;
    Z lo = 7, hi = 0, square_root = 0;
    std::unique_ptr<GmpFermatCursor> cursor;
    std::map<std::size_t, std::unique_ptr<NormalizedProductEngine>> engines;
    std::vector<AverageStream> streams;
    bool streams_initialized = false;

    auto record = [&](const std::string& event, std::map<std::string, Z> extras = {}) {
        if (!options.trace) return;
        TraceEvent entry;
        entry.event = event;
        entry.fields["lo"] = lo; entry.fields["hi"] = hi;
        entry.fields["fermat_next"] = cursor ? std::optional<Z>(cursor->next_midpoint) : std::nullopt;
        entry.fields["fermat_covered"] = cursor ? std::optional<Z>(cursor->covered_through) : std::nullopt;
        entry.fields["fermat_cap"] = lo <= square_root ? std::optional<Z>(fermat_midpoint_cap(n,lo)) : std::nullopt;
        entry.fields["work_units"] = Z(std::to_string(budget.used()-work_before));
        for (const auto& item : extras) entry.fields[item.first] = item.second;
        result.trace.push_back(std::move(entry));
    };
    auto finish = [&](const Z& factor, const std::string& status, const std::string& phase,
                      const std::string& reason = "") -> Result {
        if (factor != 0 && (factor <= 1 || factor >= n || !mpz_divisible_p(n.get_mpz_t(),factor.get_mpz_t())))
            throw std::logic_error("invalid factor produced");
        result.factor = factor;
        result.status = factor != 0 ? "factor" : status;
        result.phase = phase; result.reason = reason;
        result.lo = lo; result.hi = hi;
        if (cursor) {
            result.fermat_next = cursor->next_midpoint;
            result.fermat_covered = cursor->covered_through;
        }
        counts.counts["work_units"] = budget.used()-work_before;
        result.elapsed_seconds = seconds_since(started);
        return std::move(result);
    };
    auto update_upper = [&]() {
        if (options.feedback && cursor) {
            const Z updated = factor_upper_bound(n,cursor->covered_through);
            if (updated < hi) {
                hi = updated;
                counts.inc("upper_bound_updates");
                counts.integers["fermat_bound_reduction"] = square_root-hi;
            }
        }
    };
    auto advance_fermat = [&](std::uint64_t maximum, double seconds) -> Z {
        budget.reserve(0);
        maximum = std::min(maximum,budget.available_work());
        if (!maximum) throw BudgetExceeded("represented-work budget");
        const double remaining = budget.remaining();
        if (remaining >= 0) seconds = seconds < 0 ? remaining : std::min(seconds,remaining);
        const Z cap = fermat_midpoint_cap(n, options.feedback ? lo : Z(7));
        const auto tick = Clock::now();
        const std::uint64_t used = cursor->advance(maximum,seconds,cap);
        counts.add_seconds("fermat_seconds",seconds_since(tick));
        budget.charge(used);
        counts.inc("fermat_tests",used); counts.inc("fermat_chunks");
        if (cursor->factor != 0) {
            record("fermat_factor",{{"factor",cursor->factor}});
            return cursor->factor;
        }
        update_upper();
        record("fermat",{{"chunk_tests",Z(std::to_string(used))},{"applied_cap",cap}});
        budget.check();
        return Z(0);
    };

    // These shared checks do not consume represented search work. Their time
    // is part of the caller's single budget. Full recursion checks its deadline
    // before starting another cofactor, including another set of prechecks.
    for (unsigned p : {2U,3U,5U}) {
        if (n == p) return finish(Z(0),"prime","precheck");
        if (mpz_divisible_ui_p(n.get_mpz_t(),p)) return finish(Z(p),"factor","precheck");
    }
    if (mpz_perfect_square_p(n.get_mpz_t())) return finish(floor_sqrt(n),"factor","precheck");
    const std::string primality = prime_status(n,options.probable_prime_check);
    if (primality == "prime" || primality == "probable_prime")
        return finish(Z(0),primality,"precheck");
    square_root = floor_sqrt(n);
    hi = square_root;

    try {
        budget.reserve(1);
        if (options.fermat_enabled) {
            cursor = std::make_unique<GmpFermatCursor>(n,options.native64);
            result.fermat_backend = cursor->backend;
            update_upper();
        }
        record("init");
        if (cursor && options.initial_fermat_tests) {
            const Z found = advance_fermat(options.initial_fermat_tests,budget.remaining());
            if (found != 0) return finish(found,"factor","fermat");
        }
        if (lo > hi) return finish(Z(0),"prime","frontiers_met");
        const unsigned prefix_end = hi < 119 ? static_cast<unsigned>(hi.get_ui()) : 119;
        constexpr std::array<unsigned,8> residues = {1,7,11,13,17,19,23,29};
        for (unsigned d = 7; d <= prefix_end; ++d) {
            if (std::find(residues.begin(),residues.end(),d%30) == residues.end()) continue;
            budget.reserve(1); budget.charge(1);
            counts.inc("scalar_tests"); counts.inc("divisor_candidates"); counts.inc("gcd_calls");
            const unsigned long g = mpz_gcd_ui(nullptr,n.get_mpz_t(),d);
            if (g > 1 && Z(g) < n) {
                record("prefix_factor",{{"factor",Z(g)}});
                return finish(Z(g),"factor","divisor_prefix");
            }
        }
        lo = prefix_end+1;
        counts.inc("lower_bound_updates");
        record("divisor_prefix");
        std::size_t degree = options.initial_degree;

        while (lo <= hi) {
            budget.reserve(1);
            counts.inc("rounds");
            const auto tick = Clock::now();
            bool timer_recorded = false;
            auto record_time = [&]() {
                if (!timer_recorded) {
                    counts.add_seconds("product_seconds",seconds_since(tick));
                    timer_recorded = true;
                }
            };
            try {
                Z available = (hi-lo+1)*8/30;
                if (available < 1) available = 1;
                const Z degree_squared = Z(static_cast<unsigned long>(degree))*static_cast<unsigned long>(degree);
                const std::size_t m = available >= degree_squared ? degree :
                    std::max<std::size_t>(1,floor_sqrt(available).get_ui());
                const std::size_t blocks_per_lane = std::max<std::size_t>(1,m/8);
                const Z base = (lo/30)*30;
                const Z nominal_end = base + Z(30)*static_cast<unsigned long>(m)*static_cast<unsigned long>(blocks_per_lane)-1;
                const Z end = std::min(hi,nominal_end);
                std::vector<Request> requests;
                std::vector<Z> tails;
                requests.reserve(8*blocks_per_lane+6);
                for (unsigned residue : residues) {
                    Z first = base+residue;
                    if (first < lo) first += 30;
                    std::size_t count = 0;
                    if (end >= first) {
                        const Z count_z = (end-first)/30+1;
                        count = count_z.get_ui(); // bounded by the configured strip width
                    }
                    const std::size_t full = count/m, tail = count%m;
                    for (std::size_t j=0; j<full; ++j)
                        requests.push_back({Z(first+Z(30)*static_cast<unsigned long>(m)*static_cast<unsigned long>(j)),30,m,"divisor"});
                    for (std::size_t j=0; j<tail; ++j)
                        tails.push_back(first+Z(30)*static_cast<unsigned long>(m)*static_cast<unsigned long>(full)+30*static_cast<unsigned long>(j));
                }
                const std::size_t divisor_blocks = requests.size();
                if (options.vf_probes) {
                    if (!streams_initialized) {
                        auto made = make_average_streams(n);
                        result.vf_levels.push_back(made.first);
                        streams = std::move(made.second);
                        streams_initialized = true;
                    }
                    for (auto& stream : streams) {
                        if (stream.cursor+static_cast<unsigned long>(m) <= stream.count) {
                            requests.push_back({Z(stream.start+5*stream.cursor),5,m,"vf"});
                            stream.cursor += static_cast<unsigned long>(m);
                        }
                    }
                }
                const std::size_t vf_blocks = requests.size()-divisor_blocks;
                budget.reserve(static_cast<std::uint64_t>(requests.size())*m+tails.size());
                const std::map<std::string,Z> plan = {
                    {"degree",Z(static_cast<unsigned long>(m))},{"strip_end",end},{"nominal_strip_end",nominal_end},
                    {"divisor_blocks",Z(static_cast<unsigned long>(divisor_blocks))},
                    {"vf_blocks",Z(static_cast<unsigned long>(vf_blocks))},{"tail_count",Z(static_cast<unsigned long>(tails.size()))}
                };
                record("product_plan",plan);
                if (!requests.empty()) {
                    auto& engine = engines[m];
                    if (!engine) engine = std::make_unique<NormalizedProductEngine>(n,m,&budget,&counts,options.shared_polynomial);
                    std::vector<Z> starts; std::vector<long> steps;
                    starts.reserve(requests.size()); steps.reserve(requests.size());
                    for (const auto& request : requests) { starts.push_back(request.start); steps.push_back(request.step); }
                    const auto represented_before = counts.get("represented_candidates");
                    std::vector<Z> values;
                    try { values = engine->products(starts,steps); }
                    catch (...) {
                        budget.charge(counts.get("represented_candidates")-represented_before);
                        throw;
                    }
                    budget.charge(counts.get("represented_candidates")-represented_before);
                    counts.inc("divisor_candidates",static_cast<std::uint64_t>(divisor_blocks)*m);
                    counts.inc("vf_candidates",static_cast<std::uint64_t>(vf_blocks)*m);
                    const Z found = recover_mixed_batch(n,requests,values,&budget,&counts);
                    if (found != 0) {
                        record_time();
                        auto extra = plan; extra["factor"] = found;
                        record("product_factor",extra);
                        return finish(found,"factor",vf_blocks ? "shared_product" : "divisor_product");
                    }
                }
                for (const Z& d : tails) {
                    budget.reserve(1); budget.charge(1);
                    counts.inc("scalar_tests"); counts.inc("divisor_candidates"); counts.inc("gcd_calls");
                    Z found;
                    mpz_gcd(found.get_mpz_t(),d.get_mpz_t(),n.get_mpz_t());
                    if (found > 1 && found < n) {
                        record_time(); record("tail_factor",{{"factor",found}});
                        return finish(found,"factor","divisor_tail");
                    }
                }
                // Certification is atomic: all residue lanes, tails and GCD-N
                // recoveries have completed before excluding this divisor strip.
                lo = end+1;
                counts.inc("lower_bound_updates");
                if (end < nominal_end) {
                    const Z clipped = (nominal_end-end)*8/30;
                    counts.inc("clipped_strip_positions",clipped.get_ui());
                }
                const double elapsed = seconds_since(tick);
                record_time(); record("product",plan);
                budget.check();
                if (lo > hi) break;
                if (cursor) {
                    const double allocation = std::min(.02,std::max(.00005,options.fermat_ratio*elapsed));
                    const Z found = advance_fermat(std::numeric_limits<std::uint64_t>::max(),allocation);
                    if (found != 0) return finish(found,"factor","fermat");
                }
                degree = std::min(2*degree,options.max_degree);
            } catch (...) { record_time(); throw; }
        }
        return finish(Z(0),"prime","frontiers_met");
    } catch (const BudgetExceeded& error) {
        return finish(Z(0),"budget","budget",error.what());
    }
}

inline Result factor(const Z& n, const Options& options = Options{}) {
    options.validate();
    Budget budget(options.timeout,options.max_steps);
    return find_factor(n,options,budget);
}

struct FullResult {
    Z input;
    std::vector<Z> prime_factors, probable_prime_factors, unresolved;
    std::vector<Result> runs;
    double elapsed_seconds = 0;
    std::uint64_t work_units = 0;
};

inline FullResult factor_full(const Z& n, const Options& options = Options{}) {
    options.validate();
    if (n < 2) throw std::invalid_argument("N must be an integer >= 2");
    Budget budget(options.timeout,options.max_steps);
    FullResult result; result.input = n;
    std::vector<Z> pending{n};
    while (!pending.empty()) {
        // Unlike a succession of independent calls, recursive prechecks cannot
        // start new work after the overall full-factorization deadline expires.
        if (budget.remaining() == 0) {
            result.unresolved.insert(result.unresolved.end(),pending.begin(),pending.end());
            break;
        }
        Z current = std::move(pending.back()); pending.pop_back();
        Result run = find_factor(current,options,budget);
        if (run.factor != 0) {
            pending.push_back(run.factor); pending.push_back(current/run.factor);
        } else if (run.status == "prime") result.prime_factors.push_back(current);
        else if (run.status == "probable_prime") result.probable_prime_factors.push_back(current);
        else result.unresolved.push_back(current);
        result.runs.push_back(std::move(run));
    }
    for (auto* values : {&result.prime_factors,&result.probable_prime_factors,&result.unresolved})
        std::sort(values->begin(),values->end());
    result.elapsed_seconds = budget.elapsed();
    result.work_units = budget.used();
    return result;
}

} // namespace vf
