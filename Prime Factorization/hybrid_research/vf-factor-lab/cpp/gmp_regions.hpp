#pragma once

#include "gmp_arithmetic.hpp"
#include <map>
#include <optional>
#include <utility>

extern "C" {
unsigned fusion_region_abi_version();
void* fusion_fermat_create_at(std::uint64_t n, std::uint64_t begin, unsigned* error);
}

namespace vf {

struct ClosedInterval {
    Z lo = 0, hi = -1;
    bool empty() const { return lo > hi; }
    Z size() const { return empty() ? Z(0) : Z(hi-lo+1); }
};

// For odd N every actual divisor pair has an integer midpoint. A complete
// failed midpoint interval [A,B] therefore excludes exactly the integer d
// satisfying A <= (d+N/d)/2 <= B in 1 <= d <= floor(sqrt(N)).
// These bounds describe a local interval, NOT a certified prefix from sqrt(N).
inline std::optional<ClosedInterval> excluded_factor_interval(
        const Z& n, const Z& a, const Z& b, const Z& lower = Z(7)) {
    if (n < 3 || mpz_even_p(n.get_mpz_t()))
        throw std::invalid_argument("midpoint certificates require odd N >= 3");
    if (lower < 1) throw std::invalid_argument("factor-domain lower bound must be positive");
    const Z root = floor_sqrt(n);
    const Z start = std::max(a,ceil_sqrt(n));
    const Z endpoint = (n+1)/2;
    const Z end = std::min(b,endpoint);
    if (start > end || lower > root) return std::nullopt;
    const Z raw_lo = end-floor_sqrt(Z(end*end-n));
    const Z raw_hi = start-ceil_sqrt(Z(start*start-n));
    const Z lo = std::max(lower,raw_lo), hi = std::min(root,raw_hi);
    if (lo > hi) return std::nullopt;
    return ClosedInterval{lo,hi};
}

// Stronger certificate for ACTUAL factors of odd N. Their midpoints are
// integers, so excluding every integer midpoint A..B also excludes factors
// with A-1 < (d+N/d)/2 < B+1. The strict inequalities deliberately retain
// possible factor hits at A-1 and B+1. Adjacent completed midpoint windows
// consequently leave no artificial holes caused by nonintegral midpoints.
inline std::optional<ClosedInterval> excluded_integer_midpoint_interval(
        const Z& n, const Z& a, const Z& b, const Z& lower = Z(7)) {
    if (n < 3 || mpz_even_p(n.get_mpz_t()))
        throw std::invalid_argument("integer midpoint certificates require odd N >= 3");
    if (lower < 1) throw std::invalid_argument("factor-domain lower bound must be positive");
    const Z root = floor_sqrt(n), ceiling = ceil_sqrt(n);
    const Z start = std::max(a,ceiling), endpoint = (n+1)/2;
    const Z end = std::min(b,endpoint);
    if (start > end || lower > root) return std::nullopt;
    const Z after = end+1;
    const Z raw_lo = after-ceil_sqrt(Z(after*after-n))+1;
    Z raw_hi = root;
    if (start > ceiling) {
        const Z before = start-1;
        raw_hi = before-floor_sqrt(Z(before*before-n))-1;
    }
    const Z lo = std::max(lower,raw_lo), hi = std::min(root,raw_hi);
    if (lo > hi) return std::nullopt;
    return ClosedInterval{lo,hi};
}

// Every possible factor in the specified band has its midpoint in this
// interval. Only integer midpoints are relevant after the factor-2 precheck.
inline std::optional<ClosedInterval> midpoint_interval_for_factors(
        const Z& n, const Z& factor_lo, const Z& factor_hi) {
    if (n < 3 || mpz_even_p(n.get_mpz_t()))
        throw std::invalid_argument("factor-to-midpoint ranges require odd N >= 3");
    const Z lo = std::max(Z(1),factor_lo);
    const Z hi = std::min(floor_sqrt(n),factor_hi);
    if (lo > hi) return std::nullopt;
    Z start;
    const Z numerator = hi*hi+n, denominator = 2*hi;
    mpz_cdiv_q(start.get_mpz_t(),numerator.get_mpz_t(),denominator.get_mpz_t());
    const Z end = floor_div(Z(lo*lo+n),Z(2*lo));
    if (start > end) return std::nullopt;
    return ClosedInterval{start,end};
}

// This ledger is owned by the coordinator, or protected by its mutex. It is
// deliberately not internally synchronized. Only completed certificates are
// inserted. Out-of-order intervals do not move a frontier across any hole.
class CoverageLedger {
    Z lower_, upper_, length_ = 0;
    std::map<Z,Z> covered_;
public:
    CoverageLedger(const Z& lower, const Z& upper) : lower_(lower), upper_(upper) {}
    const Z& lower() const { return lower_; }
    const Z& upper() const { return upper_; }
    const Z& covered_length() const { return length_; }
    const Z& total_certified_length() const { return length_; }
    const std::map<Z,Z>& intervals() const { return covered_; }
    bool domain_empty() const { return lower_ > upper_; }

    bool covers(const Z& lo, const Z& hi) const {
        if (lo > hi) return true;
        if (domain_empty() || lo < lower_ || hi > upper_) return false;
        auto it = covered_.upper_bound(lo);
        if (it == covered_.begin()) return false;
        --it;
        return it->first <= lo && it->second >= hi;
    }

    bool add(const Z& requested_lo, const Z& requested_hi) {
        Z lo = std::max(lower_,requested_lo), hi = std::min(upper_,requested_hi);
        if (lo > hi || covers(lo,hi)) return false;
        auto it = covered_.lower_bound(lo);
        if (it != covered_.begin()) {
            auto previous = std::prev(it);
            if (previous->second+1 >= lo) it = previous;
        }
        while (it != covered_.end() && it->first <= hi+1) {
            lo = std::min(lo,it->first);
            hi = std::max(hi,it->second);
            length_ -= it->second-it->first+1;
            it = covered_.erase(it);
        }
        length_ += hi-lo+1;
        covered_.emplace(std::move(lo),std::move(hi));
        return true;
    }

    bool add(const ClosedInterval& interval) { return add(interval.lo,interval.hi); }

    std::optional<Z> next_uncovered(const Z& requested) const {
        Z point = std::max(lower_,requested);
        if (point > upper_) return std::nullopt;
        auto it = covered_.upper_bound(point);
        if (it != covered_.begin()) {
            --it;
            if (it->second >= point) point = it->second+1;
        }
        if (point > upper_) return std::nullopt;
        return point;
    }

    std::optional<Z> previous_uncovered(const Z& requested) const {
        Z point = std::min(upper_,requested);
        if (point < lower_) return std::nullopt;
        auto it = covered_.upper_bound(point);
        if (it != covered_.begin()) {
            --it;
            if (it->second >= point) point = it->first-1;
        }
        if (point < lower_) return std::nullopt;
        return point;
    }

    std::optional<Z> first_uncovered() const { return next_uncovered(lower_); }
    std::optional<Z> last_uncovered() const { return previous_uncovered(upper_); }
    bool complete() const { return !first_uncovered().has_value(); }

    std::optional<Z> first_uncovered_in(const Z& lo, const Z& hi) const {
        const auto next = next_uncovered(lo);
        return next && *next <= hi ? next : std::nullopt;
    }

    std::vector<ClosedInterval> gaps() const {
        std::vector<ClosedInterval> result;
        if (domain_empty()) return result;
        Z next = lower_;
        for (const auto& interval : covered_) {
            if (interval.first > next) result.push_back({next,Z(interval.first-1)});
            next = interval.second+1;
        }
        if (next <= upper_) result.push_back({next,upper_});
        return result;
    }
};

// A job-local cursor. covered_through certifies ONLY [start,covered_through].
// In particular, passing it to the old global factor_upper_bound is invalid.
// Use excluded_factors() or excluded_integer_midpoint_interval() instead.
class RegionFermatCursor {
public:
    Z n, start, end, next_midpoint, covered_through, factor = 0, half_gap = 0;
    Z first_midpoint = 0, last_midpoint = 0;
    std::uint64_t candidates = 0;
    bool completed = false;
    std::string status = "ready", backend = "gmp";
    double elapsed_seconds = 0, total_kernel_seconds = 0;

    RegionFermatCursor(const Z& number, const Z& requested_start,
                      const Z& requested_end, bool native64 = true) : n(number) {
        if (n < 3 || mpz_even_p(n.get_mpz_t()))
            throw std::invalid_argument("regional Fermat cursor requires odd N >= 3");
        const auto began = Clock::now();
        const Z root = ceil_sqrt(n), endpoint = (n+1)/2;
        start = std::max(root,requested_start);
        end = std::min(endpoint,requested_end);
        next_midpoint = start;
        covered_through = start-1;
        const bool use_native = native64 && mpz_sizeinbase(n.get_mpz_t(),2) <= 64;
        backend = use_native ? "native64" : "gmp";
        if (start > end) {
            completed = true; status = "completed";
            total_kernel_seconds = seconds_since(began);
            return;
        }
        if (use_native) {
            if (fusion_abi_version() != 1 || fusion_region_abi_version() != 1 ||
                fusion_snapshot_size() != sizeof(FusionSnapshot))
                throw std::runtime_error("native regional Fermat ABI mismatch");
            unsigned error = 0;
            native_ = fusion_fermat_create_at(z_to_u64(n),z_to_u64(start),&error);
            if (!native_) throw std::runtime_error("native regional Fermat initialization failed");
            FusionSnapshot snapshot{};
            if (fusion_fermat_snapshot(native_,&snapshot)) {
                fusion_fermat_destroy(native_); native_ = nullptr;
                throw std::runtime_error("native regional Fermat snapshot failed");
            }
            read_native(snapshot);
            return;
        }
        std::array<bool,20> qr20{};
        std::array<bool,64> qr64{};
        for (unsigned b=0;b<20;++b) qr20[(b*b)%20] = true;
        for (unsigned b=0;b<64;++b) qr64[(b*b)%64] = true;
        const unsigned n20 = mpz_fdiv_ui(n.get_mpz_t(),20);
        const unsigned n64 = mpz_fdiv_ui(n.get_mpz_t(),64);
        for (unsigned r=0;r<160;++r)
            if (qr20[(r*r+20-n20)%20] && qr64[(r*r+64-n64)%64]) residues_.push_back(r);
        if (residues_.empty()) throw std::logic_error("empty regional Fermat wheel");
        for (std::size_t i=0;i<residues_.size();++i)
            gaps_.push_back(i+1<residues_.size() ? residues_[i+1]-residues_[i] : 160+residues_[0]-residues_[i]);
        const unsigned remainder = mpz_fdiv_ui(start.get_mpz_t(),160);
        const auto first = std::lower_bound(residues_.begin(),residues_.end(),remainder);
        if (first == residues_.end()) { index_ = 0; next_midpoint += 160-remainder+residues_[0]; }
        else { index_ = static_cast<std::size_t>(first-residues_.begin()); next_midpoint += *first-remainder; }
        covered_through = std::min(end,Z(next_midpoint-1));
        difference_ = next_midpoint*next_midpoint-n;
        completed = next_midpoint > end;
        if (completed) status = "completed";
        total_kernel_seconds = seconds_since(began);
    }

    ~RegionFermatCursor() { if (native_) fusion_fermat_destroy(native_); }
    RegionFermatCursor(const RegionFermatCursor&) = delete;
    RegionFermatCursor& operator=(const RegionFermatCursor&) = delete;

    std::optional<ClosedInterval> certified_interval() const {
        if (start > end || covered_through < start) return std::nullopt;
        return ClosedInterval{start,std::min(end,covered_through)};
    }

    std::optional<ClosedInterval> excluded_factors(const Z& lower = Z(7)) const {
        const auto interval = certified_interval();
        if (!interval) return std::nullopt;
        return excluded_integer_midpoint_interval(n,interval->lo,interval->hi,lower);
    }

    // A finite candidate quota allows the coordinator to reserve shared work
    // before entering this cursor. A short time quantum bounds cancellation
    // latency. Negative seconds disables this call's time limit.
    std::uint64_t advance(std::uint64_t maximum, double seconds) {
        if (!std::isfinite(seconds)) throw std::invalid_argument("regional timeout must be finite");
        maximum = std::min(maximum,std::numeric_limits<std::uint64_t>::max()-candidates);
        if (factor != 0 || completed) return 0;
        if (native_) {
            FusionSnapshot snapshot{};
            if (fusion_fermat_advance(native_,maximum,seconds,z_to_u64(end),&snapshot))
                throw std::runtime_error("native regional Fermat advance failed");
            read_native(snapshot);
            return snapshot.delta_candidates;
        }
        const auto began = Clock::now();
        const auto before = candidates;
        auto finish = [&](const char* why) {
            status = why;
            elapsed_seconds = seconds_since(began);
            total_kernel_seconds += elapsed_seconds;
            return candidates-before;
        };
        if (!maximum) return finish("candidate_budget");
        if (seconds == 0) return finish("time_budget");
        while (next_midpoint <= end) {
            const auto used = candidates-before;
            if (used >= maximum) return finish("candidate_budget");
            if ((used&255U)==0 && seconds >= 0 && seconds_since(began) >= seconds)
                return finish("time_budget");
            ++candidates;
            if (candidates == 1) first_midpoint = next_midpoint;
            last_midpoint = next_midpoint;
            if (mpz_perfect_square_p(difference_.get_mpz_t())) {
                mpz_sqrt(root_buffer_.get_mpz_t(),difference_.get_mpz_t());
                const Z divisor = next_midpoint-root_buffer_;
                if (divisor > 1 && divisor < n) {
                    factor = divisor; half_gap = root_buffer_;
                    return finish("factor");
                }
                if (divisor == 1) {
                    covered_through = end;
                    ++next_midpoint;
                    completed = true;
                    return finish("completed");
                }
            }
            const unsigned gap = gaps_[index_];
            mpz_mul_ui(increment_.get_mpz_t(),next_midpoint.get_mpz_t(),2UL*gap);
            mpz_add_ui(increment_.get_mpz_t(),increment_.get_mpz_t(),gap*gap);
            difference_ += increment_;
            next_midpoint += gap;
            if (++index_ == residues_.size()) index_ = 0;
            covered_through = std::min(end,Z(next_midpoint-1));
        }
        completed = true;
        return finish("completed");
    }

private:
    void* native_ = nullptr;
    Z difference_, increment_, root_buffer_;
    std::vector<unsigned> residues_, gaps_;
    std::size_t index_ = 0;

    void read_native(const FusionSnapshot& snapshot) {
        static const char* names[] = {"ready","factor","exhausted","candidate_budget","time_budget","range_limit"};
        if (snapshot.status > 5) throw std::runtime_error("invalid native regional Fermat status");
        next_midpoint = z_from_u64(snapshot.next_midpoint);
        covered_through = std::min(end,z_from_u64(snapshot.covered_through));
        factor = z_from_u64(snapshot.factor);
        half_gap = z_from_u64(snapshot.half_gap);
        first_midpoint = z_from_u64(snapshot.first_midpoint);
        last_midpoint = z_from_u64(snapshot.last_midpoint);
        candidates = snapshot.candidates;
        elapsed_seconds = snapshot.elapsed_seconds;
        total_kernel_seconds = snapshot.total_kernel_seconds;
        completed = factor == 0 && next_midpoint > end;
        status = completed ? "completed" : names[snapshot.status];
    }
};

} // namespace vf
