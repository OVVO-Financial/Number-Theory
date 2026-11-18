#pragma once

#include "gmp_common.hpp"
#include "fermat_native_bridge.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace vf {

inline Z z_from_u64(std::uint64_t value) {
    Z result;
    mpz_import(result.get_mpz_t(), 1, 1, sizeof(value), 0, 0, &value);
    return result;
}

inline std::uint64_t z_to_u64(const Z& value) {
    if (value < 0 || mpz_sizeinbase(value.get_mpz_t(), 2) > 64)
        throw std::overflow_error("integer does not fit uint64_t");
    std::uint64_t result = 0;
    std::size_t words = 0;
    mpz_export(&result, &words, 1, sizeof(result), 0, 0, value.get_mpz_t());
    return result;
}

inline Z floor_sqrt(const Z& value) {
    if (value < 0) throw std::invalid_argument("square root requires a nonnegative integer");
    Z result;
    mpz_sqrt(result.get_mpz_t(), value.get_mpz_t());
    return result;
}

inline Z ceil_sqrt(const Z& value) {
    if (value < 0) throw std::invalid_argument("square root requires a nonnegative integer");
    Z result, remainder;
    mpz_sqrtrem(result.get_mpz_t(), remainder.get_mpz_t(), value.get_mpz_t());
    if (remainder != 0) ++result;
    return result;
}

inline Z mod(const Z& value, const Z& modulus) {
    if (modulus <= 0) throw std::invalid_argument("modulus must be positive");
    Z result;
    mpz_mod(result.get_mpz_t(), value.get_mpz_t(), modulus.get_mpz_t());
    return result;
}

inline Z floor_div(const Z& value, const Z& denominator) {
    if (denominator == 0) throw std::invalid_argument("division by zero");
    Z result;
    mpz_fdiv_q(result.get_mpz_t(), value.get_mpz_t(), denominator.get_mpz_t());
    return result;
}

inline Z factor_upper_bound(const Z& n, const Z& covered_through) {
    if (n < 2) throw std::invalid_argument("factor bound requires n >= 2");
    Z b = covered_through + 1;
    const Z ceiling = ceil_sqrt(n);
    if (b < ceiling) b = ceiling;
    return b - ceil_sqrt(Z(b*b - n));
}

inline Z fermat_midpoint_cap(const Z& n, const Z& lo) {
    if (n < 2 || lo <= 0)
        throw std::invalid_argument("midpoint cap requires n >= 2 and lo > 0");
    return floor_div(Z(lo*lo+n), Z(2*lo));
}

// Deterministic for 0 <= n < 2^64, using the same seven bases as the Python
// implementation. unsigned __int128 avoids overflow in modular products.
inline bool is_prime_u64(std::uint64_t n) {
    if (n < 2) return false;
    static constexpr std::array<std::uint64_t,12> small =
        {{2,3,5,7,11,13,17,19,23,29,31,37}};
    for (const auto p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    const auto multiply = [n](std::uint64_t a, std::uint64_t b) {
        return static_cast<std::uint64_t>(static_cast<unsigned __int128>(a)*b % n);
    };
    const auto power = [&multiply](std::uint64_t a, std::uint64_t e) {
        std::uint64_t result = 1;
        while (e) {
            if (e & 1) result = multiply(result,a);
            e >>= 1;
            if (e) a = multiply(a,a);
        }
        return result;
    };
    std::uint64_t d = n-1;
    unsigned s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }
    static constexpr std::array<std::uint64_t,7> bases =
        {{2,325,9375,28178,450775,9780504,1795265022}};
    for (const auto base : bases) {
        const std::uint64_t a = base % n;
        if (a == 0) continue;
        std::uint64_t x = power(a,d);
        if (x == 1 || x == n-1) continue;
        bool passed = false;
        for (unsigned r=1; r<s; ++r) {
            x = multiply(x,x);
            if (x == n-1) { passed = true; break; }
        }
        if (!passed) return false;
    }
    return true;
}

// GMP explicitly distinguishes definite primes (2), probable primes (1),
// and definite composites (0). Never promote code 1 to a prime certificate.
// https://gmplib.org/manual/Number-Theoretic-Functions
inline std::string prime_status(const Z& n, bool probable_check=true) {
    if (n < 2) return "composite";
    if (mpz_sizeinbase(n.get_mpz_t(),2) <= 64)
        return is_prime_u64(z_to_u64(n)) ? "prime" : "composite";
    if (!probable_check) return "unknown";
    const int status = mpz_probab_prime_p(n.get_mpz_t(),25);
    return status == 2 ? "prime" : status == 1 ? "probable_prime" : "composite";
}

// A persistent exact Fermat cursor for all positive integer sizes. The
// optional uint64 path reuses the existing native kernel without copying or
// changing its arithmetic. Larger values always remain in GMP integers.
// A failed square test and necessary-congruence gaps certify coverage;
// successful midpoints are never included in that excluded interval.
class GmpFermatCursor {
public:
    Z n, next_midpoint, covered_through, factor=0, half_gap=0;
    Z root_ceiling, endpoint, start_midpoint, first_midpoint=0, last_midpoint=0;
    std::uint64_t candidates=0;
    std::string status="ready", backend="gmp";
    bool exhausted=false;
    std::size_t residue_count=0;
    unsigned period=160;
    double elapsed_seconds=0.0, total_kernel_seconds=0.0;

    explicit GmpFermatCursor(const Z& number, bool native64=true) : n(number) {
        if (n < 3 || mpz_even_p(n.get_mpz_t()))
            throw std::invalid_argument("Fermat cursor requires an odd integer >= 3");
        const auto began = std::chrono::steady_clock::now();
        if (native64 && mpz_sizeinbase(n.get_mpz_t(),2) <= 64) {
            if (fusion_abi_version() != 1 || fusion_snapshot_size() != sizeof(FusionSnapshot))
                throw std::runtime_error("native Fermat ABI mismatch");
            unsigned error=0;
            native_ = fusion_fermat_create(z_to_u64(n), &error);
            if (!native_) throw std::runtime_error("native Fermat initialization failed");
            FusionSnapshot snapshot{};
            if (fusion_fermat_snapshot(native_,&snapshot)) {
                fusion_fermat_destroy(native_); native_=nullptr;
                throw std::runtime_error("native Fermat snapshot failed");
            }
            backend="native64";
            read_native(snapshot);
            return;
        }
        root_ceiling=ceil_sqrt(n);
        endpoint=(n+1)/2;
        std::array<bool,20> qr20{};
        std::array<bool,64> qr64{};
        for (unsigned i=0;i<20;++i) qr20[(i*i)%20]=true;
        for (unsigned i=0;i<64;++i) qr64[(i*i)%64]=true;
        const unsigned n20=mpz_fdiv_ui(n.get_mpz_t(),20);
        const unsigned n64=mpz_fdiv_ui(n.get_mpz_t(),64);
        for (unsigned r=0;r<160;++r)
            if (qr20[(r*r+20-n20)%20] && qr64[(r*r+64-n64)%64])
                residues_.push_back(r);
        if (residues_.empty()) throw std::logic_error("odd input has an empty Fermat wheel");
        residue_count=residues_.size();
        for (std::size_t i=0;i<residue_count;++i)
            gaps_.push_back(i+1<residue_count ? residues_[i+1]-residues_[i]
                                            :160+residues_[0]-residues_[i]);
        next_midpoint=root_ceiling;
        const unsigned remainder=mpz_fdiv_ui(next_midpoint.get_mpz_t(),160);
        const auto it=std::lower_bound(residues_.begin(),residues_.end(),remainder);
        if (it==residues_.end()) {
            index_=0;
            next_midpoint+=160-remainder+residues_[0];
        } else {
            index_=static_cast<std::size_t>(it-residues_.begin());
            next_midpoint+=*it-remainder;
        }
        start_midpoint=next_midpoint;
        covered_through=next_midpoint-1;
        if (covered_through>endpoint) covered_through=endpoint;
        difference_=next_midpoint*next_midpoint-n;
        total_kernel_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-began).count();
    }

    ~GmpFermatCursor() { if (native_) fusion_fermat_destroy(native_); }
    GmpFermatCursor(const GmpFermatCursor&)=delete;
    GmpFermatCursor& operator=(const GmpFermatCursor&)=delete;

    // Candidate and time limits are relative to this call. Negative seconds
    // disables the time limit; a temporary midpoint cap never exhausts state.
    // Coverage may extend across a congruence-only gap beyond that cap.
    std::uint64_t advance(std::uint64_t max_candidates, double timeout_seconds,
                          const Z& max_midpoint) {
        if (!std::isfinite(timeout_seconds))
            throw std::invalid_argument("Fermat timeout must be finite");
        if (max_midpoint < 0)
            throw std::invalid_argument("Fermat midpoint cap must be nonnegative");
        max_candidates=std::min(max_candidates,
                    std::numeric_limits<std::uint64_t>::max()-candidates);
        if (native_) {
            Z limit=max_midpoint;
            if (limit>endpoint) limit=endpoint;
            FusionSnapshot snapshot{};
            if (fusion_fermat_advance(native_,max_candidates,timeout_seconds,
                                      z_to_u64(limit),&snapshot))
                throw std::runtime_error("native Fermat advance failed");
            read_native(snapshot);
            return snapshot.delta_candidates;
        }

        const auto began=std::chrono::steady_clock::now();
        const auto before=candidates;
        const auto finish=[&](const char* reason) {
            status=reason;
            elapsed_seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-began).count();
            total_kernel_seconds+=elapsed_seconds;
            return candidates-before;
        };
        if (factor != 0) return finish("factor");
        if (exhausted) return finish("exhausted");
        if (max_candidates==0) return finish("candidate_budget");
        if (timeout_seconds==0.0) return finish("time_budget");
        const Z limit=max_midpoint<endpoint ? max_midpoint:endpoint;
        while (next_midpoint<=limit) {
            const auto used=candidates-before;
            if (used>=max_candidates) return finish("candidate_budget");
            if ((used&255U)==0 && timeout_seconds>=0.0 &&
                std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()
                    >=timeout_seconds)
                return finish("time_budget");
            ++candidates;
            if (candidates==1) first_midpoint=next_midpoint;
            last_midpoint=next_midpoint;
            if (mpz_perfect_square_p(difference_.get_mpz_t())) {
                mpz_sqrt(root_buffer_.get_mpz_t(),difference_.get_mpz_t());
                const Z divisor=next_midpoint-root_buffer_;
                if (divisor>1 && divisor<n) {
                    factor=divisor;
                    half_gap=root_buffer_;
                    return finish("factor");
                }
                if (divisor==1) {
                    covered_through=endpoint;
                    next_midpoint=endpoint+1;
                    exhausted=true;
                    return finish("exhausted");
                }
            }
            const unsigned gap=gaps_[index_];
            mpz_mul_ui(increment_.get_mpz_t(),next_midpoint.get_mpz_t(),2UL*gap);
            mpz_add_ui(increment_.get_mpz_t(),increment_.get_mpz_t(),gap*gap);
            difference_+=increment_;
            next_midpoint+=gap;
            if (++index_==residue_count) index_=0;
            covered_through=next_midpoint-1;
            if (covered_through>endpoint) covered_through=endpoint;
        }
        if (next_midpoint>endpoint) {
            exhausted=true;
            return finish("exhausted");
        }
        return finish("range_limit");
    }

private:
    void* native_=nullptr;
    Z difference_,increment_,root_buffer_;
    std::vector<unsigned> residues_,gaps_;
    std::size_t index_=0;

    void read_native(const FusionSnapshot& snapshot) {
        static const char* names[]={"ready","factor","exhausted","candidate_budget",
                                    "time_budget","range_limit"};
        if (snapshot.status>5) throw std::runtime_error("invalid native Fermat status");
        next_midpoint=z_from_u64(snapshot.next_midpoint);
        covered_through=z_from_u64(snapshot.covered_through);
        factor=z_from_u64(snapshot.factor);
        half_gap=z_from_u64(snapshot.half_gap);
        root_ceiling=z_from_u64(snapshot.root_ceiling);
        endpoint=z_from_u64(snapshot.endpoint);
        start_midpoint=z_from_u64(snapshot.start_midpoint);
        first_midpoint=z_from_u64(snapshot.first_midpoint);
        last_midpoint=z_from_u64(snapshot.last_midpoint);
        candidates=snapshot.candidates;
        status=names[snapshot.status];
        exhausted=snapshot.exhausted!=0;
        residue_count=snapshot.residue_count;
        period=snapshot.period;
        elapsed_seconds=snapshot.elapsed_seconds;
        total_kernel_seconds=snapshot.total_kernel_seconds;
    }
};

} // namespace vf
