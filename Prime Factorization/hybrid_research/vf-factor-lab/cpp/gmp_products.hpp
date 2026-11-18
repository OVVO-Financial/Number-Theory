#pragma once

// GMP-only polynomial products for exact searches over composite moduli.
//
// Large polynomial products use carry-free Kronecker substitution: nonnegative
// coefficients are packed into disjoint byte-aligned slots, GMP multiplies the
// resulting integers, and the slots are reduced modulo N. A coefficient before
// reduction is smaller than min(a.size(), b.size()) * 2^(2*bit_length(N)), so
// the chosen slot width prevents carries between polynomial coefficients.
//
// Multipoint evaluation uses a balanced monic product/remainder tree. Large
// monic remainders use reversal and Newton inversion of a series whose constant
// coefficient is 1. This requires no field assumption and is valid modulo a
// composite integer. Tiny operations use schoolbook arithmetic to avoid setup.

#include "gmp_common.hpp"
#include <climits>
#include <memory>
#include <set>
#include <utility>

namespace vf {
namespace gmp_poly {

using Poly = std::vector<Z>; // low coefficient first; empty means zero

inline void trim(Poly& a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

inline Z residue(const Z& a, const Z& n) {
    Z r;
    mpz_mod(r.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());
    return r;
}

inline void increment(Counters* counters, const char* key, std::uint64_t amount = 1) {
    if (counters) counters->inc(key, amount);
}

class ElapsedCounter {
    Counters* counters_;
    const char* name_;
    Clock::time_point start_ = Clock::now();
public:
    ElapsedCounter(Counters* counters, const char* name) : counters_(counters), name_(name) {}
    ~ElapsedCounter() { if (counters_) counters_->add_seconds(name_, seconds_since(start_)); }
};

class Ring {
    Z modulus_;
    Budget* budget_;
    Counters* counters_;
    std::size_t modulus_bits_;

    static std::size_t ceil_log2(std::size_t value) {
        std::size_t out = 0;
        for (--value; value; value >>= 1) ++out;
        return out;
    }

    Poly schoolbook(const Poly& a, const Poly& b) const {
        Poly out(a.size() + b.size() - 1);
        for (std::size_t i = 0; i < a.size(); ++i) {
            if ((i & 31) == 0) check_budget(budget_);
            for (std::size_t j = 0; j < b.size(); ++j)
                mpz_addmul(out[i + j].get_mpz_t(), a[i].get_mpz_t(), b[j].get_mpz_t());
        }
        for (auto& coefficient : out)
            mpz_mod(coefficient.get_mpz_t(), coefficient.get_mpz_t(), modulus_.get_mpz_t());
        trim(out);
        return out;
    }

    Poly kronecker(const Poly& a, const Poly& b) const {
        const std::size_t terms = std::min(a.size(), b.size());
        const std::size_t extra = ceil_log2(terms);
        if (modulus_bits_ > (std::numeric_limits<std::size_t>::max() - extra - 7) / 2)
            throw std::overflow_error("Kronecker slot size overflow");
        const std::size_t stride = (2 * modulus_bits_ + extra + 7) / 8;
        const std::size_t count = a.size() + b.size() - 1;
        if (count > std::numeric_limits<std::size_t>::max() / stride)
            throw std::overflow_error("Kronecker allocation size overflow");
        auto pack = [&](const Poly& coefficients) {
            std::vector<unsigned char> bytes(coefficients.size() * stride, 0);
            for (std::size_t i = 0; i < coefficients.size(); ++i) {
                if ((i & 127) == 0) check_budget(budget_);
                std::size_t written = 0;
                mpz_export(bytes.data() + i * stride, &written, -1, 1, 0, 0,
                           coefficients[i].get_mpz_t());
                // Ring operands are canonical residues, hence strictly fit.
                if (written > stride) throw std::logic_error("noncanonical polynomial coefficient");
            }
            Z encoded;
            mpz_import(encoded.get_mpz_t(), bytes.size(), -1, 1, 0, 0, bytes.data());
            return encoded;
        };
        Z left = pack(a), right = pack(b), encoded;
        check_budget(budget_);
        mpz_mul(encoded.get_mpz_t(), left.get_mpz_t(), right.get_mpz_t());
        check_budget(budget_);
        std::vector<unsigned char> bytes(count * stride, 0);
        std::size_t written = 0;
        mpz_export(bytes.data(), &written, -1, 1, 0, 0, encoded.get_mpz_t());
        if (written > bytes.size()) throw std::logic_error("Kronecker product overflow");
        Poly out(count);
        for (std::size_t i = 0; i < count; ++i) {
            if ((i & 127) == 0) check_budget(budget_);
            mpz_import(out[i].get_mpz_t(), stride, -1, 1, 0, 0, bytes.data() + i * stride);
            mpz_mod(out[i].get_mpz_t(), out[i].get_mpz_t(), modulus_.get_mpz_t());
        }
        trim(out);
        return out;
    }

public:
    explicit Ring(const Z& modulus, Budget* budget = nullptr, Counters* counters = nullptr)
        : modulus_(modulus), budget_(budget), counters_(counters),
          modulus_bits_(modulus > 0 ? mpz_sizeinbase(modulus.get_mpz_t(), 2) : 0) {
        if (modulus_ < 2) throw std::invalid_argument("polynomial modulus must be >= 2");
    }

    const Z& modulus() const { return modulus_; }

    Poly canonical(Poly a) const {
        for (auto& coefficient : a)
            mpz_mod(coefficient.get_mpz_t(), coefficient.get_mpz_t(), modulus_.get_mpz_t());
        trim(a);
        return a;
    }

    // Internal callers pass canonical nonnegative coefficient residues.
    Poly multiply(const Poly& a, const Poly& b) const {
        check_budget(budget_);
        increment(counters_, "polynomial_multiplications");
        if (a.empty() || b.empty()) return {};
        if (a.size() > std::numeric_limits<std::size_t>::max() - b.size() + 1)
            throw std::overflow_error("polynomial degree overflow");
        if (std::min(a.size(), b.size()) <= 4 ||
            (a.size() <= 32 && b.size() <= 32 && a.size() * b.size() <= 256)) {
            increment(counters_, "schoolbook_multiplications");
            return schoolbook(a, b);
        }
        increment(counters_, "kronecker_multiplications");
        return kronecker(a, b);
    }

    Poly multiply_truncated(const Poly& a, const Poly& b, std::size_t length) const {
        if (length == 0 || a.empty() || b.empty()) return {};
        Poly aa(a.begin(), a.begin() + std::min(a.size(), length));
        Poly bb(b.begin(), b.begin() + std::min(b.size(), length));
        Poly out = multiply(aa, bb);
        if (out.size() > length) out.resize(length);
        trim(out);
        return out;
    }

    Poly inverse_series_one(const Poly& f, std::size_t length) const {
        if (f.empty() || f[0] != 1)
            throw std::invalid_argument("series inversion requires constant coefficient 1");
        if (!length) return {};
        Poly inverse{Z(1)};
        std::size_t precision = 1;
        while (precision < length) {
            check_budget(budget_);
            const std::size_t next = precision > length - precision ? length : 2 * precision;
            Poly error = multiply_truncated(f, inverse, next);
            error.resize(next);
            error[0] = residue(Z(2) - error[0], modulus_);
            for (std::size_t j = 1; j < next; ++j)
                if (error[j] != 0) error[j] = modulus_ - error[j];
            trim(error);
            inverse = multiply_truncated(inverse, error, next);
            precision = next;
        }
        return inverse;
    }

    Poly remainder(const Poly& dividend, const Poly& divisor) const {
        check_budget(budget_);
        increment(counters_, "polynomial_remainders");
        if (divisor.empty() || divisor.back() != 1)
            throw std::invalid_argument("polynomial divisor must be monic");
        if (dividend.size() < divisor.size()) return dividend;
        const std::size_t degree = divisor.size() - 1;
        if (!degree) return {};
        const std::size_t quotient_length = dividend.size() - degree;
        if (degree <= 8 || quotient_length <= 4 ||
            (degree <= 128 && quotient_length <= 128 && degree * quotient_length <= 128)) {
            increment(counters_, "schoolbook_remainders");
            Poly out = dividend;
            for (std::size_t index = out.size(); index-- > degree;) {
                if ((index & 31) == 0) check_budget(budget_);
                const Z leading = out[index];
                if (leading != 0) {
                    const std::size_t offset = index - degree;
                    for (std::size_t j = 0; j < degree; ++j) {
                        mpz_submul(out[offset + j].get_mpz_t(), leading.get_mpz_t(), divisor[j].get_mpz_t());
                        mpz_mod(out[offset + j].get_mpz_t(), out[offset + j].get_mpz_t(), modulus_.get_mpz_t());
                    }
                }
            }
            out.resize(degree);
            trim(out);
            return out;
        }
        increment(counters_, "newton_remainders");
        Poly reversed_divisor;
        reversed_divisor.reserve(std::min(divisor.size(), quotient_length));
        for (std::size_t j = 0; j < std::min(divisor.size(), quotient_length); ++j)
            reversed_divisor.push_back(divisor[divisor.size() - 1 - j]);
        const Poly inverse = inverse_series_one(reversed_divisor, quotient_length);
        Poly reversed_dividend(quotient_length);
        for (std::size_t j = 0; j < quotient_length; ++j)
            reversed_dividend[j] = dividend[dividend.size() - 1 - j];
        Poly quotient = multiply_truncated(reversed_dividend, inverse, quotient_length);
        quotient.resize(quotient_length);
        std::reverse(quotient.begin(), quotient.end());
        trim(quotient);
        const Poly product = multiply_truncated(divisor, quotient, degree);
        Poly out(degree);
        for (std::size_t j = 0; j < degree; ++j)
            out[j] = residue(dividend[j] - (j < product.size() ? product[j] : Z(0)), modulus_);
        trim(out);
        return out;
    }

    Z evaluate(const Poly& polynomial, const Z& point) const {
        Z out = 0;
        for (std::size_t j = polynomial.size(); j-- > 0;) {
            if ((j & 127) == 0) check_budget(budget_);
            mpz_mul(out.get_mpz_t(), out.get_mpz_t(), point.get_mpz_t());
            mpz_add(out.get_mpz_t(), out.get_mpz_t(), polynomial[j].get_mpz_t());
            mpz_mod(out.get_mpz_t(), out.get_mpz_t(), modulus_.get_mpz_t());
        }
        return out;
    }

    Poly progression_polynomial(const Z& step, std::size_t degree) const {
        if (!degree) throw std::invalid_argument("block degree must be positive");
        if (degree > ULONG_MAX) throw std::overflow_error("block degree is too large");
        std::vector<Poly> level;
        level.reserve(degree);
        Z constant = 0;
        for (std::size_t i = 0; i < degree; ++i) {
            if ((i & 127) == 0) check_budget(budget_);
            level.push_back(Poly{constant, Z(1)});
            constant += step;
            mpz_mod(constant.get_mpz_t(), constant.get_mpz_t(), modulus_.get_mpz_t());
        }
        while (level.size() > 1) {
            check_budget(budget_);
            std::vector<Poly> following;
            following.reserve((level.size() + 1) / 2);
            for (std::size_t i = 0; i < level.size(); i += 2) {
                if (i + 1 < level.size()) following.push_back(multiply(level[i], level[i + 1]));
                else following.push_back(std::move(level[i]));
            }
            level = std::move(following);
        }
        return std::move(level[0]);
    }

    std::vector<Z> multipoint(const Poly& polynomial, const std::vector<Z>& points) const {
        check_budget(budget_);
        if (points.empty()) return {};
        if (points.size() <= 4) {
            std::vector<Z> out;
            out.reserve(points.size());
            for (const auto& point : points) out.push_back(evaluate(polynomial, point));
            return out;
        }
        std::vector<std::vector<Poly>> tree;
        tree.emplace_back();
        tree.back().reserve(points.size());
        for (const auto& point : points)
            tree.back().push_back(Poly{point == 0 ? Z(0) : modulus_ - point, Z(1)});
        while (tree.back().size() > 1) {
            check_budget(budget_);
            const auto& level = tree.back();
            std::vector<Poly> following;
            following.reserve((level.size() + 1) / 2);
            for (std::size_t i = 0; i < level.size(); i += 2) {
                if (i + 1 < level.size()) following.push_back(multiply(level[i], level[i + 1]));
                else following.push_back(level[i]);
            }
            tree.push_back(std::move(following));
        }
        std::vector<Poly> remainders{remainder(polynomial, tree.back()[0])};
        for (std::size_t level = tree.size() - 1; level-- > 0;) {
            check_budget(budget_);
            std::vector<Poly> following;
            following.reserve(tree[level].size());
            for (std::size_t i = 0; i < tree[level].size(); ++i)
                following.push_back(remainder(remainders[i / 2], tree[level][i]));
            remainders = std::move(following);
        }
        std::vector<Z> out;
        out.reserve(points.size());
        for (const auto& remainder_value : remainders)
            out.push_back(remainder_value.empty() ? Z(0) : remainder_value[0]);
        return out;
    }
};

} // namespace gmp_poly

struct Request {
    Z start;
    long step;
    std::size_t count;
    std::string kind;
};

class NormalizedProductEngine {
    Z n_;
    std::size_t degree_;
    Budget* budget_;
    Counters* counters_;
    bool shared_;
    gmp_poly::Ring ring_;
    std::map<Z, gmp_poly::Poly> polynomials_;
    std::map<Z, Z> inverses_{{Z(1), Z(1)}};
    std::set<Z> unit_steps_{Z(1)};

    Z canonical_step(long step) {
        const Z residue = gmp_poly::residue(Z(step), n_);
        if (!unit_steps_.count(residue)) {
            check_budget(budget_);
            Z factor;
            mpz_gcd(factor.get_mpz_t(), residue.get_mpz_t(), n_.get_mpz_t());
            gmp_poly::increment(counters_, "gcd_calls");
            gmp_poly::increment(counters_, "step_gcd_calls");
            if (factor != 1) throw std::invalid_argument("every progression step must be coprime to n");
            unit_steps_.insert(residue);
        }
        return residue;
    }

    const gmp_poly::Poly& polynomial(const Z& step) {
        const auto found = polynomials_.find(step);
        if (found != polynomials_.end()) return found->second;
        gmp_poly::ElapsedCounter timer(counters_, "polynomial_build_seconds");
        gmp_poly::increment(counters_, "polynomials_built");
        if (counters_ && degree_ > counters_->get("polynomial_degree"))
            counters_->inc("polynomial_degree", degree_ - counters_->get("polynomial_degree"));
        auto built = ring_.progression_polynomial(step, degree_);
        check_budget(budget_);
        return polynomials_.emplace(step, std::move(built)).first->second;
    }

    std::vector<Z> evaluate_blocks(const Z& step, const std::vector<Z>& points) {
        const auto& block = polynomial(step);
        gmp_poly::ElapsedCounter timer(counters_, "polynomial_evaluation_seconds");
        auto out = ring_.multipoint(block, points);
        gmp_poly::increment(counters_, "evaluations", points.size());
        gmp_poly::increment(counters_, "represented_candidates", points.size() * degree_);
        check_budget(budget_);
        return out;
    }

public:
    NormalizedProductEngine(const Z& n, std::size_t degree, Budget* budget = nullptr,
                            Counters* counters = nullptr, bool shared = true)
        : n_(n), degree_(degree), budget_(budget), counters_(counters), shared_(shared),
          ring_(n, budget, counters) {
        if (n_ < 2) throw std::invalid_argument("modulus must be >= 2");
        if (!degree_) throw std::invalid_argument("block degree must be positive");
    }

    std::vector<Z> products(const std::vector<Z>& starts, const std::vector<long>& steps) {
        if (starts.size() != steps.size()) throw std::invalid_argument("starts and steps must have equal lengths");
        if (starts.empty()) return {};
        if (starts.size() > std::numeric_limits<std::uint64_t>::max() / degree_)
            throw std::overflow_error("represented-work count overflow");
        check_budget(budget_);
        std::vector<Z> canonical_steps;
        canonical_steps.reserve(steps.size());
        std::set<Z> distinct_steps;
        for (long step : steps) {
            canonical_steps.push_back(canonical_step(step));
            distinct_steps.insert(canonical_steps.back());
        }
        gmp_poly::increment(counters_, "product_batches");
        if (distinct_steps.size() > 1) gmp_poly::increment(counters_, "mixed_batches");
        if (shared_) {
            std::vector<Z> points;
            points.reserve(starts.size());
            for (std::size_t i = 0; i < starts.size(); ++i) {
                if ((i & 127) == 0) check_budget(budget_);
                const auto& step = canonical_steps[i];
                auto found = inverses_.find(step);
                if (found == inverses_.end()) {
                    Z inverse;
                    if (!mpz_invert(inverse.get_mpz_t(), step.get_mpz_t(), n_.get_mpz_t()))
                        throw std::logic_error("validated step unexpectedly noninvertible");
                    found = inverses_.emplace(step, std::move(inverse)).first;
                    gmp_poly::increment(counters_, "step_inversions");
                }
                points.push_back(gmp_poly::residue(starts[i] * found->second, n_));
            }
            gmp_poly::increment(counters_, "normalized_blocks", points.size());
            return evaluate_blocks(Z(1), points);
        }
        std::map<Z, std::vector<std::size_t>> groups;
        for (std::size_t i = 0; i < starts.size(); ++i) groups[canonical_steps[i]].push_back(i);
        std::vector<Z> out(starts.size());
        for (const auto& [step, indices] : groups) {
            check_budget(budget_);
            std::vector<Z> points;
            points.reserve(indices.size());
            for (auto i : indices) points.push_back(gmp_poly::residue(starts[i], n_));
            auto values = evaluate_blocks(step, points);
            for (std::size_t j = 0; j < indices.size(); ++j) out[indices[j]] = std::move(values[j]);
        }
        check_budget(budget_);
        return out;
    }
};

inline Z recover_mixed_batch(const Z& n, const std::vector<Request>& requests,
                            const std::vector<Z>& values, Budget* budget = nullptr,
                            Counters* counters = nullptr) {
    if (n < 2) throw std::invalid_argument("modulus must be >= 2");
    if (requests.size() != values.size()) throw std::invalid_argument("requests and values must have equal lengths");
    if (requests.empty()) return Z(0);
    check_budget(budget);
    std::vector<Z> reduced;
    reduced.reserve(values.size());
    Z combined = 1;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if ((i & 127) == 0) check_budget(budget);
        reduced.push_back(gmp_poly::residue(values[i], n));
        combined *= reduced.back();
        mpz_mod(combined.get_mpz_t(), combined.get_mpz_t(), n.get_mpz_t());
    }
    check_budget(budget);
    Z factor;
    mpz_gcd(factor.get_mpz_t(), combined.get_mpz_t(), n.get_mpz_t());
    gmp_poly::increment(counters, "gcd_calls");
    gmp_poly::increment(counters, "batch_gcd_calls");
    if (factor > 1 && factor < n) return factor;
    if (factor == 1) { check_budget(budget); return Z(0); }
    gmp_poly::increment(counters, "gcd_n_batches");
    for (std::size_t i = 0; i < requests.size(); ++i) {
        check_budget(budget);
        mpz_gcd(factor.get_mpz_t(), reduced[i].get_mpz_t(), n.get_mpz_t());
        gmp_poly::increment(counters, "gcd_calls");
        gmp_poly::increment(counters, "recovery_block_gcd_calls");
        gmp_poly::increment(counters, "blocks_checked");
        if (factor > 1 && factor < n) return factor;
        if (factor == n) {
            gmp_poly::increment(counters, "gcd_n_blocks");
            Z candidate = requests[i].start;
            for (std::size_t j = 0; j < requests[i].count; ++j) {
                if ((j & 127) == 0) check_budget(budget);
                mpz_gcd(factor.get_mpz_t(), candidate.get_mpz_t(), n.get_mpz_t());
                gmp_poly::increment(counters, "gcd_calls");
                gmp_poly::increment(counters, "rescue_candidates");
                if (factor > 1 && factor < n) return factor;
                candidate += requests[i].step;
            }
        }
    }
    check_budget(budget);
    return Z(0);
}

} // namespace vf
