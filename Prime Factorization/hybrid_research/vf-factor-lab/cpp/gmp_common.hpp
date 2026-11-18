#pragma once

#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace vf {

using Z = mpz_class;
using Clock = std::chrono::steady_clock;

inline double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

class BudgetExceeded : public std::runtime_error {
public:
    explicit BudgetExceeded(const std::string& reason) : std::runtime_error(reason) {}
};

// One instance may be shared by every recursive cofactor split. Work units
// count represented progression positions, scalar divisors, and Fermat tests.
class Budget {
    Clock::time_point started_;
    Clock::time_point deadline_;
    bool timed_;
    std::optional<std::uint64_t> limit_;
    std::uint64_t used_ = 0;
public:
    explicit Budget(double timeout_seconds = -1,
                    std::optional<std::uint64_t> max_work = std::nullopt)
        : started_(Clock::now()), deadline_(started_), timed_(timeout_seconds >= 0),
          limit_(max_work) {
        if (!std::isfinite(timeout_seconds) || timeout_seconds < -1)
            throw std::invalid_argument("timeout must be finite and nonnegative, or -1 for unlimited");
        if (timed_) {
            const double capacity = std::chrono::duration<double>(Clock::time_point::max() - started_).count();
            if (timeout_seconds >= capacity) deadline_ = Clock::time_point::max();
            else deadline_ += std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(timeout_seconds));
        }
    }
    void check() const {
        if (timed_ && Clock::now() >= deadline_) throw BudgetExceeded("time budget");
    }
    void reserve(std::uint64_t amount) const {
        check();
        if (amount > std::numeric_limits<std::uint64_t>::max() - used_ ||
            (limit_ && amount > *limit_ - used_))
            throw BudgetExceeded("represented-work budget");
    }
    // Charge completed work even if its arithmetic crossed the deadline.
    void charge(std::uint64_t amount) {
        if (amount > std::numeric_limits<std::uint64_t>::max() - used_ ||
            (limit_ && amount > *limit_ - used_))
            throw std::logic_error("work charged without a sufficient reservation");
        used_ += amount;
    }
    double remaining() const {
        return timed_ ? std::max(0.0, std::chrono::duration<double>(deadline_ - Clock::now()).count()) : -1.0;
    }
    double elapsed() const { return seconds_since(started_); }
    std::uint64_t used() const { return used_; }
    std::uint64_t available_work() const {
        return limit_ ? *limit_ - used_ : std::numeric_limits<std::uint64_t>::max() - used_;
    }
};

inline void check_budget(Budget* budget) { if (budget) budget->check(); }

struct Counters {
    std::map<std::string, std::uint64_t> counts;
    std::map<std::string, double> seconds;
    std::map<std::string, Z> integers;
    void inc(const std::string& key, std::uint64_t amount = 1) {
        auto& value = counts[key];
        if (amount > std::numeric_limits<std::uint64_t>::max() - value)
            throw std::overflow_error("counter overflow: " + key);
        value += amount;
    }
    void add_seconds(const std::string& key, double amount) { seconds[key] += amount; }
    std::uint64_t get(const std::string& key) const {
        const auto found = counts.find(key);
        return found == counts.end() ? 0 : found->second;
    }
};

} // namespace vf
