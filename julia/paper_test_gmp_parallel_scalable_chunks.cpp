// C++17 GMP parallel seed-rho-chunk scalable spectrum staged multi-seed batch-GCD accelerator
//
// This version keeps the vertical-line search bottom-started, but adds
// adaptive per-k window scheduling:
//
//   low k  -> deep window
//   high k -> shallow window
//
// The default adaptive target is a geometric bottom-depth cap:
//
//   b_cap = floor(sqrt(b_max)), where b_max = (N - 9) / 6
//
// For each k <= low_k_max:
//
//   window_k = ceil(b_cap / (5 * 2^k))
//
// optionally capped by window_step_cap and never below base_window_steps.
//
// Compile in MSYS2 MINGW64:
// g++ -std=c++17 -O3 paper_test_gmp_parallel_seedrho.cpp -lgmpxx -lgmp -pthread -o paper_test_gmp_parallel_seedrho.exe
//
// Run:
// ./paper_test_gmp_parallel_seedrho.exe
//
// Arguments:
// ./paper_test_gmp_parallel_seedrho.exe <base_window_steps> <seed_count> <batch_size> <max_k> <thread_count> <seed_start> <digit_mode> <centered> <schedule_mode> <low_k_max> <b_cap_mode> <target_divisor> <manual_b_cap> <window_step_cap> <rho_center_den> <rho_radius_den> <rho_spectrum_mode> <chunk_steps> <chunk_start> <chunk_count>
//
// digit_mode:
//   0 = Fermat terminal b digits
//   1 = parity-compatible b digits
//   2 = all b digits
//
// centered:
//   0 = one-sided seeds: a_s = a0 + 2(seed_start + s)
//   1 = centered seeds:  a_s = a0 + 2 * offset, offset = 0,+1,-1,+2,-2,...
//
// schedule_mode:
//   0 = fixed window for every stream
//   1 = adaptive low-k deepening, bottom-started
//   2 = single rho-shell outward search
//   3 = rho-spectrum seed-rho tiled search [default]
//
// max_k = 0 means automatic.
// thread_count = 0 means automatic.
// b_cap_mode:
//   0 = sqrt(N) / target_divisor
//   1 = sqrt(b_max)   [default]
//   2 = manual_b_cap
//
// rho_spectrum_mode:
//   0 = beta_N scaled balanced spectrum [default]
//       rho_beta = (bits(N) / maxbits(decimal_digits(N))) * (1 - 1/sqrt(2))
//       shells: rho_beta, rho_beta/2, rho_beta/4, rho_beta/8, ...
//   1 = fixed balanced hard-shell spectrum: 0.293, 0.25, 0.20, 0.125, 0.0625, 0.03125
//   2 = general spectrum: 0.50, 0.375, 0.293, 0.25, 0.20, 0.125, 0.0625
//
// window_step_cap = 0 means no theoretical shell cap. Use with care for large N.
// chunk_steps is the executable per-stream chunk size.
// chunk_start is the first lower-first outward chunk-block index.
// chunk_count = 0 means automatically use worker thread count.
// Outward attempt order is: center, lower, upper, lower, upper, ...

#include <gmpxx.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using BigInt = mpz_class;

struct AcceleratorStream {
    std::uint64_t seed_index = 0;
    BigInt seed_a;

    int k = 0;
    int b_digit = 0;

    BigInt b_start;
    BigInt b_residue;

    BigInt b_current;
    BigInt stride_b;

    BigInt M_current;
    BigInt stride_M;

    BigInt A_num;

    std::uint64_t window_limit = 0;
    std::uint64_t full_window_limit = 0;
    std::uint64_t chunk_index = 0;
    std::uint64_t chunk_attempt_start = 0;
    bool use_outward = false;
    int shell_index = -1;
    BigInt shell_center_b;
    BigInt shell_radius_b;

    bool is_subtraction = true;
};

struct Candidate {
    std::uint64_t t = 0;

    std::uint64_t seed_index = 0;
    BigInt seed_a;

    int k = 0;
    int b_digit = 0;

    BigInt b_start;
    BigInt b_residue;

    BigInt b_current;
    BigInt stride_b;

    BigInt M_current;

    bool is_subtraction = true;
    std::uint64_t chunk_index = 0;
    std::uint64_t chunk_attempt_start = 0;
    int shell_index = -1;
    BigInt shell_center_b;
    BigInt shell_radius_b;
};

struct FactorHit {
    BigInt p;
    BigInt q;
    BigInt a;
    BigInt b;
    BigInt Y;

    int level_k = 0;
    BigInt gcd_coordinate;
    std::uint64_t hit_t = 0;

    bool has_stream_metadata = false;
    bool hit_is_subtraction = false;

    std::uint64_t seed_index = 0;
    BigInt seed_a;

    int hit_b_digit = -1;

    BigInt hit_b_residue;
    BigInt hit_b_start;
    BigInt hit_b_current;
    BigInt hit_stride_b;

    BigInt true_b_mod_stride;

    BigInt collision_value;
    BigInt abs_collision_value;
    BigInt proxy_gcd;

    BigInt gcd_result;
    BigInt gcd_multiple;

    int hit_shell_index = -1;
    std::uint64_t hit_chunk_index = 0;
    std::uint64_t hit_chunk_attempt_start = 0;
    BigInt hit_shell_center_b;
    BigInt hit_shell_radius_b;
};

struct ThreadStats {
    std::uint64_t proxy_candidates_tested = 0;
    std::uint64_t batch_gcd_checks = 0;
    std::uint64_t recovery_gcd_checks = 0;
    std::uint64_t active_streams = 0;
};

struct SearchStats {
    std::uint64_t proxy_candidates_tested = 0;
    std::uint64_t batch_gcd_checks = 0;
    std::uint64_t recovery_gcd_checks = 0;
    std::uint64_t active_streams = 0;

    std::uint64_t seed_count = 0;
    std::uint64_t seed_start = 0;

    int max_k = 0;
    int digit_mode = 0;
    bool centered = false;

    int schedule_mode = 2;
    std::uint64_t low_k_max = 8;
    int b_cap_mode = 1;
    std::uint64_t target_divisor = 16;
    BigInt b_cap;
    BigInt manual_b_cap;
    std::uint64_t rho_center_den = 16;
    std::uint64_t rho_radius_den = 64;
    int rho_spectrum_mode = 0;
    std::uint64_t decimal_digits_N = 0;
    std::uint64_t bit_length_N = 0;
    std::uint64_t maxbits_digits_N = 0;
    std::uint64_t beta_rho_num = 0;
    std::uint64_t beta_rho_den = 1;
    std::uint64_t rho_shell_count = 0;
    std::uint64_t chunk_steps = 1000000;
    std::uint64_t chunk_start = 0;
    std::uint64_t chunk_count = 0;
    std::uint64_t total_lanes = 0;
    BigInt rho_center_b;
    BigInt rho_radius_b;
    std::uint64_t window_step_cap = 0;
    std::uint64_t max_stream_window = 0;

    unsigned int thread_count = 1;

    std::uint64_t runtime_ms = 0;
};

static BigInt big_from_u64(std::uint64_t x) {
    BigInt v;
    v.set_str(std::to_string(x), 10);
    return v;
}

static std::uint64_t maxbits_for_decimal_digits(std::uint64_t digits) {
    // maxbits(D) = ceil(D * log2(10)).
    // Use a rational upper-quality approximation to avoid floating point.
    // log2(10) ~= 3.32192809489 = 332192809489 / 100000000000.
    const unsigned long long NUM = 332192809489ULL;
    const unsigned long long DEN = 100000000000ULL;

#if defined(__SIZEOF_INT128__)
    __uint128_t prod = static_cast<__uint128_t>(digits) * static_cast<__uint128_t>(NUM);
    __uint128_t val = (prod + DEN - 1) / DEN;
    if (val > std::numeric_limits<std::uint64_t>::max()) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    return static_cast<std::uint64_t>(val);
#else
    long double v = static_cast<long double>(digits) * 3.32192809488736234787L;
    if (v >= static_cast<long double>(std::numeric_limits<std::uint64_t>::max())) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    return static_cast<std::uint64_t>(v) + ((v > static_cast<std::uint64_t>(v)) ? 1 : 0);
#endif
}

static BigInt signed_seed_offset(std::uint64_t global_seed_index, bool centered) {
    if (!centered) {
        return big_from_u64(global_seed_index);
    }

    if (global_seed_index == 0) {
        return BigInt(0);
    }

    if (global_seed_index % 2 == 1) {
        return big_from_u64((global_seed_index + 1) / 2);
    }

    return -big_from_u64(global_seed_index / 2);
}

static BigInt abs_big(const BigInt& x) {
    return x < 0 ? -x : x;
}

static BigInt mod_pos(BigInt x, const BigInt& m) {
    x %= m;
    if (x < 0) {
        x += m;
    }
    return x;
}

static BigInt pow2_big(int k) {
    BigInt x = 0;
    mpz_setbit(x.get_mpz_t(), static_cast<unsigned long>(k));
    return x;
}

static BigInt div_pow2_floor(const BigInt& x, int k) {
    BigInt q;
    mpz_fdiv_q_2exp(
        q.get_mpz_t(),
        x.get_mpz_t(),
        static_cast<unsigned long>(k)
    );
    return q;
}

static BigInt gcd_int(const BigInt& a, const BigInt& b) {
    BigInt g;
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return g;
}

static BigInt isqrt_floor(const BigInt& n) {
    if (n < 0) {
        throw std::runtime_error("negative square root");
    }

    BigInt r;
    mpz_sqrt(r.get_mpz_t(), n.get_mpz_t());
    return r;
}

static BigInt isqrt_ceil(const BigInt& n) {
    BigInt r = isqrt_floor(n);
    return (r * r == n) ? r : r + 1;
}

static bool is_perfect_square(const BigInt& n) {
    if (n < 0) {
        return false;
    }

    return mpz_perfect_square_p(n.get_mpz_t()) != 0;
}

static std::optional<BigInt> mod_inverse(const BigInt& a, const BigInt& m) {
    BigInt inv;
    BigInt base = mod_pos(a, m);

    int ok = mpz_invert(
        inv.get_mpz_t(),
        base.get_mpz_t(),
        m.get_mpz_t()
    );

    if (!ok) {
        return std::nullopt;
    }

    return mod_pos(inv, m);
}

static BigInt first_geq_with_residue(
    const BigInt& low,
    const BigInt& residue,
    const BigInt& step
) {
    BigInt delta = mod_pos(residue - mod_pos(low, step), step);
    return low + delta;
}

static BigInt nearest_with_residue_to_anchor(
    const BigInt& anchor,
    const BigInt& residue,
    const BigInt& step,
    const BigInt& low,
    const BigInt& high
) {
    BigInt lower_delta = mod_pos(anchor - residue, step);
    BigInt lower = anchor - lower_delta;
    BigInt upper = lower + step;

    bool lower_ok = (lower >= low && lower <= high);
    bool upper_ok = (upper >= low && upper <= high);

    if (lower_ok && upper_ok) {
        BigInt dl = abs_big(anchor - lower);
        BigInt du = abs_big(upper - anchor);
        return (dl <= du) ? lower : upper;
    }

    if (lower_ok) {
        return lower;
    }

    if (upper_ok) {
        return upper;
    }

    return first_geq_with_residue(low, residue, step);
}

static BigInt outward_b_from_anchor(
    const BigInt& anchor,
    const BigInt& step,
    std::uint64_t index
) {
    // Correct rho-shell outward mapping, biased inward / decreasing first:
    //
    //   0 -> center
    //   1 -> center - 1 step
    //   2 -> center + 1 step
    //   3 -> center - 2 steps
    //   4 -> center + 2 steps
    //
    // This is the desired mapping for beta-balanced rho shells, because the
    // highest-priority rho shells often begin at an upper hard-boundary and
    // the productive proxy lies inward from that shell center.
    if (index == 0) {
        return anchor;
    }

    if (index % 2 == 1) {
        return anchor - big_from_u64((index + 1) / 2) * step;
    }

    return anchor + big_from_u64(index / 2) * step;
}


struct RhoSpec {
    std::uint64_t numerator;
    std::uint64_t denominator;
};

struct RhoShell {
    int shell_index = -1;
    BigInt center_b;
    BigInt radius_b;
};

static BigInt rho_to_b(
    const BigInt& sqrtN,
    std::uint64_t numerator,
    std::uint64_t denominator
) {
    if (denominator == 0) {
        denominator = 1;
    }

    BigInt b = (sqrtN * big_from_u64(numerator)) / big_from_u64(denominator);
    return b > 0 ? b : BigInt(1);
}

static std::vector<RhoSpec> rho_spectrum_specs(
    int spectrum_mode,
    std::uint64_t bit_length_N,
    std::uint64_t maxbits_digits_N
) {
    if (spectrum_mode == 0) {
        // Beta_N scaled balanced spectrum.
        // beta_N = bits(N) / maxbits(decimal_digits(N))
        // rho_beta = beta_N * (1 - 1/sqrt(2))
        // Use 1 - 1/sqrt(2) ~= 292893 / 1000000.
        if (maxbits_digits_N == 0) {
            maxbits_digits_N = bit_length_N == 0 ? 1 : bit_length_N;
        }

        const std::uint64_t SAME_BIT_RHO_NUM = 292893ULL;
        const std::uint64_t SAME_BIT_RHO_DEN = 1000000ULL;

        std::uint64_t base_num = bit_length_N * SAME_BIT_RHO_NUM;
        std::uint64_t base_den = maxbits_digits_N * SAME_BIT_RHO_DEN;

        if (base_num == 0) {
            base_num = 1;
        }
        if (base_den == 0) {
            base_den = 1;
        }

        return {
            {base_num, base_den},
            {base_num, base_den * 2ULL},
            {base_num, base_den * 4ULL},
            {base_num, base_den * 8ULL},
            {base_num, base_den * 16ULL},
            {base_num, base_den * 32ULL}
        };
    }

    if (spectrum_mode == 1) {
        // Fixed balanced hard-shell spectrum.
        // Starts near u = 1 - 1/sqrt(2) ~= 0.2929, the hardest same-bit boundary,
        // then emanates inward toward Fermat-easy territory.
        return {
            {293, 1000},
            {1, 4},
            {1, 5},
            {1, 8},
            {1, 16},
            {1, 32}
        };
    }

    if (spectrum_mode == 2) {
        // General semiprime spectrum.
        // Includes the broader hard shell around u ~= 0.5, then steps inward.
        return {
            {1, 2},
            {3, 8},
            {293, 1000},
            {1, 4},
            {1, 5},
            {1, 8},
            {1, 16}
        };
    }

    throw std::runtime_error("rho_spectrum_mode must be 0, 1, or 2");
}

static std::vector<RhoShell> build_rho_shells(
    const BigInt& sqrtN,
    const BigInt& b_max,
    int schedule_mode,
    std::uint64_t rho_center_den,
    std::uint64_t rho_radius_den,
    int rho_spectrum_mode,
    std::uint64_t bit_length_N,
    std::uint64_t maxbits_digits_N
) {
    std::vector<RhoShell> shells;

    if (schedule_mode == 2) {
        RhoShell shell;
        shell.shell_index = 0;
        shell.center_b = rho_to_b(sqrtN, 1, rho_center_den == 0 ? 16 : rho_center_den);
        shell.radius_b = rho_to_b(sqrtN, 1, rho_radius_den == 0 ? 128 : rho_radius_den);

        if (shell.center_b > b_max) {
            shell.center_b = b_max;
        }

        if (shell.radius_b <= 0) {
            shell.radius_b = 1;
        }

        shells.push_back(shell);
        return shells;
    }

    if (schedule_mode == 3) {
        auto specs = rho_spectrum_specs(rho_spectrum_mode, bit_length_N, maxbits_digits_N);
        std::uint64_t radius_den = rho_radius_den == 0 ? 64 : rho_radius_den;

        int idx = 0;
        for (const auto& spec : specs) {
            RhoShell shell;
            shell.shell_index = idx++;
            shell.center_b = rho_to_b(sqrtN, spec.numerator, spec.denominator);
            shell.radius_b = rho_to_b(sqrtN, 1, radius_den);

            if (shell.center_b > b_max) {
                shell.center_b = b_max;
            }

            if (shell.radius_b <= 0) {
                shell.radius_b = 1;
            }

            shells.push_back(shell);
        }

        return shells;
    }

    RhoShell shell;
    shell.shell_index = -1;
    shell.center_b = BigInt(1);
    shell.radius_b = BigInt(1);
    shells.push_back(shell);
    return shells;
}

static std::uint64_t big_to_u64_cap(const BigInt& x, std::uint64_t cap) {
    if (x <= 0) {
        return 0;
    }

    if (cap > 0 && x > big_from_u64(cap)) {
        return cap;
    }

    std::string s = x.get_str();

    if (s.size() > 20) {
        return std::numeric_limits<std::uint64_t>::max();
    }

    try {
        return static_cast<std::uint64_t>(std::stoull(s));
    } catch (...) {
        return std::numeric_limits<std::uint64_t>::max();
    }
}

static std::uint64_t ceil_div_big_to_u64_cap(
    const BigInt& numerator,
    const BigInt& denominator,
    std::uint64_t cap
) {
    if (numerator <= 0) {
        return 0;
    }

    BigInt q = (numerator + denominator - 1) / denominator;
    return big_to_u64_cap(q, cap);
}

static BigInt compute_b_cap(
    const BigInt& sqrtN,
    const BigInt& b_max,
    int b_cap_mode,
    std::uint64_t target_divisor,
    const BigInt& manual_b_cap
) {
    if (b_cap_mode == 0) {
        if (target_divisor == 0) {
            target_divisor = 16;
        }

        BigInt cap = sqrtN / big_from_u64(target_divisor);
        return cap > 0 ? cap : BigInt(1);
    }

    if (b_cap_mode == 1) {
        BigInt cap = isqrt_floor(b_max);
        return cap > 0 ? cap : BigInt(1);
    }

    if (b_cap_mode == 2) {
        if (manual_b_cap <= 0) {
            return BigInt(1);
        }

        return manual_b_cap;
    }

    throw std::runtime_error("b_cap_mode must be 0, 1, or 2");
}


static std::uint64_t outward_attempt_window_from_radius_steps(
    std::uint64_t one_side_steps,
    std::uint64_t window_step_cap
) {
    // An outward search visits:
    //
    //   center, lower1, upper1, lower2, upper2, ...
    //
    // To cover R stride-steps on the lower side AND R stride-steps on the
    // upper side, the number of attempts must be approximately 2R + 1.
    //
    // The earlier chunked code used R attempts directly. That only covered
    // about R/2 steps on each side, so hits inside the declared rho radius
    // could be skipped.
    if (one_side_steps >= (std::numeric_limits<std::uint64_t>::max() - 1) / 2) {
        return std::numeric_limits<std::uint64_t>::max();
    }

    std::uint64_t attempts = one_side_steps * 2 + 1;

    if (window_step_cap > 0 && attempts > window_step_cap) {
        return window_step_cap;
    }

    return attempts == 0 ? 1 : attempts;
}

static std::uint64_t adaptive_window_for_k(
    const BigInt& b_cap,
    int k,
    std::uint64_t base_window_steps,
    int schedule_mode,
    std::uint64_t low_k_max,
    std::uint64_t window_step_cap
) {
    std::uint64_t result = base_window_steps;

    if (schedule_mode == 0) {
        return result == 0 ? 1 : result;
    }

    if (k <= 0) {
        return result == 0 ? 1 : result;
    }

    if (static_cast<std::uint64_t>(k) > low_k_max) {
        return result == 0 ? 1 : result;
    }

    BigInt stride_b = BigInt(5) * pow2_big(k);

    std::uint64_t needed = ceil_div_big_to_u64_cap(
        b_cap,
        stride_b,
        window_step_cap == 0 ? std::numeric_limits<std::uint64_t>::max() : window_step_cap
    );

    if (schedule_mode == 2) {
        if (needed < (std::numeric_limits<std::uint64_t>::max() / 2)) {
            needed = needed * 2 + 1;
        }
    }

    if (needed > result) {
        result = needed;
    }

    if (window_step_cap > 0 && result > window_step_cap) {
        result = window_step_cap;
    }

    if (result == 0) {
        result = 1;
    }

    return result;
}

// Solve:
//     x = db mod 10
//     x = r  mod 2^k
//
// Returns x mod lcm(10, 2^k) = 5 * 2^k when compatible.
static std::optional<BigInt> crt_b_residue(
    int db,
    const BigInt& r,
    int k
) {
    if (k <= 0) {
        return BigInt(db);
    }

    BigInt pow2 = pow2_big(k);
    BigInt diff = r - BigInt(db);

    // gcd(10, 2^k) = 2 for k >= 1.
    if (mod_pos(diff, BigInt(2)) != 0) {
        return std::nullopt;
    }

    BigInt reduced_mod = pow2 / 2;
    BigInt t = 0;

    if (reduced_mod != 1) {
        BigInt rhs = mod_pos(diff / 2, reduced_mod);

        auto inv5 = mod_inverse(BigInt(5), reduced_mod);
        if (!inv5) {
            return std::nullopt;
        }

        t = mod_pos(rhs * (*inv5), reduced_mod);
    }

    BigInt lcm = 5 * pow2;
    BigInt x = BigInt(db) + 10 * t;

    return mod_pos(x, lcm);
}

static std::optional<bool> required_a_even_for_N(const BigInt& N) {
    int n4 = static_cast<int>(mpz_tdiv_ui(N.get_mpz_t(), 4));

    if (n4 == 3) {
        return true;
    }

    if (n4 == 1) {
        return false;
    }

    return std::nullopt;
}

static bool required_b_even_for_N(const BigInt& N) {
    auto req_a_even = required_a_even_for_N(N);
    if (!req_a_even) {
        throw std::runtime_error("N must be odd");
    }

    // N = 4k - 1: a even, b odd.
    // N = 4k + 1: a odd,  b even.
    return !(*req_a_even);
}

static BigInt parity_aligned_seed(const BigInt& N) {
    BigInt a0 = isqrt_ceil(N);

    auto req = required_a_even_for_N(N);
    if (!req) {
        throw std::runtime_error("N must be odd");
    }

    bool required_even = *req;
    bool a0_even = (mpz_tdiv_ui(a0.get_mpz_t(), 2) == 0);

    if (a0_even != required_even) {
        ++a0;
    }

    return a0;
}

// Derive Fermat terminal digit pairs directly:
//
//     p = a - b
//     q = a + b
//     p * q = N mod 10
//
// Also applies:
//
//     N = 4k - 1 gives a even, b odd
//     N = 4k + 1 gives a odd,  b even.
static std::vector<std::pair<int, int>> fermat_digit_pairs(const BigInt& N) {
    std::vector<std::pair<int, int>> pairs;

    int n10 = static_cast<int>(mpz_tdiv_ui(N.get_mpz_t(), 10));
    int n4  = static_cast<int>(mpz_tdiv_ui(N.get_mpz_t(), 4));

    bool a_even;
    bool b_even;

    if (n4 == 3) {
        a_even = true;
        b_even = false;
    } else if (n4 == 1) {
        a_even = false;
        b_even = true;
    } else {
        throw std::runtime_error("N must be odd for this Fermat lattice routine");
    }

    auto valid_odd_factor_digit = [](int d) {
        d = ((d % 10) + 10) % 10;
        return d == 1 || d == 3 || d == 7 || d == 9;
    };

    for (int da = 0; da <= 9; ++da) {
        if (((da % 2) == 0) != a_even) {
            continue;
        }

        for (int db = 0; db <= 9; ++db) {
            if (((db % 2) == 0) != b_even) {
                continue;
            }

            int p10 = ((da - db) % 10 + 10) % 10;
            int q10 = (da + db) % 10;

            if (!valid_odd_factor_digit(p10)) {
                continue;
            }

            if (!valid_odd_factor_digit(q10)) {
                continue;
            }

            if ((p10 * q10) % 10 == n10) {
                pairs.push_back({da, db});
            }
        }
    }

    return pairs;
}

static std::vector<int> b_digits_for_mode(const BigInt& N, int digit_mode) {
    std::set<int> digit_set;

    if (digit_mode == 0) {
        auto pairs = fermat_digit_pairs(N);

        for (auto [da, db] : pairs) {
            digit_set.insert(db);
        }
    } else if (digit_mode == 1) {
        bool b_even = required_b_even_for_N(N);

        for (int d = 0; d <= 9; ++d) {
            if (((d % 2) == 0) == b_even) {
                digit_set.insert(d);
            }
        }
    } else if (digit_mode == 2) {
        for (int d = 0; d <= 9; ++d) {
            digit_set.insert(d);
        }
    } else {
        throw std::runtime_error("digit_mode must be 0, 1, or 2");
    }

    return std::vector<int>(digit_set.begin(), digit_set.end());
}

static std::optional<FactorHit> build_hit_from_factor(
    const BigInt& N,
    BigInt g,
    const Candidate& candidate
) {
    if (g <= 1 || g >= N) {
        return std::nullopt;
    }

    if (N % g != 0) {
        return std::nullopt;
    }

    BigInt p = g;
    BigInt q = N / g;

    if (p > q) {
        std::swap(p, q);
    }

    BigInt a = (p + q) / 2;
    BigInt b = (q - p) / 2;

    if (a * a - b * b != N) {
        return std::nullopt;
    }

    BigInt Y = 2 * a * b;

    FactorHit hit;

    hit.p = p;
    hit.q = q;
    hit.a = a;
    hit.b = b;
    hit.Y = Y;

    hit.level_k = candidate.k;
    hit.gcd_coordinate = candidate.M_current;
    hit.hit_t = candidate.t;

    hit.has_stream_metadata = true;
    hit.hit_is_subtraction = candidate.is_subtraction;

    hit.seed_index = candidate.seed_index;
    hit.seed_a = candidate.seed_a;

    hit.hit_b_digit = candidate.b_digit;
    hit.hit_b_residue = candidate.b_residue;
    hit.hit_b_start = candidate.b_start;
    hit.hit_b_current = candidate.b_current;
    hit.hit_stride_b = candidate.stride_b;

    hit.true_b_mod_stride = mod_pos(b, candidate.stride_b);

    if (candidate.is_subtraction) {
        hit.collision_value = candidate.seed_a - candidate.b_current;
    } else {
        hit.collision_value = candidate.seed_a + candidate.b_current;
    }

    hit.abs_collision_value = abs_big(hit.collision_value);
    hit.proxy_gcd = gcd_int(hit.abs_collision_value, N);

    hit.gcd_result = g;
    hit.hit_shell_index = candidate.shell_index;
    hit.hit_chunk_index = candidate.chunk_index;
    hit.hit_chunk_attempt_start = candidate.chunk_attempt_start;
    hit.hit_shell_center_b = candidate.shell_center_b;
    hit.hit_shell_radius_b = candidate.shell_radius_b;

    if (g != 0 && candidate.M_current % g == 0) {
        hit.gcd_multiple = candidate.M_current / g;
    }

    return hit;
}

static std::optional<FactorHit> build_trivial_hit(
    const BigInt& N,
    const BigInt& g
) {
    if (g <= 1 || g >= N) {
        return std::nullopt;
    }

    if (N % g != 0) {
        return std::nullopt;
    }

    BigInt p = g;
    BigInt q = N / g;

    if (p > q) {
        std::swap(p, q);
    }

    BigInt a = (p + q) / 2;
    BigInt b = (q - p) / 2;
    BigInt Y = 2 * a * b;

    FactorHit hit;
    hit.p = p;
    hit.q = q;
    hit.a = a;
    hit.b = b;
    hit.Y = Y;
    hit.level_k = 0;
    hit.gcd_coordinate = g;
    hit.hit_t = 0;
    hit.has_stream_metadata = false;
    hit.gcd_result = g;

    return hit;
}

static std::optional<FactorHit> flush_batch(
    std::vector<Candidate>& batch,
    BigInt& product_mod_N,
    const BigInt& N,
    ThreadStats& stats
) {
    if (batch.empty()) {
        return std::nullopt;
    }

    ++stats.batch_gcd_checks;

    BigInt g = gcd_int(product_mod_N, N);

    if (g > 1) {
        for (const Candidate& candidate : batch) {
            ++stats.recovery_gcd_checks;

            BigInt cg = gcd_int(candidate.M_current, N);

            if (auto hit = build_hit_from_factor(N, cg, candidate)) {
                batch.clear();
                product_mod_N = 1;
                return hit;
            }
        }
    }

    batch.clear();
    product_mod_N = 1;

    return std::nullopt;
}

static std::vector<AcceleratorStream> build_streams_for_seed_range(
    const BigInt& N,
    const BigInt& b_cap,
    const RhoShell& rho_shell,
    const BigInt& base_seed,
    const BigInt& b_max,
    const std::vector<int>& b_digits,
    int max_k,
    std::uint64_t seed_start,
    std::uint64_t seed_begin,
    std::uint64_t seed_end,
    std::uint64_t chunk_index,
    std::uint64_t chunk_steps,
    bool centered,
    std::uint64_t base_window_steps,
    int schedule_mode,
    std::uint64_t low_k_max,
    std::uint64_t window_step_cap
) {
    std::vector<AcceleratorStream> streams;

    std::uint64_t local_seed_count = seed_end - seed_begin;
    std::size_t reserve_guess =
        static_cast<std::size_t>(local_seed_count) *
        static_cast<std::size_t>(max_k) *
        b_digits.size() *
        2;

    streams.reserve(reserve_guess);

    for (std::uint64_t local_seed_index = seed_begin;
         local_seed_index < seed_end;
         ++local_seed_index) {

        std::uint64_t global_seed_index = seed_start + local_seed_index;

        BigInt offset = signed_seed_offset(global_seed_index, centered);
        BigInt seed_a = base_seed + BigInt(2) * offset;

        for (int k = 1; k <= max_k; ++k) {
            const bool shell_mode = (schedule_mode == 2 || schedule_mode == 3);
            const BigInt& depth_for_window = shell_mode ? rho_shell.radius_b : b_cap;

            std::uint64_t one_side_window = adaptive_window_for_k(
                depth_for_window,
                k,
                base_window_steps,
                schedule_mode,
                low_k_max,
                window_step_cap
            );

            std::uint64_t full_stream_window = one_side_window;

            if (shell_mode) {
                full_stream_window = outward_attempt_window_from_radius_steps(
                    one_side_window,
                    window_step_cap
                );
            }

            std::uint64_t effective_chunk_steps = chunk_steps == 0 ? full_stream_window : chunk_steps;
            std::uint64_t chunk_attempt_start = 0;

            if (effective_chunk_steps == 0) {
                effective_chunk_steps = 1;
            }

            if (chunk_index > 0) {
                if (effective_chunk_steps > 0 &&
                    chunk_index > std::numeric_limits<std::uint64_t>::max() / effective_chunk_steps) {
                    continue;
                }
                chunk_attempt_start = chunk_index * effective_chunk_steps;
            }

            if (chunk_attempt_start >= full_stream_window) {
                continue;
            }

            std::uint64_t stream_window = std::min<std::uint64_t>(
                effective_chunk_steps,
                full_stream_window - chunk_attempt_start
            );

            if (stream_window == 0) {
                continue;
            }

            BigInt pow2 = pow2_big(k);

            BigInt A_num = (pow2 - 1) * N + seed_a;

            BigInt r_sub = mod_pos(A_num, pow2);
            BigInt r_add = mod_pos(-A_num, pow2);

            BigInt stride_b = 5 * pow2;
            BigInt stride_M = 5;

            for (int db : b_digits) {
                if (auto residue = crt_b_residue(db, r_sub, k)) {
                    BigInt b_start;
                    bool use_outward = shell_mode;

                    if (use_outward) {
                        b_start = nearest_with_residue_to_anchor(
                            rho_shell.center_b,
                            *residue,
                            stride_b,
                            BigInt(1),
                            b_max
                        );
                    } else {
                        // BOTTOM START:
                        // First matching CRT residue at or above b = 1.
                        b_start = first_geq_with_residue(BigInt(1), *residue, stride_b);
                    }

                    if (!use_outward && chunk_attempt_start > 0) {
                        BigInt shift = big_from_u64(chunk_attempt_start) * stride_b;
                        b_start += shift;
                    }

                    if (b_start <= b_max) {
                        BigInt M = div_pow2_floor(A_num - b_start, k);

                        streams.push_back({
                            global_seed_index,
                            seed_a,
                            k,
                            db,
                            b_start,
                            *residue,
                            b_start,
                            stride_b,
                            M,
                            stride_M,
                            A_num,
                            stream_window,
                            full_stream_window,
                            chunk_index,
                            chunk_attempt_start,
                            use_outward,
                            rho_shell.shell_index,
                            rho_shell.center_b,
                            rho_shell.radius_b,
                            true
                        });
                    }
                }

                if (auto residue = crt_b_residue(db, r_add, k)) {
                    BigInt b_start;
                    bool use_outward = shell_mode;

                    if (use_outward) {
                        b_start = nearest_with_residue_to_anchor(
                            rho_shell.center_b,
                            *residue,
                            stride_b,
                            BigInt(1),
                            b_max
                        );
                    } else {
                        // BOTTOM START:
                        // First matching CRT residue at or above b = 1.
                        b_start = first_geq_with_residue(BigInt(1), *residue, stride_b);
                    }

                    if (!use_outward && chunk_attempt_start > 0) {
                        BigInt shift = big_from_u64(chunk_attempt_start) * stride_b;
                        b_start += shift;
                    }

                    if (b_start <= b_max) {
                        BigInt M = div_pow2_floor(A_num + b_start, k);

                        streams.push_back({
                            global_seed_index,
                            seed_a,
                            k,
                            db,
                            b_start,
                            *residue,
                            b_start,
                            stride_b,
                            M,
                            stride_M,
                            A_num,
                            stream_window,
                            full_stream_window,
                            chunk_index,
                            chunk_attempt_start,
                            use_outward,
                            rho_shell.shell_index,
                            rho_shell.center_b,
                            rho_shell.radius_b,
                            false
                        });
                    }
                }
            }
        }
    }

    std::sort(
        streams.begin(),
        streams.end(),
        [](const AcceleratorStream& a, const AcceleratorStream& b) {
            if (a.k != b.k) {
                return a.k < b.k;
            }

            if (a.window_limit != b.window_limit) {
                return a.window_limit > b.window_limit;
            }

            if (a.seed_index != b.seed_index) {
                return a.seed_index < b.seed_index;
            }

            if (a.b_digit != b.b_digit) {
                return a.b_digit < b.b_digit;
            }

            return a.is_subtraction > b.is_subtraction;
        }
    );

    return streams;
}

static std::optional<FactorHit> parallel_multiseed_batch_accelerator(
    const BigInt& N,
    std::uint64_t base_window_steps,
    std::uint64_t seed_count,
    std::size_t batch_size,
    int max_k,
    unsigned int requested_threads,
    std::uint64_t seed_start,
    int digit_mode,
    bool centered,
    int schedule_mode,
    std::uint64_t low_k_max,
    int b_cap_mode,
    std::uint64_t target_divisor,
    const BigInt& manual_b_cap,
    std::uint64_t window_step_cap,
    std::uint64_t rho_center_den,
    std::uint64_t rho_radius_den,
    int rho_spectrum_mode,
    std::uint64_t chunk_steps,
    std::uint64_t chunk_start,
    std::uint64_t chunk_count,
    SearchStats& stats
) {
    if (N <= 1) {
        return std::nullopt;
    }

    if (N % 2 == 0) {
        return build_trivial_hit(N, BigInt(2));
    }

    if (N % 5 == 0 && N != 5) {
        return build_trivial_hit(N, BigInt(5));
    }

    BigInt sqrtN = isqrt_floor(N);

    if (is_perfect_square(N)) {
        return build_trivial_hit(N, sqrtN);
    }

    BigInt b_max = (N - 9) / 6;

    if (b_max < 1) {
        return std::nullopt;
    }

    BigInt base_seed = parity_aligned_seed(N);

    BigInt b_cap = compute_b_cap(
        sqrtN,
        b_max,
        b_cap_mode,
        target_divisor,
        manual_b_cap
    );

    if (b_cap > b_max) {
        b_cap = b_max;
    }

    if (rho_center_den == 0) {
        rho_center_den = 16;
    }

    if (rho_radius_den == 0) {
        rho_radius_den = 64;
    }

    std::uint64_t decimal_digits_N = static_cast<std::uint64_t>(N.get_str().size());
    std::uint64_t bit_length_N = static_cast<std::uint64_t>(mpz_sizeinbase(N.get_mpz_t(), 2));
    std::uint64_t maxbits_digits_N = maxbits_for_decimal_digits(decimal_digits_N);
    if (maxbits_digits_N == 0) {
        maxbits_digits_N = 1;
    }

    const std::uint64_t SAME_BIT_RHO_NUM = 292893ULL;
    const std::uint64_t SAME_BIT_RHO_DEN = 1000000ULL;
    std::uint64_t beta_rho_num = bit_length_N * SAME_BIT_RHO_NUM;
    std::uint64_t beta_rho_den = maxbits_digits_N * SAME_BIT_RHO_DEN;
    if (beta_rho_num == 0) {
        beta_rho_num = 1;
    }
    if (beta_rho_den == 0) {
        beta_rho_den = 1;
    }

    if (rho_spectrum_mode < 0 || rho_spectrum_mode > 2) {
        throw std::runtime_error("rho_spectrum_mode must be 0, 1, or 2");
    }

    std::vector<RhoShell> rho_shells = build_rho_shells(
        sqrtN,
        b_max,
        schedule_mode,
        rho_center_den,
        rho_radius_den,
        rho_spectrum_mode,
        bit_length_N,
        maxbits_digits_N
    );

    if (max_k <= 0) {
        std::size_t bits = mpz_sizeinbase(b_max.get_mpz_t(), 2);

        if (bits > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("N is too large for int-indexed k levels");
        }

        max_k = static_cast<int>(bits);
    }

    std::vector<int> b_digits = b_digits_for_mode(N, digit_mode);

    unsigned int hardware = std::thread::hardware_concurrency();
    unsigned int thread_count = requested_threads;

    if (thread_count == 0) {
        if (hardware <= 1) {
            thread_count = 1;
        } else {
            thread_count = hardware - 1;
        }
    }

    if (thread_count == 0) {
        thread_count = 1;
    }

    if (thread_count == 0) {
        thread_count = 1;
    }

    if (chunk_steps == 0) {
        chunk_steps = 1000000;
    }

    if (chunk_count == 0) {
        chunk_count = thread_count;
    }

    if (chunk_count == 0) {
        chunk_count = 1;
    }

    stats.seed_count = seed_count;
    stats.seed_start = seed_start;
    stats.max_k = max_k;
    stats.digit_mode = digit_mode;
    stats.centered = centered;
    stats.schedule_mode = schedule_mode;
    stats.low_k_max = low_k_max;
    stats.b_cap_mode = b_cap_mode;
    stats.target_divisor = target_divisor;
    stats.b_cap = b_cap;
    stats.manual_b_cap = manual_b_cap;
    stats.rho_center_den = rho_center_den;
    stats.rho_radius_den = rho_radius_den;
    stats.rho_spectrum_mode = rho_spectrum_mode;
    stats.decimal_digits_N = decimal_digits_N;
    stats.bit_length_N = bit_length_N;
    stats.maxbits_digits_N = maxbits_digits_N;
    stats.beta_rho_num = beta_rho_num;
    stats.beta_rho_den = beta_rho_den;
    stats.rho_shell_count = rho_shells.size();
    stats.chunk_steps = chunk_steps;
    stats.chunk_start = chunk_start;
    stats.chunk_count = chunk_count;
    stats.rho_center_b = rho_shells.empty() ? BigInt(0) : rho_shells.front().center_b;
    stats.rho_radius_b = rho_shells.empty() ? BigInt(0) : rho_shells.front().radius_b;
    stats.window_step_cap = window_step_cap;
    stats.thread_count = thread_count;

    std::atomic<bool> factor_found(false);
    std::mutex hit_mutex;
    std::optional<FactorHit> final_hit;

    // Seed-rho lane scheduler:
    // Work is split over the Cartesian product of rho shells and seeds.
    // Each lane is one (rho shell, local seed index). Threads pull lanes from
    // an atomic counter, so all cores are kept busy and rho shells are truly
    // parallelized across cores rather than being processed only as sequential waves.
    struct SeedRhoLane {
        std::size_t shell_pos = 0;
        std::uint64_t seed_begin = 0;
        std::uint64_t seed_end = 0;
        std::uint64_t chunk_index = 0;
    };

    std::vector<SeedRhoLane> lanes;
    lanes.reserve(static_cast<std::size_t>(seed_count) * rho_shells.size() * static_cast<std::size_t>(chunk_count));

    for (std::uint64_t c = 0; c < chunk_count; ++c) {
        std::uint64_t global_chunk_index = chunk_start + c;
        for (std::size_t shell_pos = 0; shell_pos < rho_shells.size(); ++shell_pos) {
            for (std::uint64_t s = 0; s < seed_count; ++s) {
                lanes.push_back({shell_pos, s, s + 1, global_chunk_index});
            }
        }
    }

    stats.total_lanes = lanes.size();

    std::atomic<std::size_t> next_lane(0);
    std::vector<ThreadStats> thread_stats(thread_count);
    std::vector<std::thread> threads;

    auto worker = [&](unsigned int thread_id) {
        ThreadStats local_stats;

        std::vector<Candidate> batch;
        batch.reserve(batch_size);

        BigInt product_mod_N = 1;

        auto publish_hit = [&](const FactorHit& hit) {
            bool expected = false;

            if (factor_found.compare_exchange_strong(expected, true)) {
                std::lock_guard<std::mutex> guard(hit_mutex);
                final_hit = hit;
            }
        };

        auto add_candidate_to_batch = [&](const Candidate& candidate) -> bool {
            BigInt m_mod = mod_pos(candidate.M_current, N);

            product_mod_N *= m_mod;
            product_mod_N %= N;

            batch.push_back(candidate);
            ++local_stats.proxy_candidates_tested;

            if (batch.size() >= batch_size) {
                if (auto hit = flush_batch(batch, product_mod_N, N, local_stats)) {
                    publish_hit(*hit);
                    return true;
                }
            }

            return false;
        };

        while (!factor_found.load()) {
            std::size_t lane_index = next_lane.fetch_add(1);

            if (lane_index >= lanes.size()) {
                break;
            }

            const SeedRhoLane& lane = lanes[lane_index];
            const RhoShell& shell = rho_shells[lane.shell_pos];

            std::vector<AcceleratorStream> streams = build_streams_for_seed_range(
                N,
                b_cap,
                shell,
                base_seed,
                b_max,
                b_digits,
                max_k,
                seed_start,
                lane.seed_begin,
                lane.seed_end,
                lane.chunk_index,
                chunk_steps,
                centered,
                base_window_steps,
                schedule_mode,
                low_k_max,
                window_step_cap
            );

            local_stats.active_streams += streams.size();

            for (auto& stream : streams) {
                if (factor_found.load()) {
                    break;
                }

                if (stream.use_outward) {
                    std::uint64_t attempt = stream.chunk_attempt_start;
                    std::uint64_t attempt_end = stream.chunk_attempt_start + stream.window_limit;

                    while (attempt < attempt_end && !factor_found.load()) {
                        BigInt b_value = outward_b_from_anchor(
                            stream.b_start,
                            stream.stride_b,
                            attempt
                        );

                        if (b_value < 1 || b_value > b_max) {
                            ++attempt;
                            continue;
                        }

                        BigInt M_value;
                        if (stream.is_subtraction) {
                            M_value = div_pow2_floor(stream.A_num - b_value, stream.k);
                        } else {
                            M_value = div_pow2_floor(stream.A_num + b_value, stream.k);
                        }

                        Candidate candidate;
                        candidate.t = attempt;
                        candidate.seed_index = stream.seed_index;
                        candidate.seed_a = stream.seed_a;
                        candidate.k = stream.k;
                        candidate.b_digit = stream.b_digit;
                        candidate.b_start = stream.b_start;
                        candidate.b_residue = stream.b_residue;
                        candidate.b_current = b_value;
                        candidate.stride_b = stream.stride_b;
                        candidate.M_current = M_value;
                        candidate.is_subtraction = stream.is_subtraction;
                        candidate.chunk_index = stream.chunk_index;
                        candidate.chunk_attempt_start = stream.chunk_attempt_start;
                        candidate.shell_index = stream.shell_index;
                        candidate.shell_center_b = stream.shell_center_b;
                        candidate.shell_radius_b = stream.shell_radius_b;

                        if (add_candidate_to_batch(candidate)) {
                            break;
                        }

                        ++attempt;
                    }
                } else {
                    for (std::uint64_t t = 0; t < stream.window_limit && !factor_found.load(); ++t) {
                        if (stream.b_current > b_max) {
                            break;
                        }

                        Candidate candidate;
                        candidate.t = stream.chunk_attempt_start + t;
                        candidate.seed_index = stream.seed_index;
                        candidate.seed_a = stream.seed_a;
                        candidate.k = stream.k;
                        candidate.b_digit = stream.b_digit;
                        candidate.b_start = stream.b_start;
                        candidate.b_residue = stream.b_residue;
                        candidate.b_current = stream.b_current;
                        candidate.stride_b = stream.stride_b;
                        candidate.M_current = stream.M_current;
                        candidate.is_subtraction = stream.is_subtraction;
                        candidate.chunk_index = stream.chunk_index;
                        candidate.chunk_attempt_start = stream.chunk_attempt_start;
                        candidate.shell_index = stream.shell_index;
                        candidate.shell_center_b = stream.shell_center_b;
                        candidate.shell_radius_b = stream.shell_radius_b;

                        if (add_candidate_to_batch(candidate)) {
                            break;
                        }

                        stream.b_current += stream.stride_b;

                        if (stream.is_subtraction) {
                            stream.M_current -= stream.stride_M;
                        } else {
                            stream.M_current += stream.stride_M;
                        }
                    }
                }
            }
        }

        if (!factor_found.load()) {
            if (auto hit = flush_batch(batch, product_mod_N, N, local_stats)) {
                publish_hit(*hit);
            }
        }

        thread_stats[thread_id] = local_stats;
    };

    for (unsigned int tid = 0; tid < thread_count; ++tid) {
        threads.emplace_back(worker, tid);
    }

    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }

    for (const auto& ts : thread_stats) {
        stats.proxy_candidates_tested += ts.proxy_candidates_tested;
        stats.batch_gcd_checks += ts.batch_gcd_checks;
        stats.recovery_gcd_checks += ts.recovery_gcd_checks;
        stats.active_streams += ts.active_streams;
    }

    std::uint64_t max_stream_window = 0;

    // Report the maximum adaptive window that could have been assigned.
    for (int k = 1; k <= max_k; ++k) {
        if (schedule_mode == 2 || schedule_mode == 3) {
            for (const auto& shell : rho_shells) {
                std::uint64_t one_side_w = adaptive_window_for_k(
                    shell.radius_b,
                    k,
                    base_window_steps,
                    schedule_mode,
                    low_k_max,
                    window_step_cap
                );

                std::uint64_t w = outward_attempt_window_from_radius_steps(
                    one_side_w,
                    window_step_cap
                );

                if (w > max_stream_window) {
                    max_stream_window = w;
                }
            }
        } else {
            std::uint64_t w = adaptive_window_for_k(
                b_cap,
                k,
                base_window_steps,
                schedule_mode,
                low_k_max,
                window_step_cap
            );

            if (w > max_stream_window) {
                max_stream_window = w;
            }
        }
    }

    stats.max_stream_window = max_stream_window;

    if (factor_found.load()) {
        std::lock_guard<std::mutex> guard(hit_mutex);
        return final_hit;
    }

    return std::nullopt;
}

static std::optional<FactorHit> factor_with_complex_lattice_method(
    const BigInt& N,
    std::uint64_t base_window_steps,
    std::uint64_t seed_count,
    std::size_t batch_size,
    int max_k,
    unsigned int thread_count,
    std::uint64_t seed_start,
    int digit_mode,
    bool centered,
    int schedule_mode,
    std::uint64_t low_k_max,
    int b_cap_mode,
    std::uint64_t target_divisor,
    const BigInt& manual_b_cap,
    std::uint64_t window_step_cap,
    std::uint64_t rho_center_den,
    std::uint64_t rho_radius_den,
    int rho_spectrum_mode,
    std::uint64_t chunk_steps,
    std::uint64_t chunk_start,
    std::uint64_t chunk_count,
    SearchStats& stats
) {
    auto start = std::chrono::steady_clock::now();

    auto hit = parallel_multiseed_batch_accelerator(
        N,
        base_window_steps,
        seed_count,
        batch_size,
        max_k,
        thread_count,
        seed_start,
        digit_mode,
        centered,
        schedule_mode,
        low_k_max,
        b_cap_mode,
        target_divisor,
        manual_b_cap,
        window_step_cap,
        rho_center_den,
        rho_radius_den,
        rho_spectrum_mode,
        chunk_steps,
        chunk_start,
        chunk_count,
        stats
    );

    auto stop = std::chrono::steady_clock::now();

    stats.runtime_ms = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
    );

    return hit;
}

static std::uint64_t parse_u64_arg(
    int argc,
    char** argv,
    int index,
    std::uint64_t default_value
) {
    if (argc <= index) {
        return default_value;
    }

    return static_cast<std::uint64_t>(std::stoull(argv[index]));
}

static int parse_int_arg(
    int argc,
    char** argv,
    int index,
    int default_value
) {
    if (argc <= index) {
        return default_value;
    }

    return std::stoi(argv[index]);
}

static BigInt parse_big_arg(
    int argc,
    char** argv,
    int index,
    const std::string& default_value
) {
    BigInt v;
    if (argc <= index) {
        v.set_str(default_value, 10);
        return v;
    }

    if (v.set_str(argv[index], 10) != 0) {
        throw std::runtime_error("invalid big integer argument");
    }

    return v;
}

int main(int argc, char** argv) {
    std::uint64_t base_window_steps = parse_u64_arg(argc, argv, 1, 10000);
    std::uint64_t seed_count        = parse_u64_arg(argc, argv, 2, 8);
    std::size_t batch_size          = static_cast<std::size_t>(parse_u64_arg(argc, argv, 3, 8192));
    int max_k                       = parse_int_arg(argc, argv, 4, 0);
    unsigned int thread_count       = static_cast<unsigned int>(parse_u64_arg(argc, argv, 5, 0));
    std::uint64_t seed_start        = parse_u64_arg(argc, argv, 6, 0);
    int digit_mode                  = parse_int_arg(argc, argv, 7, 1);
    bool centered                   = parse_int_arg(argc, argv, 8, 0) != 0;
    int schedule_mode               = parse_int_arg(argc, argv, 9, 3);
    std::uint64_t low_k_max         = parse_u64_arg(argc, argv, 10, 8);
    int b_cap_mode                  = parse_int_arg(argc, argv, 11, 1);
    std::uint64_t target_divisor    = parse_u64_arg(argc, argv, 12, 16);
    BigInt manual_b_cap             = parse_big_arg(argc, argv, 13, "0");
    std::uint64_t window_step_cap   = parse_u64_arg(argc, argv, 14, 0);
    std::uint64_t rho_center_den    = parse_u64_arg(argc, argv, 15, 16);
    std::uint64_t rho_radius_den    = parse_u64_arg(argc, argv, 16, 64);
    int rho_spectrum_mode           = parse_int_arg(argc, argv, 17, 0);
    std::uint64_t chunk_steps       = parse_u64_arg(argc, argv, 18, 1000000);
    std::uint64_t chunk_start       = parse_u64_arg(argc, argv, 19, 0);
    std::uint64_t chunk_count       = parse_u64_arg(argc, argv, 20, 0);

    if (base_window_steps == 0) {
        base_window_steps = 1;
    }

    if (seed_count == 0) {
        seed_count = 1;
    }

    if (batch_size == 0) {
        batch_size = 1;
    }

    if (chunk_steps == 0) {
        chunk_steps = 1000000;
    }

    if (digit_mode < 0 || digit_mode > 2) {
        std::cerr << "digit_mode must be 0, 1, or 2.\n";
        return 1;
    }

    if (schedule_mode < 0 || schedule_mode > 3) {
        std::cerr << "schedule_mode must be 0, 1, 2, or 3.\n";
        return 1;
    }

    if (b_cap_mode < 0 || b_cap_mode > 2) {
        std::cerr << "b_cap_mode must be 0, 1, or 2.\n";
        return 1;
    }

    if (low_k_max == 0) {
        low_k_max = 1;
    }

    if (target_divisor == 0) {
        target_divisor = 16;
    }

    if (rho_center_den == 0) {
        rho_center_den = 16;
    }

    if (rho_radius_den == 0) {
        rho_radius_den = 64;
    }

    if (rho_spectrum_mode < 0 || rho_spectrum_mode > 2) {
        std::cerr << "rho_spectrum_mode must be 0, 1, or 2.\n";
        return 1;
    }

    BigInt N;

    std::cout << "Complex Space Accelerator: PARALLEL STAGED BOTTOM-START multi-seed batch-GCD version\n";
    std::cout << "Settings:\n";
    std::cout << "base_window_steps = " << base_window_steps << "\n";
    std::cout << "seed_count = " << seed_count << "\n";
    std::cout << "batch_size = " << batch_size << "\n";
    std::cout << "max_k = " << max_k << " (0 means automatic)\n";
    std::cout << "thread_count = " << thread_count << " (0 means automatic)\n";
    std::cout << "seed_start = " << seed_start << "\n";
    std::cout << "vertical_start = bottom density-prior\n";
    std::cout << "digit_mode = " << digit_mode << " ";
    if (digit_mode == 0) {
        std::cout << "(Fermat digits)\n";
    } else if (digit_mode == 1) {
        std::cout << "(parity digits)\n";
    } else {
        std::cout << "(all digits)\n";
    }
    std::cout << "centered seeds = " << (centered ? "yes" : "no") << "\n";
    std::cout << "schedule_mode = " << schedule_mode << " ";
    if (schedule_mode == 0) {
        std::cout << "(fixed window)\n";
    } else if (schedule_mode == 1) {
        std::cout << "(adaptive low-k deepening)\n";
    } else if (schedule_mode == 2) {
        std::cout << "(single rho-shell outward low-k search)\n";
    } else {
        std::cout << "(rho-spectrum seed-rho tiled search)\n";
    }
    std::cout << "low_k_max = " << low_k_max << "\n";
    std::cout << "b_cap_mode = " << b_cap_mode << " ";
    if (b_cap_mode == 0) {
        std::cout << "(sqrt(N) / target_divisor)\n";
    } else if (b_cap_mode == 1) {
        std::cout << "(sqrt(b_max))\n";
    } else {
        std::cout << "(manual_b_cap)\n";
    }
    std::cout << "target_divisor = " << target_divisor << "\n";
    std::cout << "manual_b_cap = " << manual_b_cap << "\n";
    std::cout << "window_step_cap = " << window_step_cap << " (0 means uncapped)\n";
    std::cout << "rho_center_den = " << rho_center_den << "\n";
    std::cout << "rho_radius_den = " << rho_radius_den << "\n";
    std::cout << "rho_spectrum_mode = " << rho_spectrum_mode << " ";
    if (rho_spectrum_mode == 0) {
        std::cout << "(beta_N scaled balanced spectrum)\n";
    } else if (rho_spectrum_mode == 1) {
        std::cout << "(fixed balanced hard-shell spectrum)\n";
    } else {
        std::cout << "(general semiprime spectrum)\n";
    }
    std::cout << "chunk_steps = " << chunk_steps << "\n";
    std::cout << "chunk_start = " << chunk_start << "\n";
    std::cout << "chunk_count = " << chunk_count << " (0 means auto thread_count)\n";
    std::cout << "chunk_mapping = center, lower, upper, lower, upper\n";
    std::cout << "outward_window = 2 * one_side_radius_steps + 1\n\n";

    std::cout << "Enter odd composite N: ";
    std::cin >> N;

    if (N <= 1) {
        std::cout << "N must be greater than 1.\n";
        return 0;
    }

    SearchStats stats;

    auto hit = factor_with_complex_lattice_method(
        N,
        base_window_steps,
        seed_count,
        batch_size,
        max_k,
        thread_count,
        seed_start,
        digit_mode,
        centered,
        schedule_mode,
        low_k_max,
        b_cap_mode,
        target_divisor,
        manual_b_cap,
        window_step_cap,
        rho_center_den,
        rho_radius_den,
        rho_spectrum_mode,
        chunk_steps,
        chunk_start,
        chunk_count,
        stats
    );

    if (!hit) {
        std::cout << "\nNo factor found by accelerator window.\n";

        std::cout << "\nSearch diagnostics:\n";
        std::cout << "seed_start = " << stats.seed_start << "\n";
        std::cout << "seed_count = " << stats.seed_count << "\n";
        std::cout << "max_k = " << stats.max_k << "\n";
        std::cout << "digit_mode = " << stats.digit_mode << "\n";
        std::cout << "centered = " << (stats.centered ? "yes" : "no") << "\n";
        std::cout << "schedule_mode = " << stats.schedule_mode << "\n";
        std::cout << "low_k_max = " << stats.low_k_max << "\n";
        std::cout << "b_cap_mode = " << stats.b_cap_mode << "\n";
        std::cout << "target_divisor = " << stats.target_divisor << "\n";
        std::cout << "manual_b_cap = " << stats.manual_b_cap << "\n";
        std::cout << "b_cap = " << stats.b_cap << "\n";
        std::cout << "rho_center_den = " << stats.rho_center_den << "\n";
        std::cout << "rho_radius_den = " << stats.rho_radius_den << "\n";
        std::cout << "rho_spectrum_mode = " << stats.rho_spectrum_mode << "\n";
        std::cout << "decimal_digits_N = " << stats.decimal_digits_N << "\n";
        std::cout << "bit_length_N = " << stats.bit_length_N << "\n";
        std::cout << "maxbits_digits_N = " << stats.maxbits_digits_N << "\n";
        std::cout << "beta_rho_num = " << stats.beta_rho_num << "\n";
        std::cout << "beta_rho_den = " << stats.beta_rho_den << "\n";
        std::cout << "rho_shell_count = " << stats.rho_shell_count << "\n";
        std::cout << "chunk_steps = " << stats.chunk_steps << "\n";
        std::cout << "chunk_start = " << stats.chunk_start << "\n";
        std::cout << "chunk_count = " << stats.chunk_count << "\n";
        std::cout << "total_lanes = " << stats.total_lanes << "\n";
        std::cout << "rho_center_b = " << stats.rho_center_b << "\n";
        std::cout << "rho_radius_b = " << stats.rho_radius_b << "\n";
        std::cout << "window_step_cap = " << stats.window_step_cap << "\n";
        std::cout << "max_stream_window = " << stats.max_stream_window << "\n";
        std::cout << "thread_count = " << stats.thread_count << "\n";
        std::cout << "vertical_start = bottom density-prior\n";
        std::cout << "active accelerator streams = " << stats.active_streams << "\n";
        std::cout << "proxy candidates tested = " << stats.proxy_candidates_tested << "\n";
        std::cout << "batch GCD checks = " << stats.batch_gcd_checks << "\n";
        std::cout << "recovery GCD checks = " << stats.recovery_gcd_checks << "\n";
        std::cout << "runtime ms = " << stats.runtime_ms << "\n";

        return 0;
    }

    std::cout << "\nFactor found:\n";
    std::cout << "p = " << hit->p << "\n";
    std::cout << "q = " << hit->q << "\n";

    std::cout << "\nFermat lattice point:\n";
    std::cout << "a = " << hit->a << "\n";
    std::cout << "b = " << hit->b << "\n";
    std::cout << "a^2 - b^2 = " << N << "\n";

    std::cout << "\nSquared complex vertical-line point:\n";
    std::cout << "Z = " << N << " + " << hit->Y << "i\n";
    std::cout << "Y = 2ab = " << hit->Y << "\n";

    if (!hit->has_stream_metadata) {
        std::cout << "\nSearch path:\n";
        std::cout << "Trivial factor shortcut.\n";
    } else {
        std::cout << "\nSearch path:\n";
        std::cout << "Parallel staged bottom-start multi-seed iterated-average batch-GCD accelerator.\n";

        std::cout << "seed index = " << hit->seed_index << "\n";
        std::cout << "seed a_s = " << hit->seed_a << "\n";

        std::cout << "iterated average level k = " << hit->level_k << "\n";
        std::cout << "GCD coordinate M = " << hit->gcd_coordinate << "\n";
        std::cout << "hit offset t = " << hit->hit_t << "\n";
        std::cout << "rho shell index = " << hit->hit_shell_index << "\n";
        std::cout << "rho shell center b = " << hit->hit_shell_center_b << "\n";
        std::cout << "rho shell radius b = " << hit->hit_shell_radius_b << "\n";

        std::cout << "\nSuccessful accelerator stream:\n";
        std::cout << "branch = "
                  << (hit->hit_is_subtraction ? "subtraction" : "addition")
                  << "\n";

        std::cout << "b terminal digit = " << hit->hit_b_digit << "\n";
        std::cout << "CRT b residue = " << hit->hit_b_residue << "\n";
        std::cout << "b_start = " << hit->hit_b_start << "\n";
        std::cout << "b_current at hit = " << hit->hit_b_current << "\n";
        std::cout << "stride_b = " << hit->hit_stride_b << "\n";

        std::cout << "\nTrue Fermat b stream comparison:\n";
        std::cout << "true b = " << hit->b << "\n";
        std::cout << "true b mod stride_b = " << hit->true_b_mod_stride << "\n";

        std::cout << "true b lies in same CRT stream = "
                  << (hit->true_b_mod_stride == hit->hit_b_residue ? "yes" : "no")
                  << "\n";

        std::cout << "\nProxy collision diagnostic:\n";
        if (hit->hit_is_subtraction) {
            std::cout << "collision expression = seed_a - b_current\n";
        } else {
            std::cout << "collision expression = seed_a + b_current\n";
        }

        std::cout << "collision value = " << hit->collision_value << "\n";
        std::cout << "abs collision value = " << hit->abs_collision_value << "\n";
        std::cout << "gcd(abs collision value, N) = " << hit->proxy_gcd << "\n";

        std::cout << "\nGCD exposure:\n";
        std::cout << "gcd(M, N) = " << hit->gcd_result << "\n";

        if (hit->gcd_multiple != 0) {
            std::cout << "M / gcd(M, N) = " << hit->gcd_multiple << "\n";
        }
    }

    std::cout << "\nSearch diagnostics:\n";
    std::cout << "seed_start = " << stats.seed_start << "\n";
    std::cout << "seed_count = " << stats.seed_count << "\n";
    std::cout << "max_k = " << stats.max_k << "\n";
    std::cout << "digit_mode = " << stats.digit_mode << "\n";
    std::cout << "centered = " << (stats.centered ? "yes" : "no") << "\n";
    std::cout << "schedule_mode = " << stats.schedule_mode << "\n";
    std::cout << "low_k_max = " << stats.low_k_max << "\n";
    std::cout << "b_cap_mode = " << stats.b_cap_mode << "\n";
    std::cout << "target_divisor = " << stats.target_divisor << "\n";
    std::cout << "manual_b_cap = " << stats.manual_b_cap << "\n";
    std::cout << "b_cap = " << stats.b_cap << "\n";
    std::cout << "window_step_cap = " << stats.window_step_cap << "\n";
    std::cout << "max_stream_window = " << stats.max_stream_window << "\n";
    std::cout << "thread_count = " << stats.thread_count << "\n";
    std::cout << "vertical_start = bottom density-prior\n";
    std::cout << "active accelerator streams = " << stats.active_streams << "\n";
    std::cout << "proxy candidates tested = " << stats.proxy_candidates_tested << "\n";
    std::cout << "batch GCD checks = " << stats.batch_gcd_checks << "\n";
    std::cout << "recovery GCD checks = " << stats.recovery_gcd_checks << "\n";
    std::cout << "runtime ms = " << stats.runtime_ms << "\n";

    return 0;
}
