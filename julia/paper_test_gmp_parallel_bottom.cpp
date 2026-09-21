// C++17 GMP parallel multi-seed batch-GCD accelerator
// BOTTOM-STARTING vertical-line version.
//
// This version explicitly starts every CRT b-stream at the bottom of the
// complex vertical strip:
//
//     b_start = first b >= 1 with b == residue mod stride_b
//     b_t     = b_start + t * stride_b
//
// No midpoint or y-midpoint alignment is used.
//
// Compile in MSYS2 MINGW64:
// g++ -std=c++17 -O3 paper_test_gmp_parallel_bottom.cpp -lgmpxx -lgmp -pthread -o paper_test_gmp_parallel_bottom.exe
//
// Run:
// ./paper_test_gmp_parallel_bottom.exe
//
// Optional arguments:
// ./paper_test_gmp_parallel_bottom.exe <window_steps> <seed_count> <batch_size> <max_k> <thread_count> <seed_start> <digit_mode> <centered>
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
// max_k = 0 means automatic.
// thread_count = 0 means automatic.

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

    unsigned int thread_count = 1;

    std::uint64_t runtime_ms = 0;
};

static BigInt big_from_u64(std::uint64_t x) {
    BigInt v;
    v.set_str(std::to_string(x), 10);
    return v;
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

    // Correction (see Number Theory Papers/Complete_Fermat_Sieve_Verified.pdf).
    //
    // When 5 does not divide N, neither p nor q may end in 5 or 0, so both
    // must end in 1, 3, 7 or 9.  When 5 DOES divide N exactly one of them
    // carries the factor 5, and that restriction has to be lifted -- the
    // published 8-class table excludes N = 5 (mod 10) precisely because it
    // is structurally different (9 admissible classes per parity lane, 18
    // in total, rather than 4 and 8).
    //
    // Without this branch the routine returns an EMPTY pair set for
    // N = 5 (mod 10), which silently yields an empty b-digit sieve rather
    // than an error.  Callers here strip the factor 5 before reaching this
    // point, so the old behaviour was unreachable, but the function is now
    // total in its own right.
    const bool five_divides_N = (n10 == 5);

    auto valid_odd_factor_digit = [five_divides_N](int d) {
        d = ((d % 10) + 10) % 10;
        if (d % 2 == 0) return false;
        if (five_divides_N) return true;
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
    const BigInt& base_seed,
    const BigInt& b_max,
    const std::vector<int>& b_digits,
    int max_k,
    std::uint64_t seed_start,
    std::uint64_t seed_begin,
    std::uint64_t seed_end,
    bool centered
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
            BigInt pow2 = pow2_big(k);

            BigInt A_num = (pow2 - 1) * N + seed_a;

            BigInt r_sub = mod_pos(A_num, pow2);
            BigInt r_add = mod_pos(-A_num, pow2);

            BigInt stride_b = 5 * pow2;
            BigInt stride_M = 5;

            for (int db : b_digits) {
                if (auto residue = crt_b_residue(db, r_sub, k)) {
                    // BOTTOM START:
                    // First matching CRT residue at or above b = 1.
                    BigInt b_start = first_geq_with_residue(BigInt(1), *residue, stride_b);

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
                            true
                        });
                    }
                }

                if (auto residue = crt_b_residue(db, r_add, k)) {
                    // BOTTOM START:
                    // First matching CRT residue at or above b = 1.
                    BigInt b_start = first_geq_with_residue(BigInt(1), *residue, stride_b);

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
                            false
                        });
                    }
                }
            }
        }
    }

    return streams;
}

static std::optional<FactorHit> parallel_multiseed_batch_accelerator(
    const BigInt& N,
    std::uint64_t window_steps,
    std::uint64_t seed_count,
    std::size_t batch_size,
    int max_k,
    unsigned int requested_threads,
    std::uint64_t seed_start,
    int digit_mode,
    bool centered,
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
        thread_count = hardware == 0 ? 1 : hardware;
    }

    if (thread_count == 0) {
        thread_count = 1;
    }

    if (seed_count < thread_count) {
        thread_count = static_cast<unsigned int>(seed_count);
    }

    if (thread_count == 0) {
        thread_count = 1;
    }

    stats.seed_count = seed_count;
    stats.seed_start = seed_start;
    stats.max_k = max_k;
    stats.digit_mode = digit_mode;
    stats.centered = centered;
    stats.thread_count = thread_count;

    std::atomic<bool> factor_found(false);
    std::mutex hit_mutex;
    std::optional<FactorHit> final_hit;

    std::vector<ThreadStats> thread_stats(thread_count);
    std::vector<std::thread> threads;

    auto worker = [&](unsigned int thread_id,
                      std::uint64_t seed_begin,
                      std::uint64_t seed_end) {
        ThreadStats local_stats;

        std::vector<AcceleratorStream> streams = build_streams_for_seed_range(
            N,
            base_seed,
            b_max,
            b_digits,
            max_k,
            seed_start,
            seed_begin,
            seed_end,
            centered
        );

        local_stats.active_streams = streams.size();

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

        for (std::uint64_t t = 0; t < window_steps && !factor_found.load(); ++t) {
            bool any_target_active = false;

            for (auto& stream : streams) {
                if (factor_found.load()) {
                    break;
                }

                if (stream.b_current > b_max) {
                    continue;
                }

                any_target_active = true;

                Candidate candidate;
                candidate.t = t;
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

            if (!any_target_active) {
                break;
            }
        }

        if (!factor_found.load()) {
            if (auto hit = flush_batch(batch, product_mod_N, N, local_stats)) {
                publish_hit(*hit);
            }
        }

        thread_stats[thread_id] = local_stats;
    };

    std::uint64_t base_chunk = seed_count / thread_count;
    std::uint64_t remainder = seed_count % thread_count;

    std::uint64_t current = 0;

    for (unsigned int tid = 0; tid < thread_count; ++tid) {
        std::uint64_t chunk = base_chunk + (tid < remainder ? 1 : 0);
        std::uint64_t begin = current;
        std::uint64_t end = current + chunk;
        current = end;

        threads.emplace_back(worker, tid, begin, end);
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

    if (factor_found.load()) {
        std::lock_guard<std::mutex> guard(hit_mutex);
        return final_hit;
    }

    return std::nullopt;
}

static std::optional<FactorHit> factor_with_complex_lattice_method(
    const BigInt& N,
    std::uint64_t window_steps,
    std::uint64_t seed_count,
    std::size_t batch_size,
    int max_k,
    unsigned int thread_count,
    std::uint64_t seed_start,
    int digit_mode,
    bool centered,
    SearchStats& stats
) {
    auto start = std::chrono::steady_clock::now();

    auto hit = parallel_multiseed_batch_accelerator(
        N,
        window_steps,
        seed_count,
        batch_size,
        max_k,
        thread_count,
        seed_start,
        digit_mode,
        centered,
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

int main(int argc, char** argv) {
    std::uint64_t window_steps = parse_u64_arg(argc, argv, 1, 10000);
    std::uint64_t seed_count   = parse_u64_arg(argc, argv, 2, 8);
    std::size_t batch_size     = static_cast<std::size_t>(parse_u64_arg(argc, argv, 3, 8192));
    int max_k                  = parse_int_arg(argc, argv, 4, 0);
    unsigned int thread_count  = static_cast<unsigned int>(parse_u64_arg(argc, argv, 5, 0));
    std::uint64_t seed_start   = parse_u64_arg(argc, argv, 6, 0);
    int digit_mode             = parse_int_arg(argc, argv, 7, 0);
    bool centered              = parse_int_arg(argc, argv, 8, 0) != 0;

    if (window_steps == 0) {
        window_steps = 1;
    }

    if (seed_count == 0) {
        seed_count = 1;
    }

    if (batch_size == 0) {
        batch_size = 1;
    }

    if (digit_mode < 0 || digit_mode > 2) {
        std::cerr << "digit_mode must be 0, 1, or 2.\n";
        return 1;
    }

    BigInt N;

    std::cout << "Complex Space Accelerator: PARALLEL BOTTOM-START multi-seed batch-GCD version\n";
    std::cout << "Settings:\n";
    std::cout << "window_steps = " << window_steps << "\n";
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
    std::cout << "centered seeds = " << (centered ? "yes" : "no") << "\n\n";

    std::cout << "Enter odd composite N: ";
    std::cin >> N;

    if (N <= 1) {
        std::cout << "N must be greater than 1.\n";
        return 0;
    }

    SearchStats stats;

    auto hit = factor_with_complex_lattice_method(
        N,
        window_steps,
        seed_count,
        batch_size,
        max_k,
        thread_count,
        seed_start,
        digit_mode,
        centered,
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
        std::cout << "thread_count = " << stats.thread_count << "\n";
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
        std::cout << "Parallel bottom-start multi-seed iterated-average batch-GCD accelerator.\n";

        std::cout << "seed index = " << hit->seed_index << "\n";
        std::cout << "seed a_s = " << hit->seed_a << "\n";

        std::cout << "iterated average level k = " << hit->level_k << "\n";
        std::cout << "GCD coordinate M = " << hit->gcd_coordinate << "\n";
        std::cout << "hit offset t = " << hit->hit_t << "\n";

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
    std::cout << "thread_count = " << stats.thread_count << "\n";
    std::cout << "vertical_start = bottom density-prior\n";
    std::cout << "active accelerator streams = " << stats.active_streams << "\n";
    std::cout << "proxy candidates tested = " << stats.proxy_candidates_tested << "\n";
    std::cout << "batch GCD checks = " << stats.batch_gcd_checks << "\n";
    std::cout << "recovery GCD checks = " << stats.recovery_gcd_checks << "\n";
    std::cout << "runtime ms = " << stats.runtime_ms << "\n";

    return 0;
}
