// =====================================================================
//  complex_space_factorization.cpp
//
//  A unified factorization engine for the complex-space program:
//
//    * "Complex Space Factorization" (Viole)        -- strip geometry,
//      combined trial-division / vertical-Fermat / horizontal-Fermat
//      traversal (Sec. 4), the iterated-average GCD method (Sec. 5) and
//      its remainder sieve (Sec. 5.2).
//    * "The Complete Fermat Sieve by Terminal Digit -- a corrected and
//      verified reference"                          -- the corrected
//      terminal-digit sieve, both parity lanes, including the
//      N = 5 (mod 10) case that the 8-class table excludes.
//    * "Fermat Sieve Using Complex Numbers"         -- the original
//      per-lane (R,i) terminal-digit tables.
//    * Complex space.md                             -- the squared
//      complex space identity  sqrt(N^2 + (2ab)^2) = a^2 + b^2,
//      used here as an exact integer certificate.
//
//  ---------------------------------------------------------------
//  THE CENTRAL OBSERVATION THIS FILE IS BUILT ON
//  ---------------------------------------------------------------
//
//  Write N = a^2 - b^2 = (a-b)(a+b), so a = (p+q)/2 and b = (q-p)/2.
//  The terminal-digit sieve asks which (a mod 10, b mod 10) classes can
//  occur for a given N mod 10.  That is exactly the statement
//
//        a^2 - N  must be a square modulo 10,
//
//  together with the lane condition, which is "a^2 - N is a square
//  modulo 4".  In other words the corrected terminal-digit sieve is
//  precisely the quadratic-residue sieve at the single modulus 20.
//
//  Once seen that way it generalises: for ANY modulus m,
//
//        a admissible  <=>  (a^2 - N) mod m is a square mod m.
//
//  Measured survival fractions (see --selftest):
//
//        modulus 10 alone          ->  2.5x reduction
//        wheel 2^6*3^3*5^2*7*11*13*17
//                                  ->  ~4,000x reduction
//
//  So the digit sieve, generalised to a CRT wheel, is worth three to
//  four orders of magnitude rather than a factor of 2.5.  Everything
//  below is organised around that wheel.
//
//  ---------------------------------------------------------------
//  ALGORITHM
//  ---------------------------------------------------------------
//
//  Stage 0   Reduction: powers of 2, perfect powers, wheel-30 trial
//            division to B1, strong probable-prime test.
//
//  Stage 1   A race of deterministic streams over all worker threads,
//            stopping at the first split.  This is the paper's Sec. 4
//            "combined method", made concrete:
//
//              TD   trial division climbing from B1        (top of strip)
//              VF   wheel-sieved vertical Fermat from      (bottom of strip)
//                   ceil(sqrt(N)) upward
//              HF   wheel-sieved horizontal Fermat         (middle of strip,
//                   from b = 1 upward                       optional, see below)
//              IA   iterated-average GCD accelerator       (Sec. 5 / 5.2)
//                   with the every-5th-integer sieve        (optional)
//              LM   Lehman multiplier sweep: the same       (completeness)
//                   wheel-sieved scan applied to 4kN
//
//  Stage 2   Recurse on both cofactors until everything is prime.
//
//  Worst-case cost is O(N^(1/3)) arithmetic operations: TD covers every
//  prime below N^(1/3) and LM covers every remaining case, which is
//  Lehman's theorem.  Racing TD against VF alone would only give
//  O(N^(3/8)).  The wheel reduces the constant on VF/LM by ~4,000x.
//
//  Fidelity note: HF and IA are implemented faithfully to the paper and
//  are available via flags, but they are OFF by default because VF
//  strictly dominates them.  With p,q = sqrt(N)(1 -+ e), the vertical
//  scan length is j ~ b^2 / (2 sqrt(N)), i.e. quadratically shorter than
//  the horizontal scan length b.  See the benchmark table in
//  Prime Factorization/Unified Complex Space Factorization.md.
//
//  ---------------------------------------------------------------
//  BUILD
//  ---------------------------------------------------------------
//    g++ -std=c++17 -O3 complex_space_factorization.cpp -lgmpxx -lgmp -pthread -o csf
//
//  MSYS2 / MinGW-w64 is a supported target.  uint64_t there is
//  unsigned long long and unsigned long is only 32 bits, so every 64-bit
//  value that meets an mpz goes through big_from_u64 / u64_from_big /
//  divisible_by_u64 rather than gmpxx's unsigned long overloads.
//
//  RUN THIS BEFORE COMMITTING any code that puts a uint64_t near a BigInt.
//  It needs no Windows box:
//
//    sed 's/std::uint64_t/unsigned long long/g' complex_space_factorization.cpp > probe.cpp
//    g++ -std=c++17 -fsyntax-only probe.cpp        # must report 0 errors
//
//  On LP64 unsigned long long is a distinct type from unsigned long, so
//  that substitution reproduces the MinGW overload set exactly.  It is not
//  advisory: ctm_stream shipped with seven ambiguous-overload sites --
//  BigInt(CH * STEP), BigInt(chunk_index) and three `R += STEP` in the mpz
//  walk -- which built clean on Linux, broke every MinGW build, and which
//  this probe reproduces exactly, at the same seven lines.  The trap is
//  that LP64 hides the whole class of bug, so a clean local build proves
//  nothing here.
//
//  RUN
//    ./csf <N> [options]
//    ./csf --selftest
//
//  Options:
//    --threads=T     worker threads (0 = hardware concurrency)
//    --b1=B          stage-0 trial division bound      (default 1000000)
//    --no-lehman     disable the Lehman multiplier stream
//    --hf            enable the horizontal (b-driven) Fermat stream
//    --ia            enable the iterated-average GCD stream (Sec. 5/5.2)
//    --digit-only    restrict the sieve to modulus 20, i.e. the corrected
//                    terminal-digit sieve alone (for comparison)
//    --no-sieve      disable sieving entirely (for comparison)
//    --quiet         print only the factorization
//    --verbose       print stream statistics
// =====================================================================

#include <gmpxx.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using BigInt = mpz_class;

// ---------------------------------------------------------------------
// Basic helpers
// ---------------------------------------------------------------------

static inline BigInt isqrt_floor(const BigInt& n) {
    BigInt r;
    mpz_sqrt(r.get_mpz_t(), n.get_mpz_t());
    return r;
}

static inline BigInt isqrt_ceil(const BigInt& n) {
    BigInt r = isqrt_floor(n);
    if (r * r < n) ++r;
    return r;
}

static inline bool is_perfect_square(const BigInt& x, BigInt* root = nullptr) {
    if (x < 0) return false;
    if (mpz_perfect_square_p(x.get_mpz_t()) == 0) return false;
    if (root) *root = isqrt_floor(x);
    return true;
}

static inline bool is_probable_prime(const BigInt& n) {
    return mpz_probab_prime_p(n.get_mpz_t(), 30) != 0;
}

static inline std::uint64_t umod(const BigInt& n, std::uint64_t m) {
    // Every modulus used here is < 2^32 (build_wheel caps W), so the cast to
    // unsigned long is exact on LLP64 as well as LP64.
    return static_cast<std::uint64_t>(
        mpz_fdiv_ui(n.get_mpz_t(), static_cast<unsigned long>(m)));
}

// ---------------------------------------------------------------------
// Portable 64-bit <-> mpz helpers.
//
// On LP64 (Linux, macOS) uint64_t IS unsigned long, so gmpxx's unsigned
// long overloads cover everything and plain BigInt(x) compiles.  On LLP64
// (MinGW-w64, MSVC) unsigned long is 32 bits and uint64_t is unsigned long
// long, which gmpxx has no constructor for: BigInt(x) is then ambiguous,
// and mpz_*_ui / mpz_get_ui silently truncate to 32 bits.
//
// These three helpers are used everywhere a 64-bit value meets an mpz, so
// the same source is correct on both models.
// ---------------------------------------------------------------------

static inline BigInt big_from_u64(std::uint64_t v) {
    BigInt r;
    mpz_import(r.get_mpz_t(), 1, -1, sizeof(v), 0, 0, &v);
    return r;
}

// Saturating export: returns `cap` when z does not fit in 64 bits.
static inline std::uint64_t u64_from_big(const BigInt& z, std::uint64_t cap) {
    if (mpz_sgn(z.get_mpz_t()) <= 0) return 0;
    if (mpz_sizeinbase(z.get_mpz_t(), 2) > 64) return cap;
    std::uint64_t v = 0;
    std::size_t words = 0;
    mpz_export(&v, &words, -1, sizeof(v), 0, 0, z.get_mpz_t());
    return words ? v : 0;
}

// Exact export of a value known to fit 128 bits.
static inline __uint128_t u128_from_big(const BigInt& z) {
    __uint128_t v = 0;
    std::size_t words = 0;
    mpz_export(&v, &words, -1, sizeof(v), 0, 0, z.get_mpz_t());
    return v;
}

// Least non-negative representative of x mod m, for small positive m.
static inline BigInt mod_pos_small(const BigInt& x, long m) {
    BigInt r = x % m;
    if (mpz_sgn(r.get_mpz_t()) < 0) r += m;
    return r;
}

static inline bool divisible_by_u64(const BigInt& n, std::uint64_t d) {
    if (d <= static_cast<std::uint64_t>(std::numeric_limits<unsigned long>::max()))
        return mpz_divisible_ui_p(n.get_mpz_t(), static_cast<unsigned long>(d)) != 0;
    return mpz_divisible_p(n.get_mpz_t(), big_from_u64(d).get_mpz_t()) != 0;
}

// Set of squares modulo m.  For a perfect square x, x mod m must lie in
// this set; that is the whole content of the sieve.
static std::vector<std::uint8_t> squares_mod(std::uint64_t m) {
    std::vector<std::uint8_t> s(m, 0);
    for (std::uint64_t r = 0; r < m; ++r) s[(r * r) % m] = 1;
    return s;
}

// ---------------------------------------------------------------------
// Corrected complete Fermat terminal-digit sieve
//
//   N = R^2 - i^2,  p = R - i,  q = R + i.
//   Lane A: N = 1 (mod 4)  <=>  R odd,  i even
//   Lane B: N = 3 (mod 4)  <=>  R even, i odd
//
// Indexed by N mod 10 alone the admissible (R,i) mod 10 classes number 8
// for N ending in 1,3,7,9 (four from each lane) and 18 for N ending in 5.
// Because we always know N mod 4 we use the lane directly: 4 classes for
// N coprime to 10, 9 classes for 5 | N.
//
// The 5 | N branch is the one the published 8-class table excludes.  It
// is included here so the routine is total: any caller that has not
// already stripped the factor 5 still gets a sound (not empty) sieve.
// ---------------------------------------------------------------------

struct DigitPair { int a10; int b10; };

static std::vector<DigitPair> fermat_digit_pairs(const BigInt& N) {
    const int n10 = static_cast<int>(umod(N, 10));
    const int n4  = static_cast<int>(umod(N, 4));

    bool a_even;
    if (n4 == 3)      a_even = true;   // Lane B
    else if (n4 == 1) a_even = false;  // Lane A
    else throw std::runtime_error("fermat_digit_pairs: N must be odd");

    const bool b_even = !a_even;
    const bool five_divides_N = (n10 == 5);

    std::vector<DigitPair> pairs;
    for (int da = 0; da <= 9; ++da) {
        if (((da % 2) == 0) != a_even) continue;
        for (int db = 0; db <= 9; ++db) {
            if (((db % 2) == 0) != b_even) continue;

            const int p10 = ((da - db) % 10 + 10) % 10;
            const int q10 = (da + db) % 10;

            // p and q are odd.  When 5 does not divide N neither p nor q
            // may end in 5 or 0, so both must end in 1, 3, 7 or 9.  When
            // 5 does divide N exactly one of them carries the 5, so that
            // restriction must be lifted -- this is the correction.
            auto odd_digit   = [](int d) { return d % 2 == 1; };
            auto coprime10   = [](int d) { return d == 1 || d == 3 || d == 7 || d == 9; };

            if (!odd_digit(p10) || !odd_digit(q10)) continue;
            if (!five_divides_N && (!coprime10(p10) || !coprime10(q10))) continue;

            if ((p10 * q10) % 10 == n10) pairs.push_back({da, db});
        }
    }
    return pairs;
}

// ---------------------------------------------------------------------
// Generalised quadratic-residue wheel
//
// The terminal-digit sieve is the m = 20 case of:
//
//     a is admissible  <=>  (a^2 - M) mod m is a square mod m
//
// We build a CRT wheel from several coprime prime powers and materialise
// the sorted list of admissible residues modulo W = prod(m_i).  The list
// is built by CRT composition, so construction costs O(output), not O(W).
//
// A second tier of small prime moduli is applied per-candidate; those are
// cheaper to test than to fold into W.
// ---------------------------------------------------------------------

struct Wheel {
    std::uint64_t W = 1;
    std::vector<std::uint32_t> residues;        // admissible a mod W, sorted
    std::vector<std::uint32_t> sec_mod;         // second-tier moduli
    std::vector<std::vector<std::uint8_t>> sec_ok;

    // Division-free walk of the residue list.
    //
    // The secondary test is ok[s][(base_mod[s] + residues[i]) % p_s].  Doing
    // that modulo per residue per modulus was the whole cost of the scan:
    // an integer division is ~20-40 cycles against ~2 for an add, and it is
    // far worse on a GPU, where integer division is emulated.
    //
    // residues[] is fixed for a given N, so the per-step increment
    //     step[s][i] = (residues[i] - residues[i-1]) mod p_s
    // can be precomputed once and amortised over every block.  Walking then
    // costs one add and one conditional subtract: since step < p and
    // cur < p, the sum is below 2p and a single subtract renormalises.
    //
    // A chunk starting at an arbitrary index pays one division to seed
    // cur[s], then none for the rest of the chunk.  This is also exactly the
    // shape a GPU kernel wants: a coalesced byte read and an add.
    // sec_res[i * sec_mod.size() + s] = residues[i] mod sec_mod[s].
    //
    // Absolute, not incremental. Since base_mod[s] < p and sec_res < p, the
    // live residue is (base_mod[s] + sec_res) with one conditional subtract
    // -- no division anywhere in the scan, and no loop-carried dependency,
    // so residues are independent and the loop vectorises.
    //
    // Interleaved so one 64-byte line carries every modulus for several
    // consecutive residues. Independence + coalesced reads + the branchless
    // fold below are exactly the three properties a GPU kernel needs; the
    // CUDA port in complex_space_factorization_cuda.cu reuses this table
    // verbatim.
    std::vector<std::uint8_t> sec_res;

    double density = 1.0;

    bool empty() const { return residues.empty(); }
};

// Admissible residues of a modulo m for a^2 - M = square.
static std::vector<std::uint32_t> admissible_mod(const BigInt& M, std::uint64_t m) {
    const std::vector<std::uint8_t> sq = squares_mod(m);
    const std::uint64_t Mm = umod(M, m);
    std::vector<std::uint32_t> out;
    for (std::uint64_t r = 0; r < m; ++r) {
        const std::uint64_t v = ((r * r) % m + m - Mm) % m;
        if (sq[v]) out.push_back(static_cast<std::uint32_t>(r));
    }
    return out;
}

// Extended gcd based CRT for two coprime moduli.
static std::uint64_t crt_pair(std::uint64_t r1, std::uint64_t m1,
                              std::uint64_t r2, std::uint64_t m2) {
    // Return the unique x in [0, m1*m2) with x = r1 (mod m1), x = r2 (mod m2).
    // Requires gcd(m1, m2) = 1.
    const long long m = static_cast<long long>(m2);

    // m1^{-1} (mod m2) by extended Euclid.
    long long old_r = static_cast<long long>(m1 % m2), r = m;
    long long old_s = 1, sc = 0;
    while (r != 0) {
        const long long q = old_r / r;
        long long tmp = old_r - q * r; old_r = r; r = tmp;
        tmp = old_s - q * sc; old_s = sc; sc = tmp;
    }
    const long long inv = ((old_s % m) + m) % m;

    const long long diff = static_cast<long long>((r2 + m2 - r1 % m2) % m2);
    long long t = static_cast<long long>((static_cast<__int128>(diff) * inv) % m);
    if (t < 0) t += m;
    return r1 + m1 * static_cast<std::uint64_t>(t);
}

static Wheel build_wheel(const BigInt& M,
                         std::uint64_t residue_budget,
                         bool digit_only,
                         bool disabled) {
    Wheel w;
    if (disabled) {
        w.W = 1;
        w.residues = {0};
        w.density = 1.0;
        return w;
    }

    // Primary moduli, most selective first.  2^k and 5^k subsume the
    // corrected terminal-digit sieve (modulus 20 = 4 * 5).
    std::vector<std::uint64_t> primary;
    if (digit_only) {
        primary = {4, 5};                        // exactly the digit sieve
    } else {
        primary = {64, 27, 25, 7, 11, 13, 17};
    }

    std::vector<std::uint32_t> cur = {0};
    std::uint64_t curW = 1;

    for (std::uint64_t m : primary) {
        std::vector<std::uint32_t> am = admissible_mod(M, m);
        if (am.empty()) { w.W = curW; w.residues.clear(); return w; }

        // Stop growing if the composed list would exceed the budget or
        // overflow 32-bit residues.
        // Two independent caps. residue_budget bounds the materialised list;
        // 0xFFFFFFFF bounds W itself, because residues are stored as uint32.
        // The modulus cap binds first with the default moduli, so raising
        // --wheel past ~2^18 has no effect -- W is already maximal.
        const __uint128_t newW = static_cast<__uint128_t>(curW) * m;
        const __uint128_t newN = static_cast<__uint128_t>(cur.size()) * am.size();
        if (newW > 0xFFFFFFFFull || newN > residue_budget) break;

        std::vector<std::uint32_t> next;
        next.reserve(static_cast<std::size_t>(newN));
        for (std::uint32_t r1 : cur)
            for (std::uint32_t r2 : am)
                next.push_back(static_cast<std::uint32_t>(crt_pair(r1, curW, r2, m)));

        cur.swap(next);
        curW = static_cast<std::uint64_t>(newW);
    }

    std::sort(cur.begin(), cur.end());
    w.W = curW;
    w.residues = std::move(cur);
    w.density = static_cast<double>(w.residues.size()) / static_cast<double>(w.W);

    if (!digit_only) {
        // Second tier: applied per surviving candidate.
        for (std::uint32_t p : {19u, 23u, 29u, 31u, 37u, 41u, 43u, 47u}) {
            if (w.W % p == 0) continue;
            std::vector<std::uint8_t> ok(p, 0);
            const std::vector<std::uint8_t> sq = squares_mod(p);
            const std::uint64_t Mp = umod(M, p);
            std::size_t cnt = 0;
            for (std::uint32_t r = 0; r < p; ++r) {
                const std::uint64_t v = ((1ull * r * r) % p + p - Mp) % p;
                if (sq[v]) { ok[r] = 1; ++cnt; }
            }
            w.sec_mod.push_back(p);
            w.sec_ok.push_back(std::move(ok));
            w.density *= static_cast<double>(cnt) / static_cast<double>(p);
        }

        // Precompute the per-residue increments (see Wheel::sec_step).
        const std::size_t R = w.residues.size();
        const std::size_t NS = w.sec_mod.size();
        w.sec_res.assign(R * NS, 0);
        for (std::size_t s2 = 0; s2 < NS; ++s2) {
            const std::uint32_t p = w.sec_mod[s2];
            for (std::size_t i = 0; i < R; ++i)
                w.sec_res[i * NS + s2] = static_cast<std::uint8_t>(w.residues[i] % p);
        }
    }
    return w;
}

// ---------------------------------------------------------------------
// Squared complex space certificate (Complex space.md, footnote 1)
//
// A complex factor (a + bi) squares to n + 2abi, and the right triangle
// it represents gives the exact integer identity
//
//        sqrt(N^2 + (2ab)^2) = a^2 + b^2.
//
// We use it as an independent check on every split the engine reports:
// it must hold exactly, in integers, or the split is rejected.
// ---------------------------------------------------------------------

static bool certify_complex_square(const BigInt& N, const BigInt& a, const BigInt& b) {
    if (a <= b || b < 0) return false;
    if ((a - b) * (a + b) != N) return false;

    const BigInt two_ab = 2 * a * b;
    const BigInt hyp2   = N * N + two_ab * two_ab;
    BigInt hyp;
    if (!is_perfect_square(hyp2, &hyp)) return false;
    return hyp == a * a + b * b;
}

// ---------------------------------------------------------------------
// Shared result slot
// ---------------------------------------------------------------------

struct Found {
    std::atomic<bool> done{false};
    std::mutex mu;
    BigInt factor;
    std::string source;
    BigInt cert_a, cert_b;
    bool certified = false;

    bool ready() const { return done.load(std::memory_order_relaxed); }

    void submit(const BigInt& N, const BigInt& f, const char* src,
                const BigInt* a = nullptr, const BigInt* b = nullptr) {
        if (f <= 1 || f >= N) return;
        if (N % f != 0) return;
        std::lock_guard<std::mutex> lk(mu);
        if (done.load()) return;
        factor = f;
        source = src;
        if (a && b) {
            cert_a = *a; cert_b = *b;
            certified = certify_complex_square(N, *a, *b);
        }
        done.store(true, std::memory_order_release);
    }
};

struct Stats {
    std::atomic<std::uint64_t> vf_candidates{0};   // survived the wheel
    std::atomic<std::uint64_t> vf_scanned{0};      // residues visited
    std::atomic<std::uint64_t> td_tested{0};
    std::atomic<std::uint64_t> hf_candidates{0};
    std::atomic<std::uint64_t> ia_gcds{0};
    std::atomic<std::uint64_t> yp_steps{0};
    std::atomic<std::uint64_t> ctm_steps{0};
    std::atomic<std::uint64_t> lm_candidates{0};
};

// ---------------------------------------------------------------------
// Generic wheel-sieved Fermat scan
//
// Finds x in [x_lo, x_hi] with x^2 - M a perfect square, then offers
// gcd(|x-y|, N) and gcd(x+y, N) as splits of N.
//
//   vertical Fermat    M = N,     x = a, y = b
//   horizontal Fermat  M = -N,    x = b, y = a
//   Lehman multiplier  M = 4kN,   x = a, y = b
//
// All three are the same scan; only M and the interval differ.  This is
// the paper's Sec. 4 "run every method simultaneously" made literal --
// the methods are one routine pointed at different parts of the strip.
// ---------------------------------------------------------------------

// Chunked, work-stealing variant for the unbounded vertical/horizontal
// streams: residue indices are handed out by an atomic counter so that
// low x (short scans, the balanced-semiprime case) is always covered
// first regardless of thread count.
static void fermat_scan_chunked(const BigInt& N,
                                const BigInt& M,
                                const Wheel& w,
                                const BigInt& x_lo,
                                const BigInt& x_hi,
                                Found& found,
                                std::atomic<std::uint64_t>& next_chunk,
                                std::atomic<std::uint64_t>& cand_counter,
                                std::atomic<std::uint64_t>& scan_counter,
                                const char* src) {
    if (w.empty() || x_hi < x_lo) return;

    const std::uint64_t W = w.W;
    const std::uint64_t R = w.residues.size();
    const std::uint64_t CHUNK = std::min<std::uint64_t>(R, 4096);
    const std::uint64_t chunks_per_block = (R + CHUNK - 1) / CHUNK;

    const BigInt base0 = x_lo - big_from_u64(umod(x_lo, W));

    BigInt x, t, y, g, base;
    std::vector<std::uint64_t> base_mod(w.sec_mod.size());
    std::uint64_t cur_blk = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t local_cand = 0, local_scan = 0;
    bool exhausted = false;

    while (!found.ready() && !exhausted) {
        const std::uint64_t c = next_chunk.fetch_add(1, std::memory_order_relaxed);
        const std::uint64_t blk = c / chunks_per_block;
        const std::uint64_t k   = c % chunks_per_block;

        if (blk != cur_blk) {
            cur_blk = blk;
            base = base0 + big_from_u64(blk) * big_from_u64(W);
            if (base > x_hi) break;
            for (std::size_t i = 0; i < w.sec_mod.size(); ++i)
                base_mod[i] = umod(base, w.sec_mod[i]);
        }

        const std::uint64_t i0 = k * CHUNK;
        const std::uint64_t i1 = std::min(R, i0 + CHUNK);

        // Bounds as base-relative OFFSETS, so the hot loop never touches an
        // mpz.  base is the largest multiple of W at or below x_lo, so
        // x_lo - base < W < 2^32; and when x_hi - base exceeds W no residue
        // in this block can reach it.
        std::uint64_t lo_off = 0, hi_off = std::numeric_limits<std::uint64_t>::max();
        if (base < x_lo) lo_off = u64_from_big(x_lo - base, W);
        {
            const BigInt span = x_hi - base;
            if (mpz_sgn(span.get_mpz_t()) < 0) { exhausted = true; break; }
            hi_off = u64_from_big(span, std::numeric_limits<std::uint64_t>::max());
        }

        // Seed the division-free walk: one modulo per modulus per chunk,
        // then add-and-conditionally-subtract for every residue after that.
        const std::size_t NS = w.sec_mod.size();
        std::uint32_t smod[16], bmod[16];
        const std::uint8_t* okp[16];
        for (std::size_t s = 0; s < NS; ++s) {
            smod[s] = w.sec_mod[s];
            okp[s]  = w.sec_ok[s].data();
            bmod[s] = static_cast<std::uint32_t>(base_mod[s] % smod[s]);
        }

        for (std::uint64_t i = i0; i < i1; ++i) {
            const std::uint32_t r = w.residues[i];
            if (r < lo_off) continue;
            if (r > hi_off) { exhausted = true; break; }

            ++local_scan;

            // Branchless fold. The early-exit form tested ~2 moduli per
            // residue, but each test is a coin-flip branch and a mispredict
            // costs more than the check itself. Folding all NS lookups with
            // & leaves one well-predicted branch (survivors are ~1 in 300).
            // On a GPU the same rewrite removes warp divergence.
            const std::uint8_t* rm = &w.sec_res[i * NS];
            unsigned okm = 1u;
            for (std::size_t s = 0; s < NS; ++s) {
                std::uint32_t c = bmod[s] + rm[s];
                const std::uint32_t p = smod[s];
                c -= (c >= p ? p : 0);
                okm &= okp[s][c];
            }
            if (!okm) continue;

            // Survivor: only now is an mpz built.  This is the host/device
            // split -- everything above is machine words and maps directly
            // onto a GPU kernel that emits surviving offsets.
            ++local_cand;
            x = base + r;
            t = x * x - M;
            if (t < 0) continue;
            if (!is_perfect_square(t, &y)) continue;

            BigInt d = x > y ? x - y : y - x;
            mpz_gcd(g.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());
            if (g > 1 && g < N) {
                if (M == N)       found.submit(N, g, src, &x, &y);
                else if (M == -N) found.submit(N, g, src, &y, &x);
                else              found.submit(N, g, src);
                break;
            }
            BigInt s2 = x + y;
            mpz_gcd(g.get_mpz_t(), s2.get_mpz_t(), N.get_mpz_t());
            if (g > 1 && g < N) { found.submit(N, g, src); break; }
        }
        if ((local_scan & 0xFFFF) == 0 && found.ready()) break;
    }
    scan_counter += local_scan;
    cand_counter += local_cand;
}

// ---------------------------------------------------------------------
// Trial division stream -- the top of the factor strip.
// Wheel-30 residues, climbing from `start` to `limit`.
// ---------------------------------------------------------------------

// Trial division over p in [start, limit], chunked by an atomic counter so a
// thread that exhausts its own arc of the strip can migrate here rather than
// idle. Wheel-30 residues.
static void trial_division_stream(const BigInt& N,
                                  std::uint64_t start,
                                  std::uint64_t limit,
                                  Found& found,
                                  std::atomic<std::uint64_t>& next_block,
                                  std::atomic<std::uint64_t>& counter) {
    static const std::uint64_t W30[8] = {1, 7, 11, 13, 17, 19, 23, 29};
    const std::uint64_t TURNS = 1u << 14;
    const std::uint64_t base0 = (start / 30) * 30;
    std::uint64_t tested = 0;

    while (!found.ready()) {
        const std::uint64_t blk = next_block.fetch_add(1, std::memory_order_relaxed);
        if (blk > limit / (TURNS * 30) + 1) break;
        const std::uint64_t lo = base0 + blk * TURNS * 30;
        if (lo > limit) break;
        const std::uint64_t hi = std::min(limit, lo + TURNS * 30 - 1);

        for (std::uint64_t base = lo; base <= hi; base += 30) {
            for (std::uint64_t off : W30) {
                const std::uint64_t d = base + off;
                if (d < start || d < 7) continue;
                if (d > limit) { counter += tested; return; }
                ++tested;
                if (divisible_by_u64(N, d)) {
                    found.submit(N, big_from_u64(d), "trial-division");
                    counter += tested;
                    return;
                }
            }
            if ((tested & 0xFFFFF) == 0 && found.ready()) { counter += tested; return; }
        }
    }
    counter += tested;
}

// ---------------------------------------------------------------------
// Lehman multiplier stream -- the completeness guarantee.
//
// If N has no prime factor below N^(1/3), then for some k in
// [1, N^(1/3)] there are a, b with a^2 - b^2 = 4kN and
//
//     ceil(sqrt(4kN)) <= a <= ceil(sqrt(4kN)) + N^(1/6) / (4 sqrt(k)).
//
// Each k is just another wheel-sieved Fermat scan, this time on 4kN.
// Geometrically: multiplying by k re-scales the complex space so that a
// different part of the factor strip is brought down to the near-vertical
// region where Fermat's method is efficient.
// ---------------------------------------------------------------------

// A small-window sieved scan.  The Lehman windows are short
// (N^(1/6) / (4 sqrt k)), so a residue-list wheel would cost far more to
// build and to skip past than the window itself.  Here the same
// quadratic-residue condition is applied incrementally, one modulus at a
// time, directly over the window.
struct SmallSieve {
    static constexpr std::uint32_t MODS[9] = {64, 27, 25, 7, 11, 13, 17, 19, 23};
    std::vector<std::uint8_t> ok[9];
    std::uint64_t start_mod[9];

    // k-independent, so built once per thread rather than once per k.
    std::vector<std::uint8_t> sq[9];
    std::vector<std::uint32_t> sqr[9];   // r^2 mod m, also k-independent

    SmallSieve() {
        for (int i = 0; i < 9; ++i) {
            const std::uint32_t m = MODS[i];
            sq[i] = squares_mod(m);
            sqr[i].resize(m);
            for (std::uint32_t r = 0; r < m; ++r) sqr[i][r] = (r * r) % m;
            ok[i].assign(m, 0);
        }
    }

    // Per k: only the M-dependent shift changes. This used to rebuild
    // squares_mod for all nine moduli on every k -- allocating nine vectors
    // and redoing O(m^2) work per multiplier, with k running to N^(1/3).
    void rebuild(const BigInt& M, const BigInt& lo) {
        for (int i = 0; i < 9; ++i) {
            const std::uint32_t m = MODS[i];
            const std::uint32_t Mm = static_cast<std::uint32_t>(umod(M, m));
            for (std::uint32_t r = 0; r < m; ++r)
                ok[i][r] = sq[i][(sqr[i][r] + m - Mm) % m];
            start_mod[i] = umod(lo, m);
        }
    }
};
constexpr std::uint32_t SmallSieve::MODS[9];

static void lehman_stream(const BigInt& N,
                          std::uint64_t k_max,
                          Found& found,
                          std::atomic<std::uint64_t>& next_k,
                          std::atomic<std::uint64_t>& counter) {
    // Integer 6th root, rounded up: the Lehman window width scale.
    BigInt n16;
    mpz_root(n16.get_mpz_t(), N.get_mpz_t(), 6);
    ++n16;

    SmallSieve ss;
    BigInt x, t, y, g, d;
    std::uint64_t local = 0;

    while (!found.ready()) {
        const std::uint64_t k = next_k.fetch_add(1, std::memory_order_relaxed);
        if (k < 1 || k > k_max) break;

        const BigInt M  = 4 * big_from_u64(k) * N;
        const BigInt lo = isqrt_ceil(M);

        // span = n16 / (4 sqrt k), rounded up.  Using floor(sqrt(k)) makes
        // the divisor no larger than the true one, so the window is a
        // superset of Lehman's -- never a miss.
        std::uint64_t sk = static_cast<std::uint64_t>(std::sqrt(static_cast<double>(k)));
        if (sk < 1) sk = 1;
        BigInt span = n16 / (4 * big_from_u64(sk)) + 2;
        const std::uint64_t span_u = u64_from_big(span, std::numeric_limits<std::uint64_t>::max());

        ss.rebuild(M, lo);

        std::uint64_t cur[9];
        for (int i = 0; i < 9; ++i) cur[i] = ss.start_mod[i];

        for (std::uint64_t off = 0; off <= span_u; ++off) {
            bool pass = true;
            for (int i = 0; i < 9; ++i) {
                if (!ss.ok[i][cur[i]]) { pass = false; break; }
            }

            if (pass) {
                ++local;
                x = lo + big_from_u64(off);
                t = x * x - M;
                if (t >= 0 && is_perfect_square(t, &y)) {
                    d = x > y ? x - y : y - x;
                    mpz_gcd(g.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());
                    if (g > 1 && g < N) { found.submit(N, g, "lehman"); counter += local; return; }
                    d = x + y;
                    mpz_gcd(g.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());
                    if (g > 1 && g < N) { found.submit(N, g, "lehman"); counter += local; return; }
                }
            }

            for (int i = 0; i < 9; ++i) {
                if (++cur[i] == SmallSieve::MODS[i]) cur[i] = 0;
            }
            if ((off & 0xFFFF) == 0 && found.ready()) break;
        }
    }
    counter += local;
}

// ---------------------------------------------------------------------
// Yellow-path stream  (Complex Space Factorization, Sec. 4, Figure 8)
//
// The "key complex number" is (r + (r-3)i) with r = ceil(sqrt(N)), and the
// traversal runs down the 135-degree line through it.  Writing m = r - 3
// and S = r + m, the path is
//
//     R_j = r + j     i_j = m - j     P_j = R - i = 2j + 3     R + i = S
//
// so one index j drives all three of the paper's methods at once:
//
//     vertical Fermat    is V_j = R_j^2 - N a perfect square?
//     horizontal Fermat  is H_j = i_j^2 + N a perfect square?
//     trial division     does P_j divide N?
//
// Two facts make the step cheap (both verified in --selftest):
//
//   1. V and H advance by ADDITION alone:
//          V_{j+1} = V_j + 2 R_j + 1        H_{j+1} = H_j - 2 i_j + 1
//
//   2. V and H are linearly linked:
//          H_j = V_j - P_j * S + 2N
//      so a single state variable carries both tests.
//
// The path stops when P_j * S > N, i.e. at j_max = (N/S - 3)/2, which is
// asymptotically sqrt(N)/4 (measured: 0.2500 * sqrt(N)).  That bound is
// exactly complete: a factor p is caught by the trial-division leg when
// p <= 2 j_max + 3 ~ sqrt(N)/2, and otherwise p > sqrt(N)/2 forces
// j_true = (sqrt(q) - sqrt(p))^2 / 2 <= sqrt(N)/4, so the vertical leg
// catches it.  The two legs cover each other precisely -- which is what
// Figure 8 shows geometrically.
//
// Cost is sqrt(N)/4 iterations of a very tight loop: a few additions, two
// wheel lookups, one small division.  No multiplication, no square root in
// the common case.  Note what fusing costs, though: because a = r + j is
// locked to p = 2j + 3, the wheel can suppress the isqrt but cannot skip
// the iteration, so this stream cannot exploit the ~420,000x the decoupled
// vertical scan gets.  It is a bounded, tight-constant fallback, not a
// replacement for the race.
// ---------------------------------------------------------------------

static void yellow_path_stream(const BigInt& N,
                               Found& found,
                               std::atomic<std::uint64_t>& next_chunk,
                               std::atomic<std::uint64_t>& counter) {
    const BigInt r = isqrt_ceil(N);
    if (r < 4) return;
    const BigInt m = r - 3;
    if (m < 1) return;
    const BigInt S = r + m;

    BigInt jmax_big = (N / S - 3) / 2;
    if (jmax_big < 0) return;
    const std::uint64_t jmax = u64_from_big(jmax_big, std::numeric_limits<std::uint64_t>::max());

    // Admissible-residue tables: R for V = R^2 - N, i for H = i^2 + N.
    static const std::uint32_t MODS[7] = {64, 27, 25, 7, 11, 13, 17};
    std::vector<std::uint8_t> okR[7], okI[7];
    for (int t = 0; t < 7; ++t) {
        const std::uint32_t mm = MODS[t];
        const std::vector<std::uint8_t> sq = squares_mod(mm);
        const std::uint64_t Np = umod(N, mm);
        const std::uint64_t Nn = umod(-N, mm);
        okR[t].assign(mm, 0); okI[t].assign(mm, 0);
        for (std::uint32_t x = 0; x < mm; ++x) {
            const std::uint64_t v1 = ((1ull*x*x) % mm + mm - Np) % mm;   // x^2 - N
            const std::uint64_t v2 = ((1ull*x*x) % mm + mm - Nn) % mm;   // x^2 + N
            okR[t][x] = sq[v1] ? 1 : 0;
            okI[t][x] = sq[v2] ? 1 : 0;
        }
    }

    // Chunk size sets how many independent bracket starts the path is cut
    // into. Nothing here is a ladder: R_j, i_j, V_j, H_j and P_j are all
    // closed forms in j, so a bracket can begin at ANY j and every j is
    // visited exactly once whatever the cut. The only variable cost is the
    // per-chunk start -- two squares plus fourteen modulos -- and measured
    // against the step it drives (mpz, as this stream actually runs):
    //
    //      bits      start     step    ratio
    //        40    211.3ns   87.1ns     2.4
    //        62    211.6ns   87.2ns     2.4
    //        96    211.6ns   87.9ns     2.4
    //       112    213.5ns   87.0ns     2.5
    //
    // The ratio is flat in N -- both ends grow together -- so a chunk of
    // ~120 already holds the start under 2%, at every size.
    //
    // 1<<16 was a CPU tuning choice and caps the path at j_max/65536 starts:
    // only ~3,800 at 60 bits, nowhere near enough to fill a device. 1<<8 is
    // free on CPU -- min-of-11 wall times are 4/4, 54/54 and 126/127 ms on
    // the three cases below -- and gives ~1e6 starts at 60 bits, ~1e9 at 96.
    // The per-chunk atomic costs 16 ns fully contended against ~22,000 ns of
    // chunk work, so it stays far under 1%.
    const std::uint64_t CHUNK = 1u << 8;

    // P_j = R_j - i_j = 2j + 3 is a MACHINE WORD, not a bignum: it runs to
    // 2*j_max+3 = N/S ~ sqrt(N)/2, so it fits 64 bits for any N below 2^128.
    // Forming it as an mpz (P = R - I) and calling mpz_divisible_p costs a
    // subtraction, two increments and a general division every step; passing
    // the machine word to divisible_by_u64 reaches mpz_divisible_ui_p, which
    // divides by a PRECOMPUTED reciprocal because the divisor is one limb.
    // Measured per step, trial-division leg only:
    //
    //                                        96-bit N   112-bit N
    //      P = R - I (mpz), mpz_divisible_p    28.77ns     44.01ns
    //      P = 2j+3 (u64),  divisible_by_u64    5.90ns     10.03ns
    //      ... and skipping 3 | P or 5 | P      3.44ns      3.52ns
    //
    // The skip is sound only when N itself is coprime to 3 and 5: then no P
    // sharing those factors can divide N. When it is not, those P must still
    // be tested, so the skip is switched off rather than assumed.
    const bool p_fits64 =
        jmax <= (std::numeric_limits<std::uint64_t>::max() - 3) / 2;
    const bool skip35 = !divisible_by_u64(N, 3) && !divisible_by_u64(N, 5);
    const bool m_fits64 = mpz_sizeinbase(m.get_mpz_t(), 2) <= 64;
    const std::uint64_t m_w =
        m_fits64 ? u64_from_big(m, std::numeric_limits<std::uint64_t>::max())
                 : std::numeric_limits<std::uint64_t>::max();

    BigInt R, I, V, H, y, g, P;
    std::uint64_t local = 0;

    while (!found.ready()) {
        const std::uint64_t c = next_chunk.fetch_add(1, std::memory_order_relaxed);
        const std::uint64_t j0 = c * CHUNK;
        if (j0 > jmax) break;
        std::uint64_t j1 = std::min(jmax, j0 + CHUNK - 1);
        // i_j = m - j must stay non-negative.  The incremental walk used to
        // notice that via `if (I < 0) break`; with I gone from the hot loop
        // the bound is applied to j up front instead.
        if (m_fits64 && j1 > m_w) j1 = m_w;

        R = r + big_from_u64(j0);
        I = m - big_from_u64(j0);
        if (I < 0) break;

        std::uint64_t rm[7], im[7];
        for (int t = 0; t < 7; ++t) { rm[t] = umod(R, MODS[t]); im[t] = umod(I, MODS[t]); }

        // The bracket start hands the trial-division leg its seed for free:
        // P is just 2*j0+3, and it advances by 2 alongside the traversal.
        std::uint64_t Pw = p_fits64 ? 2 * j0 + 3 : 0;

        for (std::uint64_t j = j0; j <= j1; ++j) {
            ++local;

            // trial-division leg
            if (p_fits64) {
                if (Pw > 1 && !(skip35 && (Pw % 3 == 0 || Pw % 5 == 0)) &&
                    divisible_by_u64(N, Pw)) {
                    found.submit(N, big_from_u64(Pw), "yellow-path/trial-division");
                    counter += local; return;
                }
            } else {
                P = R - I;
                if (P > 1 && mpz_divisible_p(N.get_mpz_t(), P.get_mpz_t())) {
                    found.submit(N, P, "yellow-path/trial-division");
                    counter += local; return;
                }
            }

            // Vertical Fermat leg.  The sieve decides this on machine words
            // alone, so V is formed only for a survivor -- about 1 in 300.
            // Carrying V incrementally instead costs a bignum add, and R and
            // I two increments, on EVERY step to serve that one; materialising
            // R = r + j and V = R^2 - N on demand is one add and one multiply
            // in the rare case and nothing in the common one.  It is also what
            // makes this loop the same shape as csf_yp_kernel, which never
            // forms V or H either.
            bool ok = true;
            for (int t = 0; t < 7; ++t) if (!okR[t][rm[t]]) { ok = false; break; }
            if (ok) {
                R = r + big_from_u64(j);
                V = R * R - N;
                if (V >= 0 && is_perfect_square(V, &y)) {
                    BigInt d = R > y ? R - y : y - R;
                    mpz_gcd(g.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());
                    if (g > 1 && g < N) {
                        found.submit(N, g, "yellow-path/vertical-fermat", &R, &y);
                        counter += local; return;
                    }
                }
            }

            // horizontal Fermat leg, same treatment
            ok = true;
            for (int t = 0; t < 7; ++t) if (!okI[t][im[t]]) { ok = false; break; }
            if (ok) {
                I = m - big_from_u64(j);
                H = I * I + N;
                if (is_perfect_square(H, &y)) {
                    BigInt d = y > I ? y - I : I - y;
                    mpz_gcd(g.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());
                    if (g > 1 && g < N) {
                        found.submit(N, g, "yellow-path/horizontal-fermat", &y, &I);
                        counter += local; return;
                    }
                }
            }

            // advance: machine words only, no mpz touched
            Pw += 2;
            for (int t = 0; t < 7; ++t) {
                if (++rm[t] == MODS[t]) rm[t] = 0;
                im[t] = (im[t] == 0) ? MODS[t] - 1 : im[t] - 1;
            }
            if ((local & 0xFFFF) == 0 && found.ready()) break;
        }
    }
    counter += local;
}

// ---------------------------------------------------------------------
// Complex Trial Multiplication  (Prime Factorization/Complex Trial
// Multiplication.md, and the reference C++ implementation)
//
// The reference runs TWO walks from ONE start point, simultaneously:
//
//     ascending :  if p*q > n  ->  i += 10   else  R += 10
//     descending:  if p*q < n  ->  i -= 10   else  R -= 10
//
// The moves are on R and i, by 10 -- not on p and q independently. Both
// branches of the ascending walk raise q = R + i by 10, and both branches
// of the descending walk lower it by 10, so q is monotone in each: up in
// one, down in the other. That is the ascend and the descend.
//
// START: the first RED point of the 135-degree traversal,
//
//     (R, i) = (r + j_max + 1,  m - j_max - 1)
//
// which for N = 8051 is 112 + 65i -- p = 47, q = 177, p*q = 8319 > N, the
// first point past the crossing. Both walks begin there, and because
// R + i = r + m = S on the whole traversal, both begin at q = S.
//
// The Fermat bracket and trial division start at the OTHER end of the same
// traversal, (r, m) = 90 + 87i for N = 8051, and walk toward the crossing.
// So the two methods are indexed to the same line from opposite ends.
//
// SIEVE: the (R, i) digit classes are the paired Fermat sieve for this N's
// classification -- 4k+1 or 4k-1 crossed with the last digit -- and both
// the increases and the decreases preserve them, since every move is +-10
// on R or on i. A stream synced into a class at the start stays in it.
//
// Worked, N = 8051, class (0,7), start 112+65i synced to (120, 67):
//
//     120+67i  p=53 q=187  9911 > N   ->  R -= 10
//     110+67i  p=43 q=177  7611 < N   ->  i -= 10
//     ...
//      90+ 7i  p=83 q= 97  8051 = N   found, 9 steps
//
// An earlier version of this routine had only the ascending walk. From the
// crossing, ascending can never reach a factor with q below S -- 8051's
// q = 97 against S = 177 -- so it could not factor 8051 at all.
// ---------------------------------------------------------------------

static void ctm_stream(const BigInt& N,
                       Found& found,
                       std::atomic<std::uint64_t>& next_chunk,
                       std::atomic<std::uint64_t>& counter) {
    const BigInt r = isqrt_ceil(N);
    if (r < 4) return;
    const BigInt m = r - 3;
    if (m < 1) return;
    const BigInt S = r + m;
    if (S <= 0) return;
    const BigInt jm = (N / S - 3) / 2;
    if (jm < 0) return;

    // The first red point. Both walks start here, at q = S.
    const BigInt R0 = r + jm + 1;
    const BigInt I0 = m - jm - 1;
    if (I0 < 0) return;

    const auto pairs = fermat_digit_pairs(N);      // paired (R, i) classes
    if (pairs.empty()) return;

    // Ascending runs q up to N/3 (p >= 3); descending runs q down to
    // ceil(sqrt N), where p and q meet.
    const BigInt q_top = N / 3;
    const BigInt q_bot = r;

    const std::uint64_t STEP = 10;
    const std::uint64_t CH   = 1u << 13;
    const BigInt STEPB  = big_from_u64(STEP);
    const BigInt CHSPAN = big_from_u64(CH * STEP);

    const std::uint64_t n_up =
        (q_top > S) ? u64_from_big((q_top - S) / CHSPAN + 1,
                                   std::numeric_limits<std::uint64_t>::max()) : 0;
    const std::uint64_t n_dn =
        (S > q_bot) ? u64_from_big((S - q_bot) / CHSPAN + 1,
                                   std::numeric_limits<std::uint64_t>::max()) : 0;
    const std::uint64_t ndraw = (n_up > n_dn ? n_up : n_dn) * 2 + 2;

    BigInt R, I, P, Q, prod, qa, qb, ps;
    std::uint64_t local = 0;

    while (!found.ready()) {
        const std::uint64_t draw = next_chunk.fetch_add(1, std::memory_order_relaxed);
        if (draw >= ndraw) break;

        const bool up = (draw & 1u) == 0;          // alternate the two walks
        const std::uint64_t c = draw >> 1;
        if (up ? (c >= n_up) : (c >= n_dn)) continue;

        // Chunk bounds in q, measured out from S in the walk's direction.
        if (up) {
            qa = S + big_from_u64(c) * CHSPAN;                 // from
            qb = qa + CHSPAN; if (qb > q_top) qb = q_top;      // to (higher)
            if (qa >= q_top) continue;
        } else {
            qa = S - big_from_u64(c) * CHSPAN;                 // from
            qb = qa - CHSPAN; if (qb < q_bot) qb = q_bot;      // to (lower)
            if (qa <= q_bot) continue;
        }

        for (const auto& pr : pairs) {
            if (found.ready()) break;

            // Seed on the hyperbola at this chunk's q, then sync the pair
            // into its digit class. Both +-10 moves preserve the class.
            ps = N / qa;
            R = (qa + ps) / 2;
            I = (qa - ps) / 2;
            if (I < 0) I = 0;
            R += mod_pos_small(BigInt(pr.a10) - R, 10);
            I += mod_pos_small(BigInt(pr.b10) - I, 10);

            const bool native = mpz_sizeinbase(N.get_mpz_t(), 2) <= 126 &&
                                mpz_sizeinbase(qb.get_mpz_t(), 2) <= 62 &&
                                mpz_sizeinbase(qa.get_mpz_t(), 2) <= 62;

            if (native) {
                const __uint128_t Nn = u128_from_big(N);
                std::int64_t Rn = (std::int64_t)u64_from_big(R, 0);
                std::int64_t In = (std::int64_t)u64_from_big(I, 0);
                const std::int64_t lim = (std::int64_t)u64_from_big(qb, ~0ull);
                while (In >= 0 && Rn > In) {
                    const std::int64_t q = Rn + In;
                    if (up ? (q > lim) : (q < lim)) break;
                    // No p < 3 break. The reference has none, and it would
                    // be wrong: p*q < N when p is tiny, so the rule moves i
                    // and RAISES p by 10 -- the walk recovers. Breaking
                    // discards the class first. N = 143 fails that way: its
                    // seed passes through p = 1 one step before p = 11.
                    const std::int64_t p = Rn - In;
                    ++local;
                    const __uint128_t z = (__uint128_t)(std::uint64_t)p *
                                          (std::uint64_t)q;
                    if (z == Nn) {
                        found.submit(N, big_from_u64((std::uint64_t)p), "ctm");
                        counter += local; return;
                    }
                    if (up) { if (z > Nn) In += STEP; else Rn += STEP; }
                    else    { if (z < Nn) In -= STEP; else Rn -= STEP; }
                    if ((local & 0xFFFF) == 0 && found.ready()) break;
                }
            } else {
                while (I >= 0 && R > I) {
                    Q = R + I;
                    if (up ? (Q > qb) : (Q < qb)) break;
                    P = R - I;              // no P < 3 break; see above
                    ++local;
                    prod = P * Q;
                    if (prod == N) {
                        found.submit(N, P, "ctm");
                        counter += local; return;
                    }
                    if (up) { if (prod > N) I += STEPB; else R += STEPB; }
                    else    { if (prod < N) I -= STEPB; else R -= STEPB; }
                    if ((local & 0xFFFF) == 0 && found.ready()) break;
                }
            }
        }
    }
    counter += local;
}

// ---------------------------------------------------------------------
// Iterated-average GCD stream  (Complex Space Factorization, Sec. 5)
// with the remainder sieve of Sec. 5.2.
//
//   a0    = ceil(sqrt(N))
//   IA_k  = (IA_{k-1} + N) / 2,  IA_0 = a0
//         = (N*(2^k - 1) + a0) / 2^k  =  Num_k / 2^k
//
// The paper's test is GCD(IA_k - b/2^k, N).  Writing T = (Num_k - b)/2^k,
// T is an integer exactly when b = Num_k (mod 2^k).  Intersecting that
// with the terminal-digit condition b = db (mod 10) gives, by CRT,
//
//        b = r (mod 5 * 2^k),
//
// and since b advances by 5*2^k while T advances by 5, the scan visits
// every 5th integer below IA_k -- which is precisely the Sec. 5.2 claim,
// and the reason the admissible b per series is 2^(k-1).
//
// IMPORTANT -- the scope of the underlying identity.
//
// T_k is only defined when 2^k divides (N - a0 + b) exactly.  N is odd,
// so gcd(2^k, N) = 1 and that division cannot change a gcd with N:
//
//     gcd(T_k, N) = gcd(N - a0 + b, N) = gcd(a0 - b, N)   for EVERY k
//
// Every level of the iterated average returns the identical gcd.  The k
// iterations are not k independent chances at a factor -- they are one
// arithmetic fact re-tested k times, stopping at k > v2(q-1) where T_k
// ceases to be an integer.  (N = 8051: q-1 = 96 = 2^5*3, so k = 1..5.)
//
// Substituting a0 = a - j, where j = a - ceil(sqrt(N)):
//
//     a0 - b = p - j        a0 + b = q - j
//
// so the descending branch asks whether p - j divides N and the ascending
// branch asks the same of q - j.  Equivalently, via N - p = p(q-1),
//
//     p | T_k  <=>  p | j          q | T_k  <=>  q | (p - j)
//
// With 0 <= j < p -- true unless q/p >= 3 + 2*sqrt(2) = 5.828 -- both
// collapse to j = 0.  Measured over 300 random semiprimes: 51/51 when
// j = 0, 0/249 when j > 0.
//
// And j = 0 means b <~ sqrt(2) * N^(1/4), which is exactly the condition
// for the vertical scan to succeed on its FIRST trial: a = ceil(sqrt(N))
// is candidate number one.  So the identity fires only where one isqrt
// would already have finished.
//
// The Sec. 5.2 sieve below is correct and is implemented, but what it
// accelerates is a blind GCD scan.  This stream is therefore OFF by
// default and is NOT exhaustive: when it is selected alone and fails, the
// caller reports INCOMPLETE rather than presenting N as its own factor.
// ---------------------------------------------------------------------

static std::optional<BigInt> crt_b_residue(int db, const BigInt& r2, int k) {
    // b = db (mod 10), b = r2 (mod 2^k)  ->  b mod 5*2^k
    if (k <= 0) return BigInt(db);
    BigInt pow2 = BigInt(1) << k;
    if (umod(r2, 2) != static_cast<std::uint64_t>(db % 2)) return std::nullopt;

    // b = r2 + 2^k * s,  need r2 + 2^k s = db (mod 5)
    const std::uint64_t r5  = umod(r2, 5);
    const std::uint64_t p5  = umod(pow2, 5);           // 2^k mod 5, never 0
    const std::uint64_t tgt = (static_cast<std::uint64_t>(db % 5) + 5 - r5) % 5;

    std::uint64_t inv = 0;
    for (std::uint64_t i = 1; i < 5; ++i) if ((p5 * i) % 5 == 1) { inv = i; break; }
    const std::uint64_t s = (tgt * inv) % 5;

    BigInt b = r2 + pow2 * big_from_u64(s);
    BigInt mod = 5 * pow2;
    b %= mod; if (b < 0) b += mod;
    return b;
}

static void iterated_average_stream(const BigInt& N,
                                    std::uint64_t steps_per_level,
                                    Found& found,
                                    std::atomic<std::uint64_t>& next_job,
                                    std::atomic<std::uint64_t>& counter) {
    const BigInt a0 = isqrt_ceil(N);
    const BigInt b_max = (N - 9) / 6;
    if (b_max < 1) return;

    const auto pairs = fermat_digit_pairs(N);
    std::vector<int> b_digits;
    for (const auto& pr : pairs)
        if (std::find(b_digits.begin(), b_digits.end(), pr.b10) == b_digits.end())
            b_digits.push_back(pr.b10);
    if (b_digits.empty()) return;

    const int k_max = std::max(1, static_cast<int>(mpz_sizeinbase(b_max.get_mpz_t(), 2)));
    const std::uint64_t jobs = static_cast<std::uint64_t>(k_max) * b_digits.size() * 2;

    BigInt g, T;
    std::uint64_t gcds = 0;

    while (!found.ready()) {
        const std::uint64_t j = next_job.fetch_add(1, std::memory_order_relaxed);
        if (j >= jobs) break;

        const int k        = static_cast<int>(j / (b_digits.size() * 2)) + 1;
        const std::size_t t = static_cast<std::size_t>(j % (b_digits.size() * 2));
        const int db       = b_digits[t / 2];
        const bool descend = (t % 2) == 0;

        const BigInt pow2  = BigInt(1) << k;
        const BigInt Num   = N * (pow2 - 1) + a0;          // IA_k = Num / 2^k

        // descending: T = (Num - b)/2^k  -> b = Num (mod 2^k)
        // ascending : T = (Num + b)/2^k  -> b = -Num (mod 2^k)
        BigInt r2 = Num % pow2;
        if (!descend) r2 = (pow2 - r2) % pow2;
        if (r2 < 0) r2 += pow2;

        auto b0 = crt_b_residue(db, r2, k);
        if (!b0) continue;

        const BigInt stride = 5 * pow2;
        BigInt b = *b0;
        if (b == 0) b = stride;

        for (std::uint64_t s = 0; s < steps_per_level && b <= b_max; ++s) {
            if (descend) T = (Num - b) / pow2;
            else         T = (Num + b) / pow2;
            if (T <= 1) break;
            ++gcds;
            mpz_gcd(g.get_mpz_t(), T.get_mpz_t(), N.get_mpz_t());
            if (g > 1 && g < N) {
                found.submit(N, g, "iterated-average");
                counter += gcds;
                return;
            }
            b += stride;
            if ((s & 0x3FF) == 0 && found.ready()) break;
        }
    }
    counter += gcds;
}

// ---------------------------------------------------------------------
// Options and orchestration
// ---------------------------------------------------------------------

struct Options {
    unsigned threads   = 0;
    std::uint64_t b1   = 1000000;
    bool lehman        = true;
    bool hf            = false;
    bool ia            = false;
    bool yp            = true;    // on by default; --no-yp to ablate
    bool ctm           = false;
    bool parallel      = false;   // yellow path across the whole pool
    bool digit_only    = false;
    bool no_sieve      = false;
    bool quiet         = false;
    bool verbose       = false;
    std::string only;          // "", "vf", "hf", "ia", "td", "lm"
    std::uint64_t wheel = 0;   // 0 = pick from bit-length
};

// Residue budget for the sieve wheel.
//
// Scaling this with the bit-length of N was wrong: what the wheel has to
// pay for is the SCAN DEPTH, and the sieve's returns saturate long before
// the table does. Measured (4 threads):
//
//   budget      shallow 112-bit (j~1e7)   deep 96-bit (j~2.9e11)
//    4,096            6 ms                      662 ms
//   65,536           11 ms                      372 ms
//  262,144           11 ms                      191 ms   <- saturated
// 1,048,576          52 ms                      190 ms
// 2,097,152          50 ms                      196 ms   <- was the default
//
// Past 2^18 the extra residues buy nothing on a deep scan and cost 8x on a
// shallow one, because building sec_res is O(budget * moduli) modulos up
// front. Capping there is a strict improvement at both ends. Override with
// --wheel=R; batch work over many shallow N wants a smaller value still.
static std::uint64_t wheel_budget_for(const BigInt& N) {
    const std::size_t bits = mpz_sizeinbase(N.get_mpz_t(), 2);
    if (bits < 48)  return 1u << 12;
    if (bits < 72)  return 1u << 16;
    return 1u << 18;
}

// Integer cube root, rounded up.
static BigInt icbrt_ceil(const BigInt& n) {
    BigInt r;
    mpz_root(r.get_mpz_t(), n.get_mpz_t(), 3);
    if (r * r * r < n) ++r;
    return r;
}

// Find one nontrivial split of a composite N, or nullopt if N is prime.
static std::optional<BigInt> split_once(const BigInt& N,
                                        const Options& opt,
                                        Stats& stats,
                                        std::string& how) {
    if (N <= 3) return std::nullopt;

    if (N % 2 == 0) { how = "even"; return BigInt(2); }

    // Perfect power.
    {
        BigInt root;
        if (mpz_perfect_power_p(N.get_mpz_t()) != 0) {
            for (unsigned e = 2; e <= mpz_sizeinbase(N.get_mpz_t(), 2); ++e) {
                if (mpz_root(root.get_mpz_t(), N.get_mpz_t(), e) != 0) {
                    how = "perfect-power";
                    return root;
                }
            }
        }
    }

    // Stage 0: trial division to b1.
    {
        std::uint64_t d_found = 0;
        static const std::uint64_t small[3] = {3, 5, 7};
        for (std::uint64_t d : small)
            if (divisible_by_u64(N, d)) { d_found = d; break; }
        if (!d_found) {
            static const std::uint64_t W30[8] = {1, 7, 11, 13, 17, 19, 23, 29};
            for (std::uint64_t base = 0; base <= opt.b1 && !d_found; base += 30)
                for (std::uint64_t off : W30) {
                    const std::uint64_t d = base + off;
                    if (d < 11 || d > opt.b1) continue;
                    if (divisible_by_u64(N, d)) { d_found = d; break; }
                }
        }
        if (d_found) { how = "trial-division(stage0)"; return big_from_u64(d_found); }
    }

    if (is_probable_prime(N)) return std::nullopt;

    // Stage 1: the race.
    Found found;
    unsigned T = opt.threads ? opt.threads : std::thread::hardware_concurrency();
    if (T == 0) T = 1;

    const BigInt sqrtN = isqrt_floor(N);
    const BigInt cbrtN = icbrt_ceil(N);

    // -----------------------------------------------------------------
    // Partition the strip at the 135-degree traversal's own boundary, so
    // that no two streams ever cover the same p.
    //
    // The Sec. 4 bracket runs down the line R + i = S through the key
    // complex number (r + m*i), with r = ceil(sqrt N), m = r - 3,
    // S = r + m, halting at j_max = (N/S - 3)/2. Across that traversal the
    // two Fermat legs sweep DISJOINT arcs of the strip that abut exactly at
    // the collapse: the measured gap between them is +0.28 to +0.60 in R,
    // never negative, from N = 309 to 1e12. Its trial-division leg reaches
    // p = 2*j_max + 3.
    //
    // Take that as the crossover, c = 2*j_max + 3:
    //
    //     trial division owns  p in (b1, c]
    //     vertical Fermat owns a in [ceil(sqrt N), (c + N/c)/2]
    //
    // Disjoint, and complete because any p > c has a <= (c + N/c)/2. With
    // c ~ sqrt(N)/2 the vertical arc is a in [sqrt N, 1.25 sqrt N].
    //
    // The old bound used b1 in place of c, so the vertical scan ran to
    // a = (b1 + N/b1)/2 -- the entire strip, exactly the p range trial
    // division was already sweeping. The two duplicated each other and the
    // vertical scan could never terminate. Bounded properly it exhausts its
    // arc in finite time and hands its threads to trial division.
    //
    // Note the legs need NOT advance in lockstep on j to stay disjoint --
    // only to stay inside their own arc. That distinction matters: lockstep
    // would forfeit the wheel, which skips ~934 of every 935 values of a
    // and can only do so if the scan advances at its own rate.
    // -----------------------------------------------------------------
    BigInt cross;
    {
        const BigInt rr = isqrt_ceil(N);
        const BigInt mm = rr - 3;
        if (mm >= 1) {
            const BigInt SS = rr + mm;
            cross = 2 * ((N / SS - 3) / 2) + 3;
        } else {
            cross = sqrtN;
        }
        const BigInt b1b = big_from_u64(opt.b1 < 3 ? 3 : opt.b1);
        if (cross < b1b)  cross = b1b;
        if (cross > sqrtN) cross = sqrtN;
        if (cross < 3)    cross = 3;
    }
    BigInt vf_hi = (cross + N / cross) / 2 + 1;
    const BigInt vf_lo = isqrt_ceil(N);
    if (vf_hi < vf_lo) vf_hi = vf_lo;

    const std::uint64_t budget = opt.wheel ? opt.wheel : wheel_budget_for(N);
    Wheel wv = build_wheel(N, budget, opt.digit_only, opt.no_sieve);
    Wheel wh;
    if (opt.hf || opt.only == "hf") wh = build_wheel(-N, budget, opt.digit_only, opt.no_sieve);

    std::uint64_t k_max = 0;
    if (opt.lehman || opt.only == "lm") {
        k_max = u64_from_big(cbrtN, std::numeric_limits<std::uint64_t>::max());
    }

    std::atomic<std::uint64_t> vf_chunk{0}, hf_chunk{0}, lm_k{1}, ia_job{0}, yp_chunk{0};
    std::atomic<std::uint64_t> td_block{0}, ctm_chunk{0};

    // Thread budget: one for TD, optionally one each for HF and IA, the
    // rest split between vertical Fermat and the Lehman sweep.
    // --only isolates a single stream (used for the method comparison in
    // the write-up; the default is the full race).
    const bool only_mode = !opt.only.empty();
    const bool want_td = !only_mode || opt.only == "td";
    const bool want_vf = !only_mode || opt.only == "vf";
    (void)0;
    const bool want_hf = only_mode ? (opt.only == "hf") : opt.hf;
    const bool want_ia = only_mode ? (opt.only == "ia") : opt.ia;
    const bool want_lm = only_mode ? (opt.only == "lm") : opt.lehman;
    const bool want_yp = only_mode ? (opt.only == "yp") : opt.yp;
    const bool want_ctm = only_mode ? (opt.only == "ctm") : opt.ctm;

    unsigned reserved = 1 + (opt.hf ? 1u : 0u) + (opt.ia ? 1u : 0u) + (opt.yp ? 1u : 0u)
                      + (opt.ctm ? 1u : 0u);
    unsigned rest = (T > reserved) ? T - reserved : 1;
    // Thread split. Lehman is the O(N^(1/3)) completeness guarantee, not
    // the workhorse -- it never won a race in any measured case. Giving it
    // half the pool starved the vertical scan, which IS the workhorse for
    // balanced N: on a 107-bit case with j = 4.09e13, --no-lehman ran 3.2x
    // faster (132 s -> 41.5 s) purely because vf went from 1 thread to 3.
    //
    // One Lehman thread still sweeps every k, so the asymptotic bound is
    // unchanged; only its constant is. Everything else goes to vf.
    unsigned nLM = (want_lm && rest >= 2) ? 1u : 0u;
    unsigned nVF = rest - nLM;
    if (nVF == 0) { nVF = 1; nLM = 0; }

    // Streams that carve their range with an atomic chunk counter take as
    // many threads as they are given. Until measured, only vf and lm were
    // ever handed the pool: yp and ctm were written chunk-parallel but
    // pinned to one thread each, even under --only, so --only=yp used a
    // single core while --only=vf used all of them.
    unsigned nYP = want_yp ? 1u : 0u;
    unsigned nCTM = want_ctm ? 1u : 0u;
    if (only_mode) {
        nVF  = want_vf  ? T : 0;
        nLM  = want_lm  ? T : 0;
        nYP  = want_yp  ? T : 0;
        nCTM = want_ctm ? T : 0;
    }

    // Trial division owns p in (b1, cross]; the vertical scan owns the rest.
    // Lehman still races as the O(N^(1/3)) worst-case guarantee.
    const std::uint64_t td_limit =
        u64_from_big(cross, std::numeric_limits<std::uint64_t>::max());

    std::vector<std::thread> pool;
    auto run_td = [&] {
        trial_division_stream(N, opt.b1 + 1, td_limit, found, td_block, stats.td_tested);
    };
    if (want_td) pool.emplace_back(run_td);

    for (unsigned i = 0; i < nVF; ++i)
        pool.emplace_back([&] {
            fermat_scan_chunked(N, N, wv, vf_lo, vf_hi, found, vf_chunk,
                                stats.vf_candidates, stats.vf_scanned, "vertical-fermat");
            // Arc exhausted with no hit: the factor must be in trial
            // division's region, so migrate there rather than returning.
            if (!found.ready() && want_td) run_td();
        });

    if (want_hf)
        pool.emplace_back([&] {
            std::atomic<std::uint64_t> dummy{0};
            fermat_scan_chunked(N, -N, wh, BigInt(0), (N - 9) / 6, found, hf_chunk,
                                stats.hf_candidates, dummy, "horizontal-fermat");
        });

    if (want_ia)
        pool.emplace_back([&] {
            iterated_average_stream(N, 1u << 22, found, ia_job, stats.ia_gcds);
        });

    for (unsigned i = 0; i < nYP; ++i)
        pool.emplace_back([&] {
            yellow_path_stream(N, found, yp_chunk, stats.yp_steps);
        });

    // Launched at the collapse SIMULTANEOUSLY with the bracket, never after.
    for (unsigned i = 0; i < nCTM; ++i)
        pool.emplace_back([&] {
            ctm_stream(N, found, ctm_chunk, stats.ctm_steps);
        });

    for (unsigned i = 0; i < nLM; ++i)
        pool.emplace_back([&] {
            lehman_stream(N, k_max, found, lm_k, stats.lm_candidates);
        });

    for (auto& t : pool) t.join();

    if (found.ready() && found.factor > 1 && found.factor < N) {
        how = found.source;
        if (found.certified) how += " [complex-square certified]";
        return found.factor;
    }
    return std::nullopt;
}

// Full factorization.
static void factor_recursive(const BigInt& N,
                             const Options& opt,
                             Stats& stats,
                             std::vector<BigInt>& out,
                             std::vector<std::string>& hows,
                             std::vector<BigInt>& unfactored) {
    if (N <= 1) return;
    if (is_probable_prime(N)) { out.push_back(N); return; }

    std::string how;
    auto f = split_once(N, opt, stats, how);
    if (!f) {
        // N is composite but this configuration's streams were exhausted
        // without finding a split.  That happens only when the engine is
        // run with an incomplete stream selected (--only=ia, for example).
        // Emitting N as if it were a factor would be a silent wrong answer,
        // so record it separately and let the caller report it.
        out.push_back(N);
        unfactored.push_back(N);
        return;
    }

    hows.push_back(how);
    factor_recursive(*f, opt, stats, out, hows, unfactored);
    factor_recursive(N / *f, opt, stats, out, hows, unfactored);
}

// ---------------------------------------------------------------------
// Self-test
// ---------------------------------------------------------------------

static int selftest() {
    int failures = 0;
    std::cout << "=== 1. Corrected terminal-digit sieve vs brute force ===\n";

    // Ground truth by direct enumeration of N = (R-i)(R+i).
    std::vector<std::vector<std::vector<char>>> truth(
        20, std::vector<std::vector<char>>(10, std::vector<char>(10, 0)));
    std::vector<BigInt> rep(20, 0);

    const long LIM = 2000;
    for (long R = 1; R < LIM; ++R)
        for (long i = 0; i < R; ++i) {
            const long p = R - i, q = R + i;
            if (p < 1 || p % 2 == 0 || q % 2 == 0) continue;
            const long Nv = p * q;
            const int n20 = static_cast<int>(Nv % 20);
            truth[n20][R % 10][i % 10] = 1;
            if (rep[n20] == 0) rep[n20] = Nv;
        }

    for (int n20 = 1; n20 < 20; n20 += 2) {
        if (rep[n20] == 0) continue;
        const BigInt Nv = rep[n20];
        auto got = fermat_digit_pairs(Nv);

        std::vector<std::vector<char>> mark(10, std::vector<char>(10, 0));
        for (const auto& pr : got) mark[pr.a10][pr.b10] = 1;

        int missing = 0, extra = 0, tcount = 0;
        for (int a = 0; a < 10; ++a)
            for (int b = 0; b < 10; ++b) {
                if (truth[n20][a][b]) ++tcount;
                if (truth[n20][a][b] && !mark[a][b]) ++missing;
                if (!truth[n20][a][b] && mark[a][b]) ++extra;
            }

        const bool ok = (missing == 0 && extra == 0);
        if (!ok) ++failures;
        std::cout << "  N=" << n20 << " (mod 20)  [N%4=" << (n20 % 4)
                  << ", N%10=" << (n20 % 10) << "]  classes: got " << got.size()
                  << ", truth " << tcount
                  << (ok ? "   OK" : "   FAIL") ;
        if (!ok) std::cout << " (missing " << missing << ", extra " << extra << ")";
        std::cout << "\n";
    }

    std::cout << "\n=== 2. Is the digit sieve exactly the QR sieve at modulus 20? ===\n";
    for (int n20 = 1; n20 < 20; n20 += 2) {
        if (rep[n20] == 0) continue;
        const BigInt Nv = rep[n20];

        // digit-sieve admissible a mod 10
        std::vector<char> A_digit(10, 0);
        for (const auto& pr : fermat_digit_pairs(Nv)) A_digit[pr.a10] = 1;

        // QR-sieve admissible a mod 20, projected to mod 10
        std::vector<char> A_qr(10, 0);
        const auto sq20 = squares_mod(20);
        const std::uint64_t M20 = umod(Nv, 20);
        for (std::uint64_t a = 0; a < 20; ++a) {
            const std::uint64_t v = ((a * a) % 20 + 20 - M20) % 20;
            if (sq20[v]) A_qr[a % 10] = 1;
        }

        bool same = true;
        for (int d = 0; d < 10; ++d) if (A_digit[d] != A_qr[d]) same = false;
        std::cout << "  N=" << n20 << " (mod 20): digit {";
        for (int d = 0; d < 10; ++d) if (A_digit[d]) std::cout << d << " ";
        std::cout << "}  QR20 {";
        for (int d = 0; d < 10; ++d) if (A_qr[d]) std::cout << d << " ";
        std::cout << "}  " << (same ? "identical" : "DIFFER") << "\n";
    }

    std::cout << "\n=== 3. Sieve strength (fraction of a surviving) ===\n";
    {
        const BigInt Nv("104729", 10);
        const BigInt Ntest = BigInt("104729") * BigInt("104723");
        (void)Nv;
        Wheel dw = build_wheel(Ntest, 1u << 21, true,  false);
        Wheel fw = build_wheel(Ntest, 1u << 21, false, false);
        std::cout << "  digit sieve  (mod " << dw.W << "): density "
                  << dw.density << "  -> " << (1.0 / dw.density) << "x\n";
        std::cout << "  full wheel   (mod " << fw.W << ", +" << fw.sec_mod.size()
                  << " secondary): density " << fw.density
                  << "  -> " << (1.0 / fw.density) << "x\n";
    }

    std::cout << "\n=== 4. Wheel soundness: the true a is never sieved out ===\n";
    {
        // For a spread of semiprimes, the genuine a = (p+q)/2 must appear
        // in the wheel's admissible residue list.  A sieve that ever drops
        // it would make the engine silently incomplete.
        std::mt19937_64 rng(12345);
        int checked = 0, dropped = 0;
        for (int trial = 0; trial < 400; ++trial) {
            BigInt p = 3 + BigInt(static_cast<unsigned long>(rng() % 1000000));
            BigInt q = 3 + BigInt(static_cast<unsigned long>(rng() % 1000000));
            mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
            mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
            if (p == q) continue;
            const BigInt Nv = p * q;
            if (umod(Nv, 2) == 0) continue;

            Wheel w = build_wheel(Nv, 1u << 16, false, false);
            if (w.empty()) { ++dropped; continue; }

            const BigInt a = (p + q) / 2;
            const std::uint64_t ra = umod(a, w.W);
            const bool in_primary =
                std::binary_search(w.residues.begin(), w.residues.end(),
                                   static_cast<std::uint32_t>(ra));
            bool in_secondary = true;
            for (std::size_t i = 0; i < w.sec_mod.size(); ++i)
                if (!w.sec_ok[i][umod(a, w.sec_mod[i])]) in_secondary = false;

            ++checked;
            if (!in_primary || !in_secondary) ++dropped;
        }
        if (dropped) ++failures;
        std::cout << "  checked " << checked << " semiprimes, true a dropped by the sieve: "
                  << dropped << (dropped ? "   FAIL" : "   OK") << "\n";
    }

    std::cout << "\n=== 5. Factorization correctness ===\n";
    struct Case { const char* n; const char* label; };
    const std::vector<Case> cases = {
        {"8051",                     "paper example 83*97"},
        {"309",                      "paper example 3*103"},
        {"143",                      "paper example 11*13"},
        {"10963",                    "N=3 mod 4"},
        {"1000003",                  "prime"},
        {"10967535067",              "balanced semiprime"},
        {"1000003000000000000000000000000000000021", "40-digit, small factors"},
        {"10000000000000000000000000000000000000121", "41-digit prime"},
        {"1099511627791",            "prime (2^40 + 15)"},
        {"3825123056546413051",      "hard 62-bit"},
        {"152415787532388367501905199875019052100",  "highly composite"},
        {"9999999999999999999999999999999999999999", "5 | N and 3 | N (repunit-like)"},
        {"59393383530518219041528236287",  "96-bit balanced, vertical depth j ~ 1e6"},
        {"3168987219233877513774136225800517", "112-bit balanced, vertical depth j ~ 1e7"},
        {"1000000016000000063",      "balanced 60-bit (1000000007*1000000009)"},
        {"25",                       "perfect square"},
        {"2147483647",               "Mersenne prime"},
    };

    Options o; o.quiet = true; o.threads = 4; o.b1 = 100000;
    for (const auto& c : cases) {
        BigInt N(c.n, 10);
        Stats st;
        std::vector<BigInt> fs; std::vector<std::string> hows; std::vector<BigInt> unf;
        auto t0 = std::chrono::steady_clock::now();
        factor_recursive(N, o, st, fs, hows, unf);
        auto t1 = std::chrono::steady_clock::now();

        BigInt prod = 1;
        for (const auto& f : fs) prod *= f;
        bool allprime = true;
        for (const auto& f : fs) if (!is_probable_prime(f)) allprime = false;

        const bool ok = (prod == N) && allprime && unf.empty();
        if (!ok) ++failures;
        std::sort(fs.begin(), fs.end());

        std::cout << "  " << (ok ? "OK  " : "FAIL") << "  " << c.n << " = ";
        for (std::size_t i = 0; i < fs.size(); ++i)
            std::cout << (i ? " * " : "") << fs[i];
        std::cout << "   [" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                  << " ms]  (" << c.label << ")\n";
    }

    // -----------------------------------------------------------------
    // 6. How the yellow strip grows, and why it stops where it does.
    //
    //   P_j = 2j+3 and Q_j = r+m = S are the two coordinates; S is CONSTANT
    //   along the 135-degree line, so the path meets the hyperbola P*Q = N
    //   at exactly
    //
    //       j_max = floor((floor(N/S) - 3) / 2),      S = 2r - 3
    //
    //   Expanding r = sqrt(N) + e with e = ceil(sqrt N) - sqrt(N) in [0,1):
    //
    //       j_max = sqrt(N)/4 - 9/8 - e/4 + O(1/sqrt N)
    //
    //   so the strip is sqrt(N)/4 long -- 2^(b/2 - 2) for a b-bit N. For
    //   N = 309 that is j = 0,1,2,3: a strip of 4.
    //
    //   The constant 1/4 is not arbitrary. Writing the smaller factor as
    //   p = sqrt(N)/t, the trial-division leg reaches 2*j_max+3 = N/S ->
    //   sqrt(N)/2 and so covers t >= 2 outright; for 1 <= t < 2 the vertical
    //   leg needs depth sqrt(N)*(t-1)^2/(2t), which increases on [1,2] and
    //   hits sqrt(N)/4 exactly at t = 2. The two legs meet at the worst
    //   case and neither is longer than it has to be.
    // -----------------------------------------------------------------
    std::cout << "\n=== 6. Yellow strip: growth law and coverage ===\n";
    {
        const struct { const char* n; long want; } exact[] = {
            {"309", 3}, {"143", 1}, {"8051", 21}, {"798607", 222},
        };
        for (const auto& c : exact) {
            BigInt n(c.n), rr = isqrt_ceil(n);
            const BigInt S = rr + (rr - 3);
            const BigInt got = (n / S - 3) / 2;
            const bool ok = (got == BigInt(c.want));
            if (!ok) ++failures;
            std::cout << "  N=" << c.n << ": j_max=" << got
                      << " (strip length " << got + 1 << ")  "
                      << (ok ? "OK" : "FAIL") << "\n";
        }

        // The growth law, stated so it can be checked EXACTLY at any size.
        //
        //   j_max = sqrt(N)/4 - 9/8 - e/4 + O(1/sqrt N)
        //
        // in floating point loses all its meaning past ~96 bits: at 128
        // bits sqrt(N)/4 is near 2^62 and a double resolves it only to
        // about 512, while the quantity being bounded is O(1). Multiplying
        // through by 4 turns it into integers, which are exact at any size:
        //
        //   4*j_max = isqrt(N) - c,   c a small bounded constant
        //
        // Sampled over 3006 values of N from 12 to 512 bits, c lands in
        // [4, 9] and never drifts -- the two floors in the exact formula
        // account for the width.
        long c_lo = 1 << 30, c_hi = -(1 << 30);
        std::mt19937_64 grng(271828);
        for (int b = 12; b <= 512; b += 4) {
            BigInt n = 1; mpz_mul_2exp(n.get_mpz_t(), n.get_mpz_t(), b);
            n += big_from_u64(grng() | 1ull);
            BigInt rr = isqrt_ceil(n);
            const BigInt S = rr + (rr - 3);
            if (S <= 0) continue;
            const BigInt jm = (n / S - 3) / 2;
            BigInt fl; mpz_sqrt(fl.get_mpz_t(), n.get_mpz_t());
            const BigInt c = fl - 4 * jm;
            const long cv = mpz_get_si(c.get_mpz_t());
            c_lo = std::min(c_lo, cv); c_hi = std::max(c_hi, cv);
        }
        const bool c_ok = (c_lo >= 0 && c_hi <= 16);
        if (!c_ok) ++failures;
        std::cout << "  4*j_max = isqrt(N) - c over 12..512 bits: c in ["
                  << c_lo << ", " << c_hi << "]  -> strip ~ sqrt(N)/4  "
                  << (c_ok ? "OK" : "FAIL") << "\n";

        // coverage: every semiprime is caught by one leg or the other, and
        // where the trial leg does NOT reach, j_true stays under sqrt(N)/4.
        std::mt19937_64 rng(31415);
        std::size_t miss = 0, beyond = 0; double worst_r = 0.0;
        for (int k = 0; k < 600; ++k) {
            const int e = 6 + (int)(rng() % 18);
            BigInt p, q, seed = big_from_u64(3 + rng() % (1ull << e));
            mpz_nextprime(p.get_mpz_t(), seed.get_mpz_t());
            seed = p + big_from_u64(rng() % (1ull << e));
            mpz_nextprime(q.get_mpz_t(), seed.get_mpz_t());
            const BigInt n = p * q;
            const BigInt rr = isqrt_ceil(n), S = rr + (rr - 3);
            if (S <= 0) continue;
            const BigInt jm = (n / S - 3) / 2;
            const BigInt reach = 2 * jm + 3;
            const BigInt j_true = (p + q) / 2 - rr;
            const bool by_trial = (p <= reach);
            const bool by_vert  = (j_true >= 0 && j_true <= jm);
            if (!by_trial && !by_vert) ++miss;
            if (p > reach) {
                ++beyond;
                const double r2 = mpz_get_d(j_true.get_mpz_t()) /
                                  std::sqrt(mpz_get_d(n.get_mpz_t()));
                worst_r = std::max(worst_r, r2);
            }
        }
        const bool cov_ok = (miss == 0) && (worst_r <= 0.2500);
        if (!cov_ok) ++failures;
        std::cout << "  600 semiprimes: " << miss << " uncovered; of the "
                  << beyond << " beyond the trial leg, worst j_true/sqrt(N) = "
                  << worst_r << " (bound 0.25)  " << (cov_ok ? "OK" : "FAIL") << "\n";

        // Zero-margin cases: j_true == j_max exactly, so the factor sits on
        // the LAST step of the strip and the path would miss it if it were
        // one step shorter. These pin the strip length from below.
        //
        // Finding them takes care on two counts. First, the tight corner is
        // NOT at t = 2: there p = sqrt(N)/2 is exactly the trial leg's reach,
        // so the trial leg still covers and the vertical leg is never asked.
        // The cases that bind sit just PAST the reach, where the vertical leg
        // is the only thing left. Second, they are rare enough that a random
        // sweep over a wide range finds none in 400,000 tries; these came
        // from an exhaustive scan of p in [500, 40000) with q near 4p.
        //
        // All three are stated for r = ceil(sqrt N), which is the convention
        // this engine uses. Under floor(sqrt N) the same strip is indexed
        // from one point lower and these margins are not zero -- a different
        // set of N is tight there.
        const struct { const char* n; const char* p; } tight[] = {
            {"1007509",  "503"},
            {"1289923",  "569"},
            {"17005309", "2063"},
        };
        for (const auto& c : tight) {
            const BigInt n(c.n), pf(c.p), qf = n / pf;
            const BigInt rr = isqrt_ceil(n), S = rr + (rr - 3);
            const BigInt jm = (n / S - 3) / 2;
            const BigInt j_true = (pf + qf) / 2 - rr;
            const bool ok = (j_true == jm);
            if (!ok) ++failures;
            std::cout << "  zero margin: N=" << c.n << " j_true=" << j_true
                      << " j_max=" << jm << " (factor on the last step)  "
                      << (ok ? "OK" : "FAIL") << "\n";
        }
    }

    std::cout << "\n" << (failures == 0 ? "ALL TESTS PASSED" : "FAILURES: " + std::to_string(failures))
              << "\n";
    return failures == 0 ? 0 : 1;
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------

static void usage() {
    std::cout <<
      "usage: csf <N> [options]\n"
      "       csf --selftest\n\n"
      "  --threads=T     worker threads (0 = hardware concurrency)\n"
      "  --b1=B          stage-0 trial division bound (default 1000000)\n"
      "  --no-lehman     disable the Lehman multiplier stream\n"
      "  --hf            enable the horizontal (b-driven) Fermat stream\n"
      "  --ia            enable the iterated-average GCD stream (Sec. 5/5.2)\n"
      "  --no-yp         disable the fused yellow-path stream (Sec. 4, Figure 8)\n"
        "  --parallel      run the yellow path across every thread, alone\n"
      "  --ctm           enable complex trial multiplication from the collapse\n"
      "  --digit-only    sieve with modulus 20 only (the terminal-digit sieve)\n"
      "  --no-sieve      disable sieving entirely\n"
      "  --only=S        run a single stream: vf|hf|ia|td|lm|yp|ctm\n"
      "  --wheel=R       residue budget for the sieve wheel (0 = auto)\n"
      "  --quiet         print only the factorization\n"
      "  --verbose       print stream statistics\n";
}

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty()) { usage(); return 1; }
    if (args[0] == "--selftest") return selftest();
    if (args[0] == "-h" || args[0] == "--help") { usage(); return 0; }

    Options opt;
    BigInt N;
    bool haveN = false;

    for (const auto& a : args) {
        if (a.rfind("--threads=", 0) == 0)      opt.threads = std::stoul(a.substr(10));
        else if (a.rfind("--b1=", 0) == 0)      opt.b1      = std::stoull(a.substr(5));
        else if (a == "--no-lehman")            opt.lehman  = false;
        else if (a == "--hf")                   opt.hf      = true;
        else if (a == "--ia")                   opt.ia      = true;
        else if (a == "--no-yp")                opt.yp      = false;
        // --yp is what enabled the stream before it became the default.
        // Kept as an accepted no-op so existing command lines still run.
        else if (a == "--yp")                   opt.yp      = true;
        else if (a == "--parallel")           { opt.parallel = true; opt.yp = true; }
        else if (a == "--ctm")                  opt.ctm     = true;
        else if (a == "--digit-only")           opt.digit_only = true;
        else if (a == "--no-sieve")             opt.no_sieve   = true;
        else if (a.rfind("--only=", 0) == 0)    opt.only    = a.substr(7);
        else if (a.rfind("--wheel=", 0) == 0)   opt.wheel   = std::stoull(a.substr(8));
        else if (a == "--quiet")                opt.quiet   = true;
        else if (a == "--verbose")              opt.verbose = true;
        else if (a.rfind("--", 0) == 0) { std::cerr << "unknown option: " << a << "\n"; return 1; }
        else {
            if (N.get_mpz_t()->_mp_size == 0 && !haveN) {
                if (mpz_set_str(N.get_mpz_t(), a.c_str(), 10) != 0) {
                    std::cerr << "bad integer: " << a << "\n"; return 1;
                }
                haveN = true;
            }
        }
    }
    if (!haveN) { usage(); return 1; }
    if (N < 2)  { std::cout << N << " has no prime factorization\n"; return 0; }

    // --parallel gives the yellow path the whole pool, and stands the other
    // streams down while it does.
    //
    // Standing them down is what makes it safe rather than oversubscribed:
    // T yellow threads ON TOP OF the usual race would be T+5 runnable
    // threads on T cores. It is also sound, because the yellow path is a
    // complete algorithm by itself -- its trial-division leg catches
    // p <= 2*j_max+3 and its vertical leg catches everything above, so the
    // two legs cover each other. That is the property the strip partition
    // rests on, used here in its own right.
    //
    // Resolved after the parse rather than inside it so that --parallel and
    // --only= compose in either order; an explicit --only= always wins.
    if (opt.parallel && opt.only.empty()) opt.only = "yp";

    Stats stats;
    std::vector<BigInt> fs;
    std::vector<std::string> hows;
    std::vector<BigInt> unfactored;

    const auto t0 = std::chrono::steady_clock::now();
    factor_recursive(N, opt, stats, fs, hows, unfactored);
    const auto t1 = std::chrono::steady_clock::now();

    std::sort(fs.begin(), fs.end());

    if (!opt.quiet) {
        std::cout << "N      = " << N << "\n";
        std::cout << "bits   = " << mpz_sizeinbase(N.get_mpz_t(), 2)
                  << ", digits = " << N.get_str().size()
                  << ", N mod 4 = " << umod(N, 4)
                  << ", N mod 10 = " << umod(N, 10) << "\n";
    }

    std::cout << (opt.quiet ? "" : "factors= ");
    for (std::size_t i = 0; i < fs.size(); ++i)
        std::cout << (i ? " * " : "") << fs[i];
    std::cout << "\n";

    if (!unfactored.empty()) {
        std::cerr << "INCOMPLETE: the following composite(s) were not split -- "
                     "the selected streams are not exhaustive:\n";
        for (const auto& u : unfactored) std::cerr << "  " << u << "\n";
    }

    if (!opt.quiet) {
        std::cout << "time   = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                  << " ms\n";
        if (!hows.empty()) {
            std::cout << "splits = ";
            for (std::size_t i = 0; i < hows.size(); ++i)
                std::cout << (i ? ", " : "") << hows[i];
            std::cout << "\n";
        }
    }

    if (opt.verbose) {
        std::cout << "-- stream statistics --\n";
        std::cout << "  vertical-fermat residues visited : " << stats.vf_scanned.load() << "\n";
        std::cout << "  vertical-fermat sqrt tests       : " << stats.vf_candidates.load() << "\n";
        std::cout << "  trial-division divisors tested   : " << stats.td_tested.load() << "\n";
        std::cout << "  horizontal-fermat sqrt tests     : " << stats.hf_candidates.load() << "\n";
        std::cout << "  lehman sqrt tests                : " << stats.lm_candidates.load() << "\n";
        std::cout << "  iterated-average gcds            : " << stats.ia_gcds.load() << "\n";
        std::cout << "  yellow-path steps                : " << stats.yp_steps.load() << "\n";
        std::cout << "  ctm multiplications              : " << stats.ctm_steps.load() << "\n";
    }
    return unfactored.empty() ? 0 : 2;
}
