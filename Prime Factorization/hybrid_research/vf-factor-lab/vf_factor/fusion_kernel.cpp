// Resumable exact Fermat cursor for the unified factorizer.
// The full terminal-digit wheel has midpoint period 10. Squares modulo 64
// have period 32, so their intersection has midpoint period lcm(10,32)=160.
// Every omitted midpoint is excluded by a necessary square congruence.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct FusionSnapshot {
    u64 n;
    u64 factor;
    u64 candidates;
    u64 square_tests;
    u64 delta_candidates;
    u64 first_midpoint;
    u64 last_midpoint;
    u64 start_midpoint;
    u64 next_midpoint;
    u64 covered_through;
    u64 half_gap;
    u64 root_ceiling;
    u64 endpoint;
    std::uint32_t residue_count;
    std::uint32_t period;
    std::uint32_t status;
    std::uint32_t exhausted;
    double elapsed_seconds;
    double total_kernel_seconds;
};

namespace {
using Clock = std::chrono::steady_clock;
enum Status : unsigned {
    READY = 0, FACTOR = 1, EXHAUSTED = 2,
    CANDIDATE_BUDGET = 3, TIME_BUDGET = 4, RANGE_LIMIT = 5
};

unsigned bit_length(u64 n) {
    return n ? 64U - static_cast<unsigned>(__builtin_clzll(n)) : 0U;
}

u64 isqrt64(u64 n) {
    if (n < 2) return n;
    // Hardware sqrt supplies only an estimate. The integer inequalities below
    // determine the answer exactly, including when double(n) rounds to 2^64.
    // The estimate is at most 2^32, so x+1 and both uint128 products are safe.
    u64 x = static_cast<u64>(std::sqrt(static_cast<double>(n)));
    while (u128(x) * x > n) --x;
    while (u128(x + 1) * (x + 1) <= n) ++x;
    return x;
}

u64 isqrt128(u128 n) {
    if (n <= std::numeric_limits<u64>::max()) return isqrt64(static_cast<u64>(n));
    const unsigned bits = 64 + bit_length(static_cast<u64>(n >> 64));
    u128 x = u128(1) << ((bits + 1) / 2);
    for (;;) {
        const u128 y = (x + n / x) / 2;
        if (y >= x) return static_cast<u64>(x);
        x = y;
    }
}

constexpr u64 square_mask64() {
    u64 mask = 0;
    for (unsigned b = 0; b < 64; ++b) mask |= u64(1) << ((b*b) & 63);
    return mask;
}

struct Cursor {
    u64 n, a, covered, root_ceiling, endpoint, start;
    u128 difference;
    u64 factor = 0, candidates = 0, first = 0, last = 0, half_gap = 0;
    unsigned residues[160] = {}, gaps[160] = {}, count = 0, index = 0;
    unsigned terminal = READY;
    double total_seconds = 0;

    explicit Cursor(u64 number) : n(number) {
        const auto started = Clock::now();
        unsigned midpoint_mask = 0;
        for (unsigned r = 0; r < 10; ++r) {
            for (unsigned b = 0; b < 10; ++b) {
                if ((r*r + 100 - b*b) % 20 == n % 20) {
                    midpoint_mask |= 1U << r;
                    break;
                }
            }
        }
        static constexpr u64 qr_mask = square_mask64();
        for (unsigned r = 0; r < 160; ++r) {
            const unsigned difference64 = (r*r + 64 - (n & 63)) & 63;
            if (((midpoint_mask >> (r % 10)) & 1U) &&
                ((qr_mask >> difference64) & 1U)) {
                residues[count++] = r;
            }
        }
        // For odd n the trivial representation guarantees a nonempty wheel.
        if (!count) throw std::bad_alloc();
        for (unsigned j = 0; j < count; ++j) {
            gaps[j] = j + 1 < count ? residues[j+1] - residues[j]
                                   : 160 + residues[0] - residues[j];
        }
        root_ceiling = isqrt64(n);
        root_ceiling += u128(root_ceiling) * root_ceiling != n;
        a = root_ceiling;
        const unsigned rem = a % 160;
        while (index < count && residues[index] < rem) ++index;
        if (index == count) {
            a += 160 - rem + residues[0];
            index = 0;
        } else {
            a += residues[index] - rem;
        }
        start = a;
        endpoint = static_cast<u64>((u128(n) + 1) / 2);
        // This initial interval contains only congruence-rejected midpoints.
        covered = std::min(a - 1, endpoint);
        difference = u128(a) * a - n;
        total_seconds = std::chrono::duration<double>(Clock::now()-started).count();
    }

    void fill(FusionSnapshot& out, unsigned status, u64 previous, double elapsed) const {
        out = {};
        out.n = n;
        out.factor = factor;
        out.candidates = candidates;
        out.square_tests = candidates;
        out.delta_candidates = candidates - previous;
        out.first_midpoint = first;
        out.last_midpoint = last;
        out.start_midpoint = start;
        out.next_midpoint = a;
        out.covered_through = covered;
        out.half_gap = half_gap;
        out.root_ceiling = root_ceiling;
        out.endpoint = endpoint;
        out.residue_count = count;
        out.period = 160;
        out.status = status;
        out.exhausted = terminal == EXHAUSTED;
        out.elapsed_seconds = elapsed;
        out.total_kernel_seconds = total_seconds;
    }

    void advance(u64 maximum, double seconds, u64 midpoint_limit, FusionSnapshot& out) {
        const auto began = Clock::now();
        const u64 before = candidates;
        auto elapsed = [&]() { return std::chrono::duration<double>(Clock::now()-began).count(); };
        auto finish = [&](unsigned reason) {
            const double duration = elapsed();
            total_seconds += duration;
            fill(out, reason, before, duration);
        };
        if (terminal != READY) { finish(terminal); return; }
        if (maximum == 0) { finish(CANDIDATE_BUDGET); return; }
        if (seconds == 0.0) { finish(TIME_BUDGET); return; }
        const u64 stop = std::min(midpoint_limit, endpoint);
        while (a <= stop) {
            const u64 delta = candidates - before;
            if (delta >= maximum) { finish(CANDIDATE_BUDGET); return; }
            if ((delta & 255U) == 0 && seconds >= 0.0 && elapsed() >= seconds) {
                finish(TIME_BUDGET); return;
            }
            ++candidates;
            if (candidates == 1) first = a;
            last = a;
            const u64 b = isqrt128(difference);
            if (u128(b) * b == difference) {
                const u64 divisor = a - b;
                if (divisor > 1 && divisor < n) {
                    factor = divisor;
                    half_gap = b;
                    terminal = FACTOR;
                    // Successful midpoint a is not part of excluded coverage.
                    finish(FACTOR);
                    return;
                }
                if (divisor == 1) {
                    // The trivial endpoint cannot conceal a proper factor.
                    covered = endpoint;
                    a = endpoint + 1;
                    terminal = EXHAUSTED;
                    finish(EXHAUSTED);
                    return;
                }
            }
            const u64 gap = gaps[index];
            difference += u128(2) * a * gap + u128(gap) * gap;
            a += gap;
            if (++index == count) index = 0;
            // The tested midpoint failed, and every following skipped point
            // fails a necessary congruence. This certificate may extend beyond
            // the caller's temporary midpoint limit without any extra tests.
            covered = std::min(a - 1, endpoint);
        }
        if (a > endpoint) {
            terminal = EXHAUSTED;
            finish(EXHAUSTED);
        } else {
            finish(RANGE_LIMIT);
        }
    }
};
} // namespace

extern "C" {
unsigned fusion_abi_version() { return 1; }
unsigned fusion_snapshot_size() { return sizeof(FusionSnapshot); }
const char* fusion_compiler() { return __VERSION__; }
u64 fusion_isqrt(u64 n) { return isqrt64(n); }

void* fusion_fermat_create(u64 n, unsigned* error) {
    if (error) *error = 0;
    if (n < 3 || !(n & 1U)) {
        if (error) *error = 1;
        return nullptr;
    }
    try { return new Cursor(n); }
    catch (...) { if (error) *error = 2; return nullptr; }
}

void fusion_fermat_destroy(void* pointer) { delete static_cast<Cursor*>(pointer); }

int fusion_fermat_snapshot(void* pointer, FusionSnapshot* out) {
    if (!pointer || !out) return 1;
    auto* cursor = static_cast<Cursor*>(pointer);
    cursor->fill(*out, cursor->terminal, cursor->candidates, 0.0);
    return 0;
}

int fusion_fermat_advance(void* pointer, u64 maximum, double seconds,
                          u64 midpoint_limit, FusionSnapshot* out) {
    if (!pointer || !out) return 1;
    try {
        static_cast<Cursor*>(pointer)->advance(maximum, seconds, midpoint_limit, *out);
        return 0;
    } catch (...) { return 2; }
}
}
