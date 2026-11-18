// Complete Fermat terminal-digit sieve, with matched native controls.
// Index the pair table by N mod 20, not N mod 10 alone. R and i have
// period ten in the table because (x+10)^2 == x^2 (mod 20).
// All square roots, square comparisons, and incremental updates are exact.
#include <chrono>
#include <cstdint>
#include <initializer_list>
#include <limits>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct WheelResult {
    u64 factor;
    u64 candidates;
    u64 square_tests;
    u64 first_midpoint;
    u64 last_midpoint;
    u64 start_midpoint;
    u64 root_ceiling;
    u64 half_gap;
    std::uint32_t residue_mask;
    std::uint32_t residue_count;
    std::uint32_t status; // 0 budget, 1 factor, 2 prime, 3 error
    std::uint32_t error;  // 1 time budget, 2 candidate budget
    double elapsed_seconds;
};

namespace {
using Clock = std::chrono::steady_clock;
struct Stopped {};

unsigned bit_length(u64 n) {
    return n ? 64U - static_cast<unsigned>(__builtin_clzll(n)) : 0U;
}

u64 isqrt64(u64 n) {
    if (n < 2) return n;
    u64 x = u64(1) << ((bit_length(n) + 1) / 2);
    for (;;) {
        const u64 y = (x + n / x) / 2;
        if (y >= x) return x;
        x = y;
    }
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
    for (unsigned b = 0; b < 64; ++b) mask |= u64(1) << ((b * b) & 63);
    return mask;
}

struct Context {
    WheelResult& out;
    Clock::time_point started;
    double seconds;
    u64 maximum;

    Context(WheelResult& result, double timeout, u64 cap)
        : out(result), started(Clock::now()), seconds(timeout), maximum(cap) {}

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - started).count();
    }

    void work() {
        if (out.candidates >= maximum) {
            out.error = 2;
            throw Stopped{};
        }
        if ((out.candidates & 1023U) == 0 && seconds >= 0.0 && elapsed() >= seconds) {
            out.error = 1;
            throw Stopped{};
        }
        ++out.candidates;
    }
};

void search(u64 n, bool wheel, bool qr64, Context& ctx) {
    unsigned residues[10] = {}, gaps[10] = {}, count = 0;
    unsigned period;
    if (wheel) {
        period = 10;
        const unsigned n20 = n % 20;
        // Construct the exact projection of the complete pair table. The
        // entire setup is inside the measured call, with no per-N cache.
        for (unsigned r = 0; r < 10; ++r) {
            bool allowed = false;
            for (unsigned b = 0; b < 10; ++b) {
                if ((r * r + 100 - b * b) % 20 == n20) allowed = true;
            }
            if (allowed) {
                residues[count++] = r;
                ctx.out.residue_mask |= 1U << r;
            }
        }
    } else {
        period = 2;
        residues[count++] = n % 4 == 1 ? 1 : 0;
        ctx.out.residue_mask = 1U << residues[0];
    }
    ctx.out.residue_count = count;
    if (count == 0) {
        ctx.out.status = 3;
        ctx.out.error = 3;
        return;
    }
    for (unsigned j = 0; j < count; ++j) {
        gaps[j] = j + 1 < count ? residues[j + 1] - residues[j]
                              : period + residues[0] - residues[j];
    }

    u64 a = isqrt64(n);
    a += u128(a) * a != n;
    ctx.out.root_ceiling = a;
    unsigned index = 0;
    const unsigned rem = a % period;
    while (index < count && residues[index] < rem) ++index;
    if (index == count) {
        a += period - rem + residues[0];
        index = 0;
    } else {
        a += residues[index] - rem;
    }
    ctx.out.start_midpoint = a;
    u128 difference = u128(a) * a - n;
    const u64 last = static_cast<u64>((u128(n) + 1) / 2);
    static constexpr u64 qr_mask = square_mask64();

    while (a <= last) {
        ctx.work();
        if (ctx.out.candidates == 1) ctx.out.first_midpoint = a;
        ctx.out.last_midpoint = a;
        if (!qr64 || ((qr_mask >> static_cast<unsigned>(difference & 63)) & 1U)) {
            ++ctx.out.square_tests;
            const u64 b = isqrt128(difference);
            if (u128(b) * b == difference) {
                const u64 factor = a - b;
                if (factor > 1 && factor < n) {
                    ctx.out.factor = factor;
                    ctx.out.half_gap = b;
                    ctx.out.status = 1;
                    return;
                }
                // The trivial representation is the endpoint. For direct C
                // callers this complete ascending search proves primality.
                if (factor == 1) {
                    ctx.out.status = 2;
                    return;
                }
            }
        }
        const u64 gap = gaps[index];
        difference += u128(2) * a * gap + u128(gap) * gap;
        a += gap;
        if (++index == count) index = 0;
    }
    ctx.out.status = 2;
}
} // namespace

extern "C" {
unsigned wheel20_abi_version() { return 1; }
unsigned wheel20_result_size() { return sizeof(WheelResult); }
const char* wheel20_compiler() { return __VERSION__; }
u64 wheel20_isqrt(u64 n) { return isqrt64(n); }

int wheel20_run(u64 n, unsigned mode, double timeout, u64 maximum, WheelResult* out) {
    if (!out) return 1;
    *out = {};
    Context ctx(*out, timeout, maximum);
    if (n < 2 || mode > 3) {
        out->status = 3;
        out->error = 3;
        return 1;
    }
    // Match the shared prechecks for direct C callers. The Python wrapper
    // additionally checks primality deterministically below 2^64.
    for (u64 p : {u64(2), u64(3), u64(5)}) {
        if (n == p) { out->status = 2; return 0; }
        if (n % p == 0) { out->factor = p; out->status = 1; return 0; }
    }
    const u64 root = isqrt64(n);
    if (u128(root) * root == n) { out->factor = root; out->status = 1; return 0; }
    try {
        // mode bits: bit 0 selects the complete wheel; bit 1 enables QR64.
        search(n, mode & 1U, mode & 2U, ctx);
    } catch (const Stopped&) {
        out->status = 0;
    } catch (...) {
        out->status = 3;
        out->error = 4;
    }
    out->elapsed_seconds = ctx.elapsed();
    return out->status == 3 ? 1 : 0;
}
}
