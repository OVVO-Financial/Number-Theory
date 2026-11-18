// Optional native64 fairness baselines. All arithmetic is exact.
// No reference factors, floating-point square roots, or probable-prime
// assumptions enter this kernel. The Python wrapper applies shared prechecks.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128;

struct NativeResult {
    u64 factor;
    u64 steps;
    u64 gcd_calls;
    u64 modular_multiplications;
    u64 square_tests;
    u64 rounds;
    u64 restarts;
    u64 trivial_gcd_hits;
    u64 auxiliary;
    u64 argument_abs;
    std::uint32_t status;             // 0 budget, 1 factor, 2 prime, 3 error
    std::uint32_t selected_level;
    std::int32_t direction;           // -1 descending, +1 ascending
    std::uint32_t offset;
    std::uint32_t argument_negative;
    std::uint32_t error;
    double elapsed_seconds;
};

namespace {
using Clock = std::chrono::steady_clock;
struct Stopped {};

unsigned bit_length(u64 x) {
    return x ? 64U - static_cast<unsigned>(__builtin_clzll(x)) : 0U;
}

u64 isqrt64(u64 n) {
    if (n < 2) return n;
    u64 x = u64(1) << ((bit_length(n) + 1) / 2);
    while (true) {
        const u64 y = (x + n / x) / 2;
        if (y >= x) return x;
        x = y;
    }
}

u64 isqrt128(u128 n) {
    if (n <= std::numeric_limits<u64>::max()) return isqrt64(static_cast<u64>(n));
    const unsigned bits = 64 + bit_length(static_cast<u64>(n >> 64));
    u128 x = u128(1) << ((bits + 1) / 2);
    while (true) {
        const u128 y = (x + n / x) / 2;
        if (y >= x) return static_cast<u64>(x);
        x = y;
    }
}

u64 ceil_sqrt(u64 n) {
    const u64 r = isqrt64(n);
    return r + (u128(r) * r != n);
}

struct Context {
    NativeResult& out;
    Clock::time_point began;
    double seconds;
    u64 maximum;

    Context(NativeResult& r, double timeout, u64 max_steps)
        : out(r), began(Clock::now()), seconds(timeout), maximum(max_steps) {}

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - began).count();
    }

    void checkpoint() {
        if (seconds >= 0.0 && elapsed() >= seconds) {
            out.error = 1; // time budget
            throw Stopped{};
        }
    }

    inline void work() {
        if (out.steps >= maximum) {
            out.error = 2; // step budget
            throw Stopped{};
        }
        if ((out.steps & 1023U) == 0) checkpoint();
        ++out.steps;
    }

    void found(u64 factor) {
        out.factor = factor;
        out.status = 1;
    }
};

void trial(u64 n, Context& ctx) {
    static constexpr u64 increments[8] = {4, 2, 4, 2, 4, 6, 2, 6};
    const u64 limit = isqrt64(n);
    u64 candidate = 7;
    unsigned index = 0;
    while (candidate <= limit) {
        ctx.work();
        if (n % candidate == 0) {
            ctx.found(candidate);
            return;
        }
        candidate += increments[index];
        index = (index + 1) & 7U;
    }
    ctx.out.status = 2;
}

constexpr u64 square_mask() {
    u64 mask = 0;
    for (unsigned i = 0; i < 64; ++i) mask |= u64(1) << ((i * i) & 63);
    return mask;
}

void fermat_search(u64 n, Context& ctx) {
    static constexpr u64 mask = square_mask();
    u64 a = ceil_sqrt(n);
    const unsigned parity = n % 4 == 1 ? 1 : 0;
    a += ((a & 1U) != parity);
    u128 b_squared = u128(a) * a - n;
    const u64 last = static_cast<u64>((u128(n) + 1) / 2);
    while (a <= last) {
        ctx.work();
        if ((mask >> static_cast<unsigned>(b_squared & 63)) & 1U) {
            ++ctx.out.square_tests;
            const u64 b = isqrt128(b_squared);
            if (u128(b) * b == b_squared) {
                const u64 factor = a - b;
                if (factor > 1 && factor < n) {
                    ctx.out.argument_abs = a;
                    ctx.out.auxiliary = b;
                    ctx.found(factor);
                    return;
                }
                if (factor == 1) {
                    ctx.out.status = 2;
                    return;
                }
            }
        }
        b_squared += 4 * u128(a) + 4;
        a += 2;
    }
    ctx.out.status = 2;
}

inline u64 advance(u64 value, u64 c, u64 n, Context& ctx) {
    ctx.work();
    ++ctx.out.modular_multiplications;
    return static_cast<u64>((u128(value) * value + c) % n);
}

void brent(u64 n, Context& ctx, u64 seed, unsigned batch_size) {
    std::mt19937_64 rng(seed);
    while (true) {
        u64 y = rng() % (n - 2) + 1;
        const u64 c = rng() % (n - 2) + 1;
        u64 cycle = 1, g = 1, x = 0, saved_y = y;
        while (g == 1) {
            x = y;
            for (u64 i = 0; i < cycle; ++i) y = advance(y, c, n, ctx);
            u64 consumed = 0;
            while (consumed < cycle && g == 1) {
                saved_y = y;
                u64 product = 1;
                const u64 take = std::min<u64>(batch_size, cycle - consumed);
                for (u64 i = 0; i < take; ++i) {
                    y = advance(y, c, n, ctx);
                    const u64 difference = x >= y ? x - y : y - x;
                    product = static_cast<u64>(u128(product) * difference % n);
                    ++ctx.out.modular_multiplications;
                }
                g = std::gcd(product, n);
                ++ctx.out.gcd_calls;
                consumed += take;
            }
            if (cycle > (u64(1) << 62)) { g = n; break; }
            cycle *= 2;
        }
        if (g == n) {
            ++ctx.out.trivial_gcd_hits;
            g = 1;
            while (g == 1) {
                saved_y = advance(saved_y, c, n, ctx);
                const u64 difference = x >= saved_y ? x - saved_y : saved_y - x;
                g = std::gcd(difference, n);
                ++ctx.out.gcd_calls;
                ++ctx.out.auxiliary; // recovery steps
            }
        }
        if (g > 1 && g < n) {
            ctx.found(g);
            return;
        }
        ++ctx.out.restarts;
    }
}

bool vf_gcd(i128 argument, u64 n, Context& ctx, u64 distance,
            int direction, unsigned offset) {
    ctx.work();
    const bool negative = argument < 0;
    const u128 magnitude = negative ? u128(-argument) : u128(argument);
    const u64 value = magnitude <= std::numeric_limits<u64>::max()
        ? static_cast<u64>(magnitude) : static_cast<u64>(magnitude % n);
    const u64 g = std::gcd(value, n);
    ++ctx.out.gcd_calls;
    if (g > 1 && g < n) {
        ctx.found(g);
        ctx.out.argument_abs = static_cast<u64>(magnitude);
        ctx.out.argument_negative = negative;
        ctx.out.direction = direction;
        ctx.out.offset = offset;
        ctx.out.auxiliary = distance;
        return true;
    }
    ctx.out.trivial_gcd_hits += g == n;
    return false;
}

void published_vf(u64 n, Context& ctx) {
    const u64 seed = ceil_sqrt(n);
    const u64 max_imaginary_floor = n >= 9 ? (n - 9) / 6 : 0;
    const unsigned average_count = std::max(1U, bit_length(max_imaginary_floor));
    const unsigned cut = std::max(2U, bit_length(average_count - 1));
    unsigned best = 1;
    u64 best_remainder = 0;
    for (unsigned k = 1; k <= cut; ++k) {
        const u64 remainder = ((seed - n) & ((u64(1) << k) - 1)) << (cut - k);
        if (k == 1 || remainder > best_remainder) {
            best = k;
            best_remainder = remainder;
        }
    }
    best = std::max(2U, best);
    ctx.out.selected_level = best;
    const u64 denominator = u64(1) << best;
    const u128 numerator = u128(denominator - 1) * n + seed;
    const u64 floor_average = static_cast<u64>(numerator / denominator);
    const u64 ceil_average = floor_average + (numerator % denominator != 0);
    const unsigned offsets = n % 10 == 1 || n % 10 == 9 ? 3 : 2;
    const u64 rounds = 2 * isqrt64(n) + 2;
    for (u64 r = 0; r < rounds; ++r) {
        ++ctx.out.rounds;
        const u64 distance = 5 * r;
        ctx.out.auxiliary = distance;
        const i128 descending = i128(floor_average) - distance;
        for (unsigned offset = 0; offset < offsets; ++offset) {
            if (vf_gcd(descending - offset, n, ctx, distance, -1, offset)) return;
        }
        const u128 ascending = u128(ceil_average) + distance;
        if (ascending < u128(n) - 2) {
            for (unsigned offset = 0; offset < offsets; ++offset) {
                if (vf_gcd(static_cast<i128>(ascending + offset), n, ctx,
                           distance, +1, offset)) return;
            }
        }
    }
    ctx.out.status = 2;
}
} // namespace

extern "C" {
unsigned vf_native_abi_version() { return 1; }
unsigned vf_native_result_size() { return sizeof(NativeResult); }
const char* vf_native_compiler() { return __VERSION__; }
u64 vf_native_isqrt(u64 n) { return isqrt64(n); }

int vf_native_run(u64 n, unsigned algorithm, double timeout_seconds,
                  u64 max_steps, u64 seed, unsigned batch_size, NativeResult* out) {
    if (!out) return 1;
    *out = {};
    Context ctx(*out, timeout_seconds, max_steps);
    if (n < 2 || batch_size == 0 || algorithm > 3) {
        out->status = 3;
        out->error = 3;
        return 1;
    }
    // These checks protect direct C callers; the Python wrapper has the same
    // checks and additionally its shared exact <=64-bit primality precheck.
    for (u64 p : {u64(2), u64(3), u64(5)}) {
        if (n == p) { out->status = 2; return 0; }
        if (n % p == 0) { ctx.found(p); return 0; }
    }
    const u64 root = isqrt64(n);
    if (u128(root) * root == n) { ctx.found(root); return 0; }
    try {
        ctx.checkpoint();
        switch (algorithm) {
            case 0: published_vf(n, ctx); break;
            case 1: brent(n, ctx, seed, batch_size); break;
            case 2: trial(n, ctx); break;
            case 3: fermat_search(n, ctx); break;
        }
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
