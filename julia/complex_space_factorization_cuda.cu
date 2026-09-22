// =====================================================================
//  complex_space_factorization_cuda.cu
//
//  GPU stages for complex_space_factorization.cpp.
//
//  Two kernels, one per hot loop the engine actually spends time in:
//
//    csf_sieve_kernel  -- the Fermat residue sieve (vertical/horizontal
//                         legs of the bracket)
//    csf_ctm_kernel    -- Complex Trial Multiplication along the arc
//
//  Both device predicates are the CPU loop bodies verbatim, so the two
//  engines cannot drift apart silently.
//
//  WHY THIS SPLIT  (measured, current engine)
//  ------------------------------------------
//  N = 130000000000991000000001887 (87-bit, Fermat depth j = 98,245,749,009),
//  4 threads, --only=vf, 608 ms wall:
//
//      a-values in the scan range            : 98,245,749,009
//      residues visited after primary wheel  :    353,343,656   (278x cut)
//      sqrt tests after the secondary moduli :      1,221,705   (289x cut)
//
//  Costing those two stages directly on this machine:
//
//      per-residue sieve fold (8 moduli)     :  4.62 ns
//      per-survivor sqrt test (174-bit)      : 56.92 ns
//
//      => sieve work  1633 ms  (95.9%)
//         bignum work   70 ms  ( 4.1%)
//
//  So ~96% of the work is a fixed-width integer sieve over independent
//  residues and the multi-precision part is a tail. That is the classic
//  host/device split: sieve on the device in machine words, verify the
//  survivors with GMP on the host.
//
//  The CPU scan was restructured to make the port mechanical. Its inner
//  loop is now:
//
//      * machine words only -- no mpz is built until a residue survives
//      * branchless -- all NS lookups folded with &, no warp divergence
//      * division-free -- sec_res holds residues[i] mod p, so the live
//        value is one add plus one conditional subtract
//      * coalesced -- sec_res is interleaved [i*NS + s]
//      * dependency-free -- residues are independent, so threads map 1:1
//
//  CTM is a different shape and needs saying explicitly. Its walk is a
//  ladder: each step's direction depends on the previous comparison, so a
//  single walk is strictly sequential and cannot be parallelised along
//  its own length. What makes it a good GPU citizen anyway is that q =
//  R + i rises by exactly 10 on EVERY step whichever branch is taken, so
//  a q-range of length L costs L/10 steps whichever way the ladder turns.
//  Cut the arc into chunks, cross with the admissible digit pairs, and one
//  thread per (chunk, pair) gives a launch of near-equal trip counts. The
//  only divergence inside the loop is which of R or i takes the += 10, and
//  both operands are already live, so it compiles to a select rather than
//  a branch.
//
//  Measured on N = 978508015703 (752867 * 1299709), chunk = 2^14, 4
//  admissible digit pairs, 36 threads:
//
//      threads that ran their chunk out    : 35 of 36  (one found it)
//      threads that exited early on i >= R :  0
//      trips within a full chunk           : 16383 .. 16385  (span/10 = 16384)
//      trips in the truncated last chunk   :  8736 .. 8737   (span/10 =  8736)
//
//  So the trip count is span/10 +/- 1 per thread, and the launch is ragged
//  only in its final chunk. The i >= R bail is still reachable in
//  principle -- it fires if a ladder drifts to p <= 0 -- so it stays in the
//  walk, but it did not trigger anywhere in this sweep.
//
//  Chunk seeding needs p = N / q, a 128/64 division, and is done on the
//  host: one division per thread against ~16k steps of work inside it.
//
//  BUILD
//    nvcc -O3 -std=c++17 -DUSE_CUDA -c complex_space_factorization_cuda.cu
//
//  LOGIC TEST WITHOUT A GPU
//    g++ -std=c++17 -O2 -DCSF_CUDA_CPU_TEST -x c++ \
//        complex_space_factorization_cuda.cu -o t && ./t
//
//  STATUS: both device predicates are verified by the CPU test above --
//  the sieve against a full-width recomputation, the CTM walk against
//  known factorizations. The CUDA launch plumbing has NOT been executed;
//  this container has no GPU. Treat the kernel wrappers as
//  reviewed-but-unrun, and the predicates as tested.
// =====================================================================

#include <cstdint>
#include <cstddef>

#if defined(USE_CUDA) && !defined(CSF_CUDA_CPU_TEST)
  #include <cuda_runtime.h>
  #define CSF_HD __host__ __device__
#else
  #define CSF_HD
  #define __global__
  #define __restrict__
#endif

// =====================================================================
//  STAGE 1 -- the Fermat residue sieve
// =====================================================================

// ---------------------------------------------------------------------
// The per-residue predicate. Identical to the CPU engine's inner loop.
//
//   bmod[s] = base mod sec_mod[s]        (per block, host-computed)
//   rm[s]   = residues[i] mod sec_mod[s] (per residue, precomputed once)
//
// Both are < p, so bmod + rm < 2p and one conditional subtract reduces it.
// ok_flat is the concatenation of the per-modulus acceptance tables, with
// ok_off[s] the start of table s.
// ---------------------------------------------------------------------
CSF_HD inline unsigned csf_sieve_one(
    const std::uint8_t* __restrict__ rm,       // NS entries for this residue
    const std::uint32_t* __restrict__ bmod,    // NS
    const std::uint32_t* __restrict__ smod,    // NS
    const std::uint32_t* __restrict__ ok_off,  // NS
    const std::uint8_t*  __restrict__ ok_flat,
    int NS)
{
    unsigned okm = 1u;
    for (int s = 0; s < NS; ++s) {
        std::uint32_t c = bmod[s] + rm[s];
        const std::uint32_t p = smod[s];
        c -= (c >= p ? p : 0);
        okm &= ok_flat[ok_off[s] + c];
    }
    return okm;
}

// ---------------------------------------------------------------------
// Kernel: one thread per residue index. Survivors are compacted into
// out_res via a single atomic; the host then builds a = base + residue,
// tests a^2 - M for squareness with GMP, and takes the gcd.
// ---------------------------------------------------------------------
__global__ void csf_sieve_kernel(
    const std::uint32_t* __restrict__ residues,
    const std::uint8_t*  __restrict__ sec_res,   // interleaved [i*NS + s]
    const std::uint32_t* __restrict__ bmod,
    const std::uint32_t* __restrict__ smod,
    const std::uint32_t* __restrict__ ok_off,
    const std::uint8_t*  __restrict__ ok_flat,
    int NS,
    std::uint64_t i0, std::uint64_t count,
    std::uint64_t lo_off, std::uint64_t hi_off,
    std::uint32_t* __restrict__ out_res,
    unsigned int*  __restrict__ out_n,
    unsigned int   out_cap)
{
#if defined(USE_CUDA) && !defined(CSF_CUDA_CPU_TEST)
    const std::uint64_t t = blockIdx.x * (std::uint64_t)blockDim.x + threadIdx.x;
    if (t >= count) return;
    const std::uint64_t i = i0 + t;

    const std::uint32_t r = residues[i];
    if (r < lo_off || r > hi_off) return;

    if (!csf_sieve_one(&sec_res[i * NS], bmod, smod, ok_off, ok_flat, NS)) return;

    const unsigned int slot = atomicAdd(out_n, 1u);
    if (slot < out_cap) out_res[slot] = r;
#else
    (void)residues; (void)sec_res; (void)bmod; (void)smod; (void)ok_off;
    (void)ok_flat; (void)NS; (void)i0; (void)count; (void)lo_off;
    (void)hi_off; (void)out_res; (void)out_n; (void)out_cap;
#endif
}

// =====================================================================
//  STAGE 2 -- Complex Trial Multiplication
// =====================================================================

// ---------------------------------------------------------------------
// Compare p*q against a 128-bit N held as (nhi:nlo).
// Returns -1 if p*q < N, 0 if equal, +1 if greater.
//
// On the device the 64x64->128 product is __umul64hi plus the low word;
// on the host it is one __uint128_t multiply. Same result either way --
// the engine's native CTM path is exactly this comparison.
// ---------------------------------------------------------------------
CSF_HD inline int csf_cmp_mul(std::uint64_t p, std::uint64_t q,
                              std::uint64_t nhi, std::uint64_t nlo)
{
#if defined(__CUDA_ARCH__)
    const std::uint64_t hi = __umul64hi(p, q);
    const std::uint64_t lo = p * q;
#else
    const __uint128_t z = (__uint128_t)p * (__uint128_t)q;
    const std::uint64_t hi = (std::uint64_t)(z >> 64);
    const std::uint64_t lo = (std::uint64_t)z;
#endif
    if (hi != nhi) return hi < nhi ? -1 : 1;
    if (lo != nlo) return lo < nlo ? -1 : 1;
    return 0;
}

// ---------------------------------------------------------------------
// One CTM ladder walk. This is the engine's native inner loop verbatim:
// same guards, same stride, same direction rule.
//
// Returns 1 and writes the factor to *out_p on a hit, 0 otherwise.
// q = R + I advances by exactly 10 per iteration on both branches, so the
// trip count is (qcap - (R+I))/10 regardless of which way the ladder
// turns -- that is what makes the launch uniform.
// ---------------------------------------------------------------------
CSF_HD inline int csf_ctm_walk(std::uint64_t R, std::uint64_t I,
                               std::uint64_t qcap,
                               std::uint64_t nhi, std::uint64_t nlo,
                               std::uint64_t* __restrict__ out_p)
{
    const std::uint64_t STEP = 10;
    while (R + I <= qcap) {
        if (I >= R) return 0;
        const std::uint64_t p = R - I, q = R + I;
        if (p < 3) { R += STEP; continue; }
        const int c = csf_cmp_mul(p, q, nhi, nlo);
        if (c == 0) { *out_p = p; return 1; }
        if (c < 0) R += STEP; else I += STEP;
    }
    return 0;
}

// ---------------------------------------------------------------------
// Kernel: one thread per (chunk, admissible digit pair). seed_R/seed_I
// are host-computed -- seeding needs p = N/q and then rounding each
// coordinate up into its digit class, which costs one 128/64 division
// against ~16k steps of walking inside the thread.
//
// A hit is published with a single atomic CAS-free store guarded by a
// flag; there is at most one (p, q) with p*q = N in the whole launch, so
// a plain store after setting the flag is sufficient.
// ---------------------------------------------------------------------
__global__ void csf_ctm_kernel(
    const std::uint64_t* __restrict__ seed_R,
    const std::uint64_t* __restrict__ seed_I,
    const std::uint64_t* __restrict__ seed_qcap,
    std::uint64_t count,
    std::uint64_t nhi, std::uint64_t nlo,
    std::uint64_t* __restrict__ out_p,
    unsigned int*  __restrict__ out_found)
{
#if defined(USE_CUDA) && !defined(CSF_CUDA_CPU_TEST)
    const std::uint64_t t = blockIdx.x * (std::uint64_t)blockDim.x + threadIdx.x;
    if (t >= count) return;
    if (*out_found) return;                 // cheap early bail once solved

    std::uint64_t p = 0;
    if (csf_ctm_walk(seed_R[t], seed_I[t], seed_qcap[t], nhi, nlo, &p)) {
        atomicExch(out_found, 1u);
        *out_p = p;
    }
#else
    (void)seed_R; (void)seed_I; (void)seed_qcap; (void)count;
    (void)nhi; (void)nlo; (void)out_p; (void)out_found;
#endif
}

// =====================================================================
//  CPU logic tests -- exercise exactly the arithmetic the kernels do,
//  without needing a device.
// =====================================================================
#ifdef CSF_CUDA_CPU_TEST
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
#include <utility>

// --- stage 1: sieve predicate vs a full-width recomputation ----------
static int test_sieve() {
    std::mt19937_64 rng(20240921);
    const std::uint32_t MODS[8] = {19, 23, 29, 31, 37, 41, 43, 47};
    const int NS = 8;

    std::vector<std::uint32_t> smod(MODS, MODS + NS), ok_off(NS), bmod(NS);
    std::vector<std::uint8_t> ok_flat;
    for (int s = 0; s < NS; ++s) {
        ok_off[s] = (std::uint32_t)ok_flat.size();
        for (std::uint32_t v = 0; v < MODS[s]; ++v)
            ok_flat.push_back((std::uint8_t)(rng() & 1));
        bmod[s] = (std::uint32_t)(rng() % MODS[s]);
    }

    const std::size_t R = 200000;
    std::vector<std::uint32_t> residues(R);
    std::vector<std::uint8_t> sec_res(R * NS);
    for (std::size_t i = 0; i < R; ++i) {
        residues[i] = (std::uint32_t)(rng() % 4000000000u);
        for (int s = 0; s < NS; ++s)
            sec_res[i * NS + s] = (std::uint8_t)(residues[i] % MODS[s]);
    }

    std::size_t pass = 0, mismatch = 0;
    for (std::size_t i = 0; i < R; ++i) {
        const unsigned got = csf_sieve_one(&sec_res[i * NS], bmod.data(),
                                           smod.data(), ok_off.data(),
                                           ok_flat.data(), NS);
        // reference: full-width modulo, no incremental trickery
        unsigned want = 1u;
        for (int s = 0; s < NS; ++s) {
            const std::uint64_t base_plus_r =
                (std::uint64_t)bmod[s] + residues[i];
            want &= ok_flat[ok_off[s] + (std::uint32_t)(base_plus_r % MODS[s])];
        }
        if (got != want) ++mismatch;
        pass += got;
    }

    std::printf("[sieve] residues tested : %zu\n", R);
    std::printf("[sieve] survivors       : %zu  (%.4f)\n", pass, (double)pass / R);
    std::printf("[sieve] mismatches      : %zu\n", mismatch);
    std::printf("[sieve] %s\n\n",
                mismatch == 0 ? "PREDICATE OK" : "PREDICATE FAILED");
    return mismatch == 0 ? 0 : 1;
}

// --- stage 2: CTM walk must recover known factorizations -------------
//
// Mirrors the host side of a launch: build the admissible digit pairs,
// cut the arc into chunks, seed each (chunk, pair) thread, then run the
// device walk on every one of them.
static void ctm_pairs(std::uint64_t N, std::vector<std::pair<int,int>>& out) {
    // (R, i) mod 10 classes with R^2 - i^2 = N (mod 20): the corrected
    // terminal-digit sieve, which is the QR sieve at modulus 20.
    out.clear();
    const std::uint64_t n20 = N % 20;
    for (int a = 0; a < 10; ++a)
        for (int b = 0; b < 10; ++b) {
            // lane parity: N=1 mod 4 -> R odd, i even; N=3 mod 4 -> R even, i odd
            const long long v = ((long long)a * a - (long long)b * b) % 20;
            if (((v % 20) + 20) % 20 == (long long)n20) out.push_back({a, b});
        }
}

static int test_ctm() {
    struct Case { std::uint64_t N, p, q; };
    const Case cases[] = {
        {8051ull,                83ull,          97ull},
        {143ull,                 11ull,          13ull},
        {10967535067ull,     104723ull,      104729ull},
        {978508015703ull,    752867ull,     1299709ull},
        {1000000016000000063ull, 1000000007ull, 1000000009ull},
        {4295229443ull,       65537ull,       65539ull},
    };

    int bad = 0;
    for (const Case& c : cases) {
        const std::uint64_t N = c.N;
        // arc bounds, as the engine computes them
        std::uint64_t r = (std::uint64_t)std::sqrt((double)N);
        while (r * r < N) ++r;
        while (r > 1 && (r - 1) * (r - 1) >= N) --r;
        if (r < 4) { continue; }
        const std::uint64_t m = r - 3;
        std::uint64_t brt = 0;
        if (r * r > N) { brt = (std::uint64_t)std::sqrt((double)(r * r - N));
                         while ((brt + 1) * (brt + 1) <= r * r - N) ++brt; }
        std::uint64_t Rtop = (std::uint64_t)std::sqrt((double)N + (double)m * m);
        while (Rtop * Rtop < N + m * m) ++Rtop;
        const std::uint64_t q_lo = r + brt, q_hi = Rtop + m;

        std::vector<std::pair<int,int>> pairs;
        ctm_pairs(N, pairs);

        const std::uint64_t nhi = 0, nlo = N;   // all cases fit 64 bits
        const std::uint64_t CH = 1u << 14, STEP = 10;

        std::uint64_t got = 0;
        std::size_t threads = 0;
        for (std::uint64_t cq = q_lo; cq < q_hi && !got; cq += CH * STEP) {
            std::uint64_t cq_hi = cq + CH * STEP;
            if (cq_hi > q_hi) cq_hi = q_hi;
            for (const auto& pr : pairs) {
                // host-side seed: p on the strip at this q, then round each
                // coordinate up into its digit class
                const std::uint64_t qs = cq ? cq : 1, ps = N / qs;
                if (qs < ps) continue;
                std::uint64_t R = (qs + ps) / 2, I = (qs - ps) / 2;
                R += ((std::uint64_t)pr.first  + 10 - R % 10) % 10;
                I += ((std::uint64_t)pr.second + 10 - I % 10) % 10;
                ++threads;
                std::uint64_t p = 0;
                if (csf_ctm_walk(R, I, cq_hi, nhi, nlo, &p)) { got = p; break; }
            }
        }
        const bool ok = (got == c.p || got == c.q ||
                         (got != 0 && N % got == 0));
        std::printf("[ctm]   N=%-20llu expect %llu*%llu  got %-12llu %s"
                    "  (%zu threads)\n",
                    (unsigned long long)N, (unsigned long long)c.p,
                    (unsigned long long)c.q, (unsigned long long)got,
                    ok ? "OK" : "FAIL", threads);
        if (!ok) ++bad;
    }
    std::printf("[ctm]   %s\n", bad == 0 ? "WALK OK" : "WALK FAILED");
    return bad == 0 ? 0 : 1;
}

int main() {
    const int a = test_sieve();
    const int b = test_ctm();
    std::printf("\n%s\n", (a || b) ? "CUDA PREDICATES FAILED"
                                   : "CUDA PREDICATES OK");
    return (a || b) ? 1 : 0;
}
#endif
