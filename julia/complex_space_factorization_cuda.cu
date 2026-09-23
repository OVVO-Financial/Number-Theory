// =====================================================================
//  complex_space_factorization_cuda.cu
//
//  GPU stages for complex_space_factorization.cpp.
//
//  Three kernels, one per hot loop the engine actually spends time in:
//
//    csf_sieve_kernel  -- the Fermat residue sieve (vertical/horizontal
//                         legs of the bracket)
//    csf_ctm_kernel    -- Complex Trial Multiplication along the arc
//    csf_yp_kernel     -- the yellow path, the 135-degree traversal
//
//  Each device predicate is the CPU loop body verbatim, so the two engines
//  cannot drift apart silently.
//
//  LAUNCH THEM CONCURRENTLY -- this is where the device beats the host.
//
//  On a CPU the three compete: with 4 cores and five streams wanted (td,
//  vf, lm, yp, ctm) there is no surplus to share, and the host engine has
//  to deal spare threads round-robin just to stop one stream starving
//  another. A device has no such problem. Measured per-thread floors, in
//  the DEVICE forms of each walk (112-bit N):
//
//    sieve   no seed at all -- residues are independent, one thread each
//    yellow  seed 42.50 ns / step 12.50 ns =  3.4 steps -> floor  31 steps
//    ctm     seed 47.04 ns / step  1.07 ns = 44.0 steps -> floor 396 steps
//
//  (The yellow figure is the device form, which never builds V or H; the
//  host's own bracket start is 211 ns against an 87 ns mpz step, so its
//  floor is ~120. Do not carry the host numbers over.)
//
//  Applying those floors at 112 bits, where sqrt(N) ~ 7.2e16:
//
//    yellow path   j_max / 31       ~ 5.8e14 brackets
//    ctm           arc / (396 * 10) ~ 1.8e13 chunks
//
//  Both are orders of magnitude past any real grid, so all three kernels
//  can be resident and wide simultaneously -- three concurrent launches on
//  separate CUDA streams, sized independently, with none of them starving
//  the others. Every ceiling scales as sqrt(N), so the headroom widens as
//  the problem gets harder.
//
//  Sizing rule per kernel: threads = min(grid you want, span / floor).
//  Going wider is not wrong, just increasingly seed-bound.
//
//  The yellow path is the most parallel of the three: every quantity is a
//  closed form in j, so a bracket starts anywhere for O(1) cost and there
//  is no dependency chain at all. CTM's walks are ladders -- each step's
//  direction depends on the previous comparison -- so they parallelise
//  across chunks but not along one chunk's length.
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
//  HOW WIDE CAN THE CTM LAUNCH GO
//  ------------------------------
//  Splitting the arc into more starting points is work-conserving, which is
//  what makes an arbitrarily wide launch legitimate rather than wasteful.
//  Sweeping a whole arc and varying only the seed count, a 16,384x widening
//  cost 1.5% more work: 64 threads did 556,846,579 multiplications,
//  1,048,576 threads did 565,422,794.
//
//  The floor is per-thread, not global, and it belongs to the CURRENT walk
//  -- the two-directional R/i form, not the old single ladder. Measured at
//  112 bits:
//
//      chunk seed (mpz divide + halves + digit sync) : 47.04 ns
//      one step   (mul + cmp + add, native)          :  1.07 ns
//
//  so a seed costs 44 steps and a thread needs ~396 of them -- about 3,960
//  of arc in q -- to hold seeding under 10% of its own time. The earlier
//  figure here was 440, taken when the seed was a native 5.8-8.7 ns; the
//  seed is an mpz divide now, so the floor moved by an order of magnitude.
//
//  At 112 bits the descending arc alone is ~sqrt(N) ~ 7.2e16, giving
//  ~1.8e13 chunks -- far past any real grid. The ceiling scales as sqrt(N),
//  so the headroom widens as the problem gets harder.
//
//  Sizing rule for a launch: threads = min(grid you want, arc_span / 3960).
//  Going wider than that is not wrong, just increasingly seed-bound.
//
//
//  The sieve kernel has no equivalent floor -- its per-chunk setup is NS
//  modulos (8) amortised over the whole residue block, so one thread per
//  residue is always the right mapping.
//
//  BUILD
//    nvcc -O3 -std=c++17 -DUSE_CUDA -c complex_space_factorization_cuda.cu
//
//  LOGIC TEST WITHOUT A GPU
//    g++ -std=c++17 -O2 -DCSF_CUDA_CPU_TEST -x c++
//        complex_space_factorization_cuda.cu -o t && ./t
//    (one line; split here only to fit the margin)
//
//  STATUS: all three device predicates are verified by the CPU test above
//  -- the sieve against a full-width recomputation (200,000 residues), the
//  CTM walk against six known factorizations, csf_mod128 against a
//  __uint128_t reference (400,009 cases, 350,000 of them exercising the
//  shift-subtract branch), and the yellow path against a from-scratch
//  recomputation of every j plus a check that the j carrying the factor
//  survives. The CUDA launch plumbing has NOT been executed; this
//  container has no GPU. Treat the kernel wrappers as reviewed-but-unrun,
//  and the predicates as tested.
// =====================================================================

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <vector>
#include <utility>

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
// One CTM walk. This is the reference implementation's rule, on R and i:
//
//     ascending :  if p*q > n  ->  i += 10   else  R += 10
//     descending:  if p*q < n  ->  i -= 10   else  R -= 10
//
// The moves are on R and i by 10 -- NOT on p and q independently. Both
// branches of ascending raise q = R + i by 10 and both branches of
// descending lower it by 10, so q is monotone in each walk: up in one,
// down in the other. That is the ascend and the descend, and it is what
// makes the trip count (|q_limit - q_start|)/10 and the launch uniform.
//
// Both walks start at the SAME point -- the seed csf_ctm_build_seeds
// handed this lane, which is a point on the red/blue boundary reached by
// climbing vertically from a real value in [sqrt N, sqrt 2N]. `up`
// selects which of the two directions this thread runs. It is NOT the
// first red point of the 135-degree traversal: that sits at
// p = N/S ~ 0.5 sqrt(N), one end of the strip, and seeding there put p
// below the answer with no way back up (8051 = 83*97 was unreachable).
//
// There is deliberately no p < 3 early exit: p*q < n when p is tiny, so
// the rule moves i and raises p by 10 again. Breaking there discards the
// class before it recovers, which loses N = 143 (its seed passes through
// p = 1 one step before p = 11).
//
// Returns 1 and writes the factor to *out_p on a hit, 0 otherwise.
// ---------------------------------------------------------------------
CSF_HD inline int csf_ctm_walk(std::int64_t R, std::int64_t I,
                               std::int64_t q_limit, int up,
                               std::uint64_t nhi, std::uint64_t nlo,
                               std::uint64_t* __restrict__ out_p)
{
    const std::int64_t STEP = 10;
    while (I >= 0 && R > I) {
        const std::int64_t q = R + I;
        if (up ? (q > q_limit) : (q < q_limit)) return 0;
        const std::int64_t p = R - I;
        const int c = csf_cmp_mul((std::uint64_t)p, (std::uint64_t)q, nhi, nlo);
        if (c == 0) { *out_p = (std::uint64_t)p; return 1; }
        if (up) { if (c > 0) I += STEP; else R += STEP; }
        else    { if (c < 0) I -= STEP; else R -= STEP; }
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
    const std::int64_t* __restrict__ seed_R,
    const std::int64_t* __restrict__ seed_I,
    const std::int64_t* __restrict__ seed_qlimit,
    const int*          __restrict__ seed_up,
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
    if (csf_ctm_walk(seed_R[t], seed_I[t], seed_qlimit[t], seed_up[t],
                     nhi, nlo, &p)) {
        atomicExch(out_found, 1u);
        *out_p = p;
    }
#else
    (void)seed_R; (void)seed_I; (void)seed_qlimit; (void)seed_up;
    (void)count; (void)nhi; (void)nlo; (void)out_p; (void)out_found;
#endif
}

// ---------------------------------------------------------------------
// Host-side CTM seed builder -- one launch's worth of starting points.
//
// A CTM seed is NOT a point on the 135-degree strip. It is a point on the
// red/blue boundary, entered from the real axis, and the stretch of real
// axis worth entering from is
//
//      R in [ ceil(sqrt N), floor(sqrt 2N) ]      (1 .. 1.41421 sqrt N)
//
// Everything about that interval is forced. With R = (p+q)/2 and
// i = (q-p)/2 on the hyperbola, i = sqrt(R^2 - N) and so
//
//      R = 1.00000 sqrt N  ->  p = 1.00000 sqrt N, q = 1.00000 sqrt N
//      R = 1.06066 sqrt N  ->  p = 0.70711 sqrt N, q = 1.41421 sqrt N
//      R = 1.41421 sqrt N  ->  p = 0.41421 sqrt N, q = 2.41421 sqrt N
//
// -- the whole balanced arc, both factors, nothing outside it, with
// sqrt(2N) as its top exactly.
//
// Standing at a real value R and climbing vertically, R fixed and i
// rising, the product R^2 - i^2 falls monotonically: red at i = 0 (since
// R >= sqrt N), blue past the crossing. The crossing
//
//      i = floor(sqrt(R^2 - N))           last red point at this R
//
// is the seed. N = 309 at R = 21: 441 - 11^2 = 320 is red, 441 - 12^2 =
// 297 is blue, so the seed is 21 + 11i. R and i ARE the coordinates, so
// it costs one isqrt and no division.
//
// The seeds cannot be spaced evenly in R. The walk out of a seed is paid
// in q, and
//
//      dq/dR = 1 + R / i          with i = sqrt(R^2 - N)
//
// diverges at the balanced corner where i -> 0. At 90 bits, moving the
// real part by ONE integer at the bottom of the interval jumps q by
// 5.2e6 -- 518,000 walk steps in that single thread -- while a seed at
// the top moves q by 2.41 and does nothing. Evenly spaced in R, one lane
// carries the arc and the rest of the warp idles on it. So: take every
// integer real value while they are few (309 -> 7 seeds, 8051 -> 37),
// and past that space them for equal work, which means equal in q. Those
// are real values on the same interval with the same chasm under them --
// R = (q + N/q)/2, i = (q - N/q)/2 is that same crossing written from
// the other side -- sampled densely where the arc turns and sparsely
// where it runs straight. Equal work per lane is what makes the launch
// uniform.
// ---------------------------------------------------------------------
struct CsfCtmSeeds {
    std::vector<std::int64_t> R, I, qlimit;
    std::vector<int>          up;
    std::uint64_t             chunks = 0;    // seeds before the class fan-out
    bool                      dense  = false;
};

static std::uint64_t csf_isqrt128(__uint128_t n) {
    if (n == 0) return 0;
    std::uint64_t x = (std::uint64_t)std::sqrt((double)(long double)n);
    if (x == 0) x = 1;
    for (int k = 0; k < 8; ++k) {            // Newton, then exact fixup
        const std::uint64_t y = (std::uint64_t)((x + (std::uint64_t)(n / x)) / 2);
        if (y == x) break;
        x = y;
    }
    while (x > 0 && (__uint128_t)x * x > n) --x;
    while ((__uint128_t)(x + 1) * (x + 1) <= n) ++x;
    return x;
}

// Emits (R, I, qlimit, up) for every (seed, admissible class, direction).
// `max_chunks` bounds one launch; the arc above the interval runs to
// q = N/3 and is far longer than any single launch, so the caller walks
// `first_chunk` forward across waves.
static void csf_ctm_build_seeds(std::uint64_t N,
                                const std::vector<std::pair<int,int>>& cls,
                                CsfCtmSeeds& s,
                                std::uint64_t first_chunk = 0,
                                std::uint64_t max_chunks  = 4096,
                                std::uint64_t chunk_q     = 81920ull)
{
    s.R.clear(); s.I.clear(); s.qlimit.clear(); s.up.clear();
    s.chunks = 0; s.dense = false;
    if (N < 16 || cls.empty()) return;

    std::uint64_t r = csf_isqrt128(N);
    if ((__uint128_t)r * r < (__uint128_t)N) ++r;            // ceil(sqrt N)
    if (r < 4) return;

    const std::uint64_t q_bot = r;
    const std::uint64_t q_top = N / 3;
    const std::uint64_t R_lo  = r;
    std::uint64_t R_hi = csf_isqrt128((__uint128_t)N * 2);   // 1.41421 sqrt N
    if (R_hi < R_lo) R_hi = R_lo;

    const std::uint64_t qA_hi =
        R_hi + csf_isqrt128((__uint128_t)R_hi * R_hi - N);

    const std::uint64_t nR = R_hi - R_lo + 1;
    const std::uint64_t nq = (qA_hi > q_bot) ? (qA_hi - q_bot) / chunk_q + 1 : 1;
    const std::uint64_t DENSE_CAP = 4096;
    const bool  dense = nR <= (nq > DENSE_CAP ? nq : DENSE_CAP);
    std::uint64_t nA  = dense ? nR : nq;
    if (nA < 1) nA = 1;
    s.dense = dense;

    const std::uint64_t nB = (q_top > qA_hi) ? (q_top - qA_hi) / chunk_q + 1 : 0;
    const std::uint64_t nC = nA + nB;

    // The seed at index c, both ways round. c = 0 is the balanced corner.
    auto seed_at = [&](std::uint64_t c, std::int64_t& Rc, std::int64_t& Ic) {
        if (dense) {
            std::uint64_t Rv = (c >= nA) ? R_hi : R_lo + c;
            if (Rv > R_hi) Rv = R_hi;
            Rc = (std::int64_t)Rv;
            Ic = (std::int64_t)csf_isqrt128((__uint128_t)Rv * Rv - N);
        } else {
            std::uint64_t qv = (c >= nA) ? qA_hi : q_bot + c * chunk_q;
            if (qv > qA_hi) qv = qA_hi;
            const std::uint64_t pv = N / qv;
            Rc = (std::int64_t)((qv + pv) / 2);
            Ic = (std::int64_t)((qv - pv) / 2);
        }
    };
    auto seed_q = [&](std::uint64_t c) -> std::uint64_t {
        std::int64_t Rc, Ic; seed_at(c, Rc, Ic);
        return (std::uint64_t)(Rc + Ic);
    };

    const std::uint64_t PAD = 20;            // absorbs the class rounding
    const std::uint64_t last =
        (first_chunk + max_chunks < nC) ? first_chunk + max_chunks : nC;

    for (std::uint64_t c = first_chunk; c < last; ++c) {
        std::int64_t Rs, Is;
        std::uint64_t qa, q_up, q_dn;
        bool want_up = true, want_dn = true;
        if (c < nA) {
            seed_at(c, Rs, Is);
            qa = (std::uint64_t)(Rs + Is);

            // ENDPOINTS WALK INWARD ONLY, and neighbours meet at the
            // MIDPOINT between them.
            //
            // Below the first seed there is nothing at all, and that is a
            // theorem: every factor pair has R = (p+q)/2 >= sqrt(pq) =
            // sqrt(N) by AM-GM, so R = ceil(sqrt N) is the smallest real
            // part any pair can have, and q rises with R along the
            // boundary. Above the last seed is the tail, which the chunks
            // below own -- and under an RSA window there is no tail.
            //
            // Between seeds each direction runs the FULL span to its
            // neighbour, not half of it. The two walks are not two halves
            // of one sweep: ascending from seed c and descending from
            // seed c+1 trace DIFFERENT (p, q) trajectories over the same
            // q range, and a trajectory only lands on the true p if it is
            // allowed to run the whole way. Cutting each at the midpoint
            // looked like a clean 1.98x saving and lost 147 of 690 cases,
            // all at q/p of 3 to 5, well inside the interval.
            want_dn = (c > 0);
            want_up = (c + 1 < nA) || nA == 1;
            q_up = seed_q(c + 1) + PAD;
            q_dn = (c > 0) ? (seed_q(c - 1) > PAD ? seed_q(c - 1) - PAD : q_bot)
                           : q_bot;
        } else {
            const std::uint64_t k  = c - nA;
            qa = qA_hi + k * chunk_q;
            if (qa == 0 || qa > q_top) continue;
            const std::uint64_t pv = N / qa;
            Rs = (std::int64_t)((qa + pv) / 2);
            Is = (std::int64_t)((qa - pv) / 2);
            q_up = qa + chunk_q + PAD;
            q_dn = (qa > chunk_q + PAD) ? qa - chunk_q - PAD : q_bot;
            (void)k;
        }
        if (q_up > q_top) q_up = q_top;
        if (q_dn < q_bot) q_dn = q_bot;
        if (qa >= q_top) want_up = false;
        if (qa <= q_bot) want_dn = false;

        for (const auto& cl : cls) {
            // Sync into the class. Both +-10 moves preserve it, so the
            // class holds for the whole walk -- that is how the sieve
            // applies to a two-directional search. The rounding lifts q
            // by at most 18, which is what PAD covers.
            const std::int64_t Rc = Rs + ((cl.first  - (int)(Rs % 10) + 10) % 10);
            const std::int64_t Ic = Is + ((cl.second - (int)(Is % 10) + 10) % 10);
            if (want_up) {
                s.R.push_back(Rc); s.I.push_back(Ic);
                s.qlimit.push_back((std::int64_t)q_up); s.up.push_back(1);
            }
            if (want_dn) {
                s.R.push_back(Rc); s.I.push_back(Ic);
                s.qlimit.push_back((std::int64_t)q_dn); s.up.push_back(0);
            }
        }
        ++s.chunks;
    }
}


// =====================================================================
//  STAGE 3 -- the yellow path (135-degree bracket traversal)
// =====================================================================
//
//  The most parallel of the three, and the reason is structural: along the
//  135-degree line every quantity is a closed form in j --
//
//      R_j = r + j     i_j = m - j     P_j = 2j + 3     R_j + i_j = S
//      V_j = R_j^2 - N                 H_j = i_j^2 + N
//
//  -- so a bracket can begin at ANY j for O(1) cost. There is no ladder to
//  preserve, unlike CTM. The engine's "additions only" advance is an
//  optimisation inside a bracket, not a constraint between them.
//
//  What the device actually needs is smaller than it looks. The two Fermat
//  legs are decided by a 7-modulus residue sieve on R and i; only survivors
//  need V or H evaluated, and those are 128-bit squares best left to GMP on
//  the host. So the kernel never forms V or H at all -- it carries 14 small
//  residues, advances them by +1 / -1, and emits the j values that survive.
//  That is pure table lookup and integer increment: no division, no
//  divergence, no wide arithmetic.
//
//  The trial-division leg rides along for free, and that is the point of the
//  bracket rather than a cost to be apologised for. P_j = 2j + 3 is a machine
//  word -- it runs only to N/S ~ sqrt(N)/2 -- so the bracket start hands the
//  leg its seed with one shift, and it advances by 2 alongside the traversal.
//  Two things follow:
//
//    * the divisor is one word, so the host engine's divisible_by_u64 reaches
//      mpz_divisible_ui_p and divides by a PRECOMPUTED reciprocal. Treating
//      P as a bignum instead, which this engine did until measured, costs a
//      general division every step: 28.77 ns against 5.90 ns at 96 bits,
//      44.01 ns against 10.03 ns at 112.
//    * 3 | P and 5 | P can be carried as two small counters and skipped
//      outright, which is sound whenever N is coprime to 3 and 5 (then no
//      such P can divide N) and drops 53% of the candidates: 3.44 ns and
//      3.52 ns at those two sizes. skip35 below is that guard, and it is a
//      guard rather than an assumption -- when N does share those factors
//      the skip is switched off, since then those P must still be tested.
//
//  What remains genuinely expensive is only the wide case: above 2^64 a
//  64-bit remainder no longer suffices and csf_mod128 falls back to a
//  64-iteration shift-subtract. do_trial is kept so that leg can still be
//  handed to the td stream or a segmented prime sieve there. Note the
//  completeness argument needs both legs (trial division catches
//  p <= 2*j_max+3, the vertical leg catches the rest), so switching it off
//  means something else must cover it.
//
//  Sizing: the engine cuts the path at CHUNK = 1<<8 because a bracket start
//  costs 211 ns against an 87 ns step in mpz -- a ratio of 2.4 that is flat
//  in N -- so ~120 steps per bracket holds the start under 2% at any size.
//  That gives ~1e6 brackets at 60 bits and ~1e9 at 96. Device-side the start
//  is cheaper still (14 modulos, no mpz), so 1<<8 is a safe floor, not a
//  tight one.

// ---------------------------------------------------------------------
// (hi:lo) mod P. Requires P < 2^63, which holds because P <= 2*j_max+3 and
// j_max ~ sqrt(N)/4, so P < 2^62 for any N this engine's native path takes.
//
// hi == 0 is the common case below 2^64 and costs one remainder. Above
// that it is a shift-subtract over the 64 bits of lo: r stays < P < 2^63,
// so r<<1 cannot overflow and a single conditional subtract re-reduces.
// ---------------------------------------------------------------------
CSF_HD inline std::uint64_t csf_mod128(std::uint64_t hi, std::uint64_t lo,
                                       std::uint64_t P)
{
    if (hi == 0) return lo % P;
    std::uint64_t r = hi % P;
    for (int b = 63; b >= 0; --b) {
        r = (r << 1) | ((lo >> b) & 1ull);
        if (r >= P) r -= P;
    }
    return r;
}

// ---------------------------------------------------------------------
// The per-j predicate: which of the three legs is still live at this j.
// Returns a bitmask -- 1 trial division, 2 vertical Fermat, 4 horizontal
// Fermat. A zero return means this j is dead and needs no host work.
//
// rm[t] / im[t] are R_j and i_j reduced mod CSF_YP_MODS[t]; the caller
// advances them by +1 / -1, which is what makes the walk division-free.
// ok_flat holds okR then okI, concatenated, with ok_off giving the start
// of each of the 14 tables.
// ---------------------------------------------------------------------
#define CSF_YP_NMOD 7
CSF_HD inline unsigned csf_yp_one(
    const std::uint8_t* __restrict__ rm,
    const std::uint8_t* __restrict__ im,
    const std::uint32_t* __restrict__ ok_off,   // 2 * CSF_YP_NMOD
    const std::uint8_t*  __restrict__ ok_flat,
    std::uint64_t P, unsigned p3, unsigned p5,  // P and P mod 3, P mod 5
    std::uint64_t nhi, std::uint64_t nlo,
    int do_trial, int skip35)
{
    unsigned out = 0;

    // p3 / p5 are carried by the caller, so the skip costs two compares and
    // removes 53% of the modulos -- see skip35 in the header.
    const int worth = !(skip35 && (p3 == 0 || p5 == 0));
    if (do_trial && worth && P > 1 && csf_mod128(nhi, nlo, P) == 0) out |= 1u;

    unsigned okv = 1u, okh = 1u;
    for (int t = 0; t < CSF_YP_NMOD; ++t) {
        okv &= ok_flat[ok_off[t] + rm[t]];
        okh &= ok_flat[ok_off[CSF_YP_NMOD + t] + im[t]];
    }
    out |= (okv << 1) | (okh << 2);
    return out;
}

// ---------------------------------------------------------------------
// Kernel: one thread per bracket. Thread t walks j in
// [j0 + t*chunk, j0 + (t+1)*chunk), seeding its 14 residues with one
// modulo each and then advancing them by +1 / -1.
//
// Survivors are packed as (j << 3) | mask and compacted through one
// atomic. The host re-derives R = r + j, i = m - j and finishes with GMP:
// the perfect-square tests on V and H, and the gcd.
// ---------------------------------------------------------------------
__global__ void csf_yp_kernel(
    std::uint64_t r, std::uint64_t m,           // r = ceil(sqrt N), m = r-3
    std::uint64_t j0, std::uint64_t chunk, std::uint64_t count,
    std::uint64_t jmax,
    const std::uint32_t* __restrict__ mods,     // CSF_YP_NMOD
    const std::uint32_t* __restrict__ ok_off,   // 2 * CSF_YP_NMOD
    const std::uint8_t*  __restrict__ ok_flat,
    std::uint64_t nhi, std::uint64_t nlo,
    int do_trial, int skip35,
    std::uint64_t* __restrict__ out,
    unsigned int*  __restrict__ out_n,
    unsigned int   out_cap)
{
#if defined(USE_CUDA) && !defined(CSF_CUDA_CPU_TEST)
    const std::uint64_t t = blockIdx.x * (std::uint64_t)blockDim.x + threadIdx.x;
    if (t >= count) return;

    const std::uint64_t a = j0 + t * chunk;
    if (a > jmax) return;
    std::uint64_t b = a + chunk; if (b > jmax + 1) b = jmax + 1;

    std::uint64_t R = r + a, I = m - a;
    std::uint8_t rm[CSF_YP_NMOD], im[CSF_YP_NMOD];
    for (int s = 0; s < CSF_YP_NMOD; ++s) {
        rm[s] = (std::uint8_t)(R % mods[s]);
        im[s] = (std::uint8_t)(I % mods[s]);
    }
    unsigned p3 = (unsigned)((2 * a + 3) % 3), p5 = (unsigned)((2 * a + 3) % 5);

    for (std::uint64_t j = a; j < b; ++j) {
        const unsigned msk = csf_yp_one(rm, im, ok_off, ok_flat,
                                        2 * j + 3, p3, p5, nhi, nlo,
                                        do_trial, skip35);
        if (msk) {
            const unsigned int slot = atomicAdd(out_n, 1u);
            if (slot < out_cap) out[slot] = (j << 3) | msk;
        }
        for (int s = 0; s < CSF_YP_NMOD; ++s) {
            const std::uint32_t p = mods[s];
            if (++rm[s] == p) rm[s] = 0;
            im[s] = (im[s] == 0) ? (std::uint8_t)(p - 1) : (std::uint8_t)(im[s] - 1);
        }
        p3 += 2; if (p3 >= 3) p3 -= 3;          // P advances by 2
        p5 += 2; if (p5 >= 5) p5 -= 5;
    }
#else
    (void)r; (void)m; (void)j0; (void)chunk; (void)count; (void)jmax;
    (void)mods; (void)ok_off; (void)ok_flat; (void)nhi; (void)nlo;
    (void)do_trial; (void)skip35; (void)out; (void)out_n; (void)out_cap;
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
// Mirrors the host side of a launch: build the paired (R mod 10, i mod 10)
// Fermat classes for this N's classification. csf_ctm_build_seeds below
// spreads seeds along the real interval and syncs each into these.
//
// Both the increases and the decreases are +-10 on R or on i, so a stream
// synced into a class at the start stays in it for the whole walk -- that
// is how the sieve applies to a two-directional search.
static void csf_ctm_classes(std::uint64_t N, std::vector<std::pair<int,int>>& out) {
    out.clear();
    for (int a = 0; a < 10; ++a)
        for (int b = 0; b < 10; ++b) {
            const long long v = (((long long)a*a - (long long)b*b) % 20 + 20) % 20;
            if (v == (long long)(N % 20)) out.push_back({a, b});
        }
}

static int test_ctm() {
    struct Case { std::uint64_t N, p, q; };
    const Case cases[] = {
        {8051ull,          83ull,       97ull},
        {143ull,           11ull,       13ull},
        {2923ull,          37ull,       79ull},
        {10963ull,         19ull,      577ull},
        {798607ull,       101ull,     7907ull},
        {1007509ull,      503ull,     2003ull},
        {288419ull,       379ull,      761ull},
        {5115191ull,     1597ull,     3203ull},
        {17821157ull,    3019ull,     5903ull},
        {10967535067ull, 104723ull, 104729ull},
    };

    int bad = 0;
    for (const Case& c : cases) {
        const std::uint64_t N = c.N;

        std::vector<std::pair<int,int>> cls;
        csf_ctm_classes(N, cls);

        // Walk the launch waves the way the GPU driver would: build a
        // bounded batch of seeds, run every lane, advance. Each case here
        // is small enough to finish in the first wave or two.
        CsfCtmSeeds seeds;
        std::uint64_t got = 0, first = 0, lanes = 0, waves = 0;
        for (; waves < 64 && !got; ++waves) {
            csf_ctm_build_seeds(N, cls, seeds, first, 256);
            if (seeds.chunks == 0) break;
            first += seeds.chunks;
            lanes += seeds.R.size();
            for (std::size_t t = 0; t < seeds.R.size(); ++t) {
                std::uint64_t hit = 0;
                if (csf_ctm_walk(seeds.R[t], seeds.I[t], seeds.qlimit[t],
                                 seeds.up[t], 0, N, &hit)) { got = hit; break; }
            }
        }

        const bool ok = (got == c.p || got == c.q);
        std::printf("[ctm]   N=%-13llu %6llu*%-8llu %-6s lanes=%-7llu waves=%-3llu"
                    " got=%-9llu %s\n",
                    (unsigned long long)N, (unsigned long long)c.p,
                    (unsigned long long)c.q, seeds.dense ? "dense" : "spaced",
                    (unsigned long long)lanes, (unsigned long long)waves,
                    (unsigned long long)got, ok ? "OK" : "FAIL");
        if (!ok) ++bad;
    }
    std::printf("[ctm]   %s\n", bad == 0 ? "WALK OK" : "WALK FAILED");
    return bad == 0 ? 0 : 1;
}

// --- csf_mod128 against a __uint128_t reference ---------------------
//
// Every semiprime in the tests above fits 64 bits, so they exercise only
// the hi == 0 fast path. The shift-subtract branch is what runs for a
// 128-bit N, and it is the only new arithmetic in this stage, so it gets
// its own check across the full width.
static int test_mod128() {
    std::mt19937_64 rng(11111);
    std::size_t bad = 0, tested = 0, wide = 0;
    for (int k = 0; k < 400000; ++k) {
        const std::uint64_t hi = (k % 8 == 0) ? 0 : rng();
        const std::uint64_t lo = rng();
        std::uint64_t P = rng() >> (1 + (rng() % 40));   // P < 2^63, varied
        P |= 1ull;
        if (P < 3) continue;
        const __uint128_t z = ((__uint128_t)hi << 64) | lo;
        const std::uint64_t want = (std::uint64_t)(z % (__uint128_t)P);
        const std::uint64_t got = csf_mod128(hi, lo, P);
        if (got != want) ++bad;
        if (hi) ++wide;
        ++tested;
    }
    // and the boundary cases the random draw will not hit
    const struct { std::uint64_t hi, lo, P; } edge[] = {
        {0, 0, 3}, {0, 1, 3}, {1, 0, 3}, {~0ull, ~0ull, 3},
        {~0ull, ~0ull, (1ull << 62) - 1}, {1, 0, (1ull << 62) - 1},
        {0, ~0ull, 5}, {12345, 0, 7}, {1, ~0ull, 0x7FFFFFFFFFFFFFFFull},
    };
    for (const auto& e : edge) {
        const __uint128_t z = ((__uint128_t)e.hi << 64) | e.lo;
        if (csf_mod128(e.hi, e.lo, e.P) !=
            (std::uint64_t)(z % (__uint128_t)e.P)) ++bad;
        ++tested;
    }
    std::printf("[mod128] cases %zu (%zu with hi != 0), mismatches %zu -> %s\n",
                tested, wide, bad, bad == 0 ? "OK" : "FAILED");
    return bad == 0 ? 0 : 1;
}

// --- stage 3: yellow path must not drop the true j ------------------
//
// Mirrors a launch: build the 14 residue tables, cut [0, j_max] into
// brackets, run the device walk on each, and collect survivors. Then check
// two things -- that every j the device keeps is one a full recomputation
// also keeps (no spurious divergence), and that the j carrying the factor
// is among them (nothing real was sieved away).
static std::vector<std::uint8_t> csf_sqmod(std::uint32_t m) {
    std::vector<std::uint8_t> s(m, 0);
    for (std::uint32_t x = 0; x < m; ++x) s[(std::uint64_t)x * x % m] = 1;
    return s;
}

static int test_yp() {
    const std::uint32_t MODS[CSF_YP_NMOD] = {64, 27, 25, 7, 11, 13, 17};
    struct Case { std::uint64_t N, p, q; };
    const Case cases[] = {
        {8051ull,             83ull,        97ull},
        {10967535067ull,  104723ull,    104729ull},
        {978508015703ull, 752867ull,   1299709ull},
        {314187ull,            3ull,    104729ull},
        {10963ull,            19ull,       577ull},
    };

    int bad = 0;
    for (const Case& c : cases) {
        const std::uint64_t N = c.N;
        std::uint64_t r = (std::uint64_t)std::sqrt((double)N);
        while (r * r < N) ++r;
        while (r > 1 && (r - 1) * (r - 1) >= N) --r;
        if (r < 4) continue;
        const std::uint64_t m = r - 3, S = r + m;
        const std::uint64_t jmax = (N / S - 3) / 2;

        // tables: okR then okI
        std::vector<std::uint32_t> mods(MODS, MODS + CSF_YP_NMOD);
        std::vector<std::uint32_t> ok_off(2 * CSF_YP_NMOD);
        std::vector<std::uint8_t> ok_flat;
        for (int pass = 0; pass < 2; ++pass)
            for (int t = 0; t < CSF_YP_NMOD; ++t) {
                const std::uint32_t mm = MODS[t];
                const std::vector<std::uint8_t> sq = csf_sqmod(mm);
                const std::uint64_t Np = N % mm, Nn = (mm - Np) % mm;
                ok_off[pass * CSF_YP_NMOD + t] = (std::uint32_t)ok_flat.size();
                for (std::uint32_t x = 0; x < mm; ++x) {
                    const std::uint64_t v =
                        (((std::uint64_t)x * x) % mm + mm - (pass ? Nn : Np)) % mm;
                    ok_flat.push_back(sq[v] ? 1 : 0);
                }
            }

        // skip35 is sound only when N is coprime to 3 and 5, exactly as the
        // engine decides it.
        const int skip35 = (N % 3 != 0 && N % 5 != 0) ? 1 : 0;

        // walk every bracket, exactly as the kernel would
        const std::uint64_t CH = 256;
        std::vector<std::uint64_t> got;
        std::uint64_t brackets = 0;
        for (std::uint64_t a = 0; a <= jmax; a += CH) {
            std::uint64_t b = a + CH; if (b > jmax + 1) b = jmax + 1;
            std::uint64_t R = r + a, I = m - a;
            std::uint8_t rm[CSF_YP_NMOD], im[CSF_YP_NMOD];
            for (int s = 0; s < CSF_YP_NMOD; ++s) {
                rm[s] = (std::uint8_t)(R % MODS[s]);
                im[s] = (std::uint8_t)(I % MODS[s]);
            }
            unsigned p3 = (unsigned)((2 * a + 3) % 3),
                     p5 = (unsigned)((2 * a + 3) % 5);
            ++brackets;
            for (std::uint64_t j = a; j < b; ++j) {
                const unsigned msk = csf_yp_one(rm, im, ok_off.data(),
                                                ok_flat.data(), 2 * j + 3,
                                                p3, p5, 0, N, 1, skip35);
                if (msk) got.push_back((j << 3) | msk);
                for (int s = 0; s < CSF_YP_NMOD; ++s) {
                    const std::uint32_t p = MODS[s];
                    if (++rm[s] == p) rm[s] = 0;
                    im[s] = (im[s] == 0) ? (std::uint8_t)(p - 1)
                                         : (std::uint8_t)(im[s] - 1);
                }
                p3 += 2; if (p3 >= 3) p3 -= 3;
                p5 += 2; if (p5 >= 5) p5 -= 5;
            }
        }

        // reference: recompute each j from scratch, no incremental state
        std::size_t mismatch = 0, k = 0;
        for (std::uint64_t j = 0; j <= jmax; ++j) {
            const std::uint64_t R = r + j, I = m - j, P = 2 * j + 3;
            unsigned want = 0;
            const int worth = !(skip35 && (P % 3 == 0 || P % 5 == 0));
            if (P > 1 && worth && N % P == 0) want |= 1u;
            unsigned okv = 1u, okh = 1u;
            for (int t = 0; t < CSF_YP_NMOD; ++t) {
                const std::uint32_t mm = MODS[t];
                okv &= ok_flat[ok_off[t] + R % mm];
                okh &= ok_flat[ok_off[CSF_YP_NMOD + t] + I % mm];
            }
            want |= (okv << 1) | (okh << 2);
            if (want) {
                if (k >= got.size() || got[k] != ((j << 3) | want)) ++mismatch;
                ++k;
            }
        }
        if (k != got.size()) ++mismatch;

        // the factor must still be reachable: p = 2j+3 on the trial leg,
        // or a = r + j with a^2 - N square on the vertical leg
        const std::uint64_t a_true = (c.p + c.q) / 2;
        const std::uint64_t j_vert = a_true >= r ? a_true - r : ~0ull;
        const std::uint64_t j_trial = c.p >= 3 ? (c.p - 3) / 2 : ~0ull;
        bool reachable = false;
        for (std::uint64_t v : got) {
            const std::uint64_t j = v >> 3;
            if ((v & 1u) && j == j_trial) reachable = true;
            if ((v & 2u) && j == j_vert)  reachable = true;
        }
        const bool ok = (mismatch == 0) && reachable;
        std::printf("[yp]    N=%-14llu j_max=%-9llu brackets=%-7llu kept=%-7zu"
                    " %s\n",
                    (unsigned long long)N, (unsigned long long)jmax,
                    (unsigned long long)brackets, got.size(),
                    ok ? "OK" : (mismatch ? "FAIL(divergence)"
                                          : "FAIL(factor lost)"));
        if (!ok) ++bad;
    }
    // The skip35 guard must be load-bearing, not decorative: on an N that
    // IS divisible by 3, forcing the skip on has to lose the factor. If this
    // check ever passes with the skip forced, the guard is doing nothing and
    // the real one above is untested.
    {
        const std::uint64_t N = 314187ull;          // 3 * 104729
        std::uint64_t r = (std::uint64_t)std::sqrt((double)N);
        while (r * r < N) ++r;
        while (r > 1 && (r - 1) * (r - 1) >= N) --r;
        const std::uint64_t m = r - 3, S = r + m, jmax = (N / S - 3) / 2;
        int found_with = 0, found_without = 0;
        for (int forced = 0; forced < 2; ++forced)
            for (std::uint64_t j = 0; j <= jmax; ++j) {
                const std::uint64_t P = 2 * j + 3;
                const int worth = !(forced && (P % 3 == 0 || P % 5 == 0));
                if (P > 1 && worth && N % P == 0) {
                    if (forced) ++found_with; else ++found_without;
                }
            }
        const bool ok = (found_without > 0) && (found_with == 0);
        std::printf("[yp]    skip35 guard: 3 | N finds %d hits with the skip "
                    "off, %d with it forced on -> %s\n",
                    found_without, found_with, ok ? "GUARD REQUIRED" : "FAIL");
        if (!ok) ++bad;
    }

    std::printf("[yp]    %s\n", bad == 0 ? "WALK OK" : "WALK FAILED");
    return bad == 0 ? 0 : 1;
}

int main() {
    const int a = test_sieve();
    const int b = test_ctm();
    std::printf("\n");
    const int c = test_mod128();
    const int d = test_yp();
    std::printf("\n%s\n", (a || b || c || d) ? "CUDA PREDICATES FAILED"
                                             : "CUDA PREDICATES OK");
    return (a || b || c || d) ? 1 : 0;
}
#endif
