// =====================================================================
//  complex_space_factorization_cuda.cu
//
//  GPU sieve stage for complex_space_factorization.cpp.
//
//  WHY THIS SPLIT
//  --------------
//  Measured on a 96-bit semiprime with Fermat depth j = 285,049,336,514
//  (4 threads, --only=vf):
//
//      residues visited after the primary wheel : 72,971,567   (3,906x cut)
//      sqrt tests after the secondary moduli    :    242,744   (300x more)
//      bignum verification share of runtime     :      < 1%
//
//  So ~99% of the work is a fixed-width integer sieve over independent
//  residues, and the multi-precision part is a rounding error. That is the
//  classic host/device split: sieve on the GPU in machine words, verify
//  survivors with GMP on the host.
//
//  The CPU scan in complex_space_factorization.cpp was restructured to make
//  this port mechanical. Its inner loop is now:
//
//      * machine words only -- no mpz is built until a residue survives
//      * branchless -- all NS lookups folded with &, no warp divergence
//      * division-free -- sec_res holds residues[i] mod p, so the live
//        value is one add plus one conditional subtract
//      * coalesced -- sec_res is interleaved [i*NS + s]
//      * dependency-free -- residues are independent, so threads map 1:1
//
//  csf_sieve_one() below IS that loop body, shared verbatim between the CPU
//  engine and this kernel.
//
//  BUILD
//    nvcc -O3 -std=c++17 -DUSE_CUDA -c complex_space_factorization_cuda.cu
//
//  LOGIC TEST WITHOUT A GPU
//    g++ -std=c++17 -O2 -DCSF_CUDA_CPU_TEST -x c++ complex_space_factorization_cuda.cu -o t
//    ./t
//
//  STATUS: the predicate below is verified by the CPU test above and is the
//  same expression the engine's passing --selftest exercises. The CUDA
//  launch plumbing has NOT been executed -- this container has no GPU.
//  Treat the kernel wrapper as reviewed-but-unrun.
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
// out_idx via a single atomic; the host then builds a = base + residue,
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

// ---------------------------------------------------------------------
// CPU logic test: run the predicate over a synthetic wheel and check it
// against a direct recomputation. Exercises exactly the arithmetic the
// kernel performs, without needing a device.
// ---------------------------------------------------------------------
#ifdef CSF_CUDA_CPU_TEST
#include <cstdio>
#include <random>
#include <vector>

int main() {
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

    std::printf("residues tested : %zu\n", R);
    std::printf("survivors       : %zu  (%.4f)\n", pass, (double)pass / R);
    std::printf("mismatches      : %zu\n", mismatch);
    std::printf("%s\n", mismatch == 0 ? "KERNEL PREDICATE OK" : "KERNEL PREDICATE FAILED");
    return mismatch == 0 ? 0 : 1;
}
#endif
