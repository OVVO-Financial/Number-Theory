#pragma once

// C ABI declarations for the unchanged, independently validated native kernel.
// Link vf_factor/fusion_kernel.cpp into the C++ executable.
#include <cstdint>

struct FusionSnapshot {
    std::uint64_t n;
    std::uint64_t factor;
    std::uint64_t candidates;
    std::uint64_t square_tests;
    std::uint64_t delta_candidates;
    std::uint64_t first_midpoint;
    std::uint64_t last_midpoint;
    std::uint64_t start_midpoint;
    std::uint64_t next_midpoint;
    std::uint64_t covered_through;
    std::uint64_t half_gap;
    std::uint64_t root_ceiling;
    std::uint64_t endpoint;
    std::uint32_t residue_count;
    std::uint32_t period;
    std::uint32_t status;
    std::uint32_t exhausted;
    double elapsed_seconds;
    double total_kernel_seconds;
};

extern "C" {
unsigned fusion_abi_version();
unsigned fusion_snapshot_size();
const char* fusion_compiler();
std::uint64_t fusion_isqrt(std::uint64_t n);
void* fusion_fermat_create(std::uint64_t n, unsigned* error);
void fusion_fermat_destroy(void* pointer);
int fusion_fermat_snapshot(void* pointer, FusionSnapshot* out);
int fusion_fermat_advance(void* pointer, std::uint64_t maximum,
                          double seconds, std::uint64_t midpoint_limit,
                          FusionSnapshot* out);
}
