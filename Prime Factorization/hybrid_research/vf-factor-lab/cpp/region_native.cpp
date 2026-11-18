// Compile this translation unit INSTEAD OF vf_factor/fusion_kernel.cpp for
// regional searches. It retains the original validated kernel byte-for-byte
// and adds a fresh-cursor entry point starting at an arbitrary midpoint.
#include "../vf_factor/fusion_kernel.cpp"

extern "C" {

unsigned fusion_region_abi_version() { return 1; }

void* fusion_fermat_create_at(u64 n, u64 begin, unsigned* error) {
    if (error) *error = 0;
    if (n < 3 || !(n & 1U)) {
        if (error) *error = 1;
        return nullptr;
    }
    const auto began = Clock::now();
    try {
        auto* cursor = new Cursor(n);
        cursor->a = std::max(begin,cursor->root_ceiling);
        if (cursor->a > cursor->endpoint) {
            cursor->a = cursor->endpoint+1;
            cursor->covered = cursor->endpoint;
            cursor->start = cursor->a;
            cursor->difference = u128(cursor->a)*cursor->a-n;
            cursor->terminal = EXHAUSTED;
        } else {
            const unsigned remainder = cursor->a % 160;
            cursor->index = 0;
            while (cursor->index < cursor->count && cursor->residues[cursor->index] < remainder)
                ++cursor->index;
            if (cursor->index == cursor->count) {
                cursor->index = 0;
                cursor->a += 160-remainder+cursor->residues[0];
            } else cursor->a += cursor->residues[cursor->index]-remainder;
            cursor->start = cursor->a;
            cursor->covered = std::min(cursor->a-1,cursor->endpoint);
            cursor->difference = u128(cursor->a)*cursor->a-n;
            cursor->terminal = READY;
        }
        cursor->total_seconds = std::chrono::duration<double>(Clock::now()-began).count();
        return cursor;
    } catch (...) {
        if (error) *error = 2;
        return nullptr;
    }
}

} // extern C
