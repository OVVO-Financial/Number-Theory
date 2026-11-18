#include "gmp_parallel.hpp"
#include <iostream>

namespace {
std::uint64_t checks=0;
void require(bool condition,const char* message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        vf::ParallelOptions options;
        options.threads=4; options.algorithm.timeout=1;
        {
            vf::Budget unlimited(-1);
            bool rejected=false;
            try { (void)vf::parallel_find_factor(vf::Z(77),options,unlimited,100); }
            catch (const std::invalid_argument&) { rejected=true; }
            require(rejected,"unlimited caller budget must be rejected");
        }
        const vf::Z hard=vf::Z("1000000007")*vf::Z("1000000009");
        options.algorithm.fermat_enabled=false; options.algorithm.vf_probes=false;
        options.algorithm.trace=true;
        for (const auto cap : {0ULL,1ULL,31ULL,64ULL,65ULL,200ULL,1000ULL,20000ULL}) {
            options.algorithm.max_steps=cap;
            const auto r=vf::factor_parallel(hard,options);
            require(r.status=="budget","bounded divisor search must remain unresolved");
            require(r.counters.get("work_units")<=cap,"global work exceeded cap");
            require(r.counters.get("reservation_high_water")<=cap,"reservations exceeded cap");
            require(r.outstanding_reservations==0,"unreleased reservations");
            require(r.counters.get("active_workers_high_water")<=4,"too many active workers");
            std::uint64_t traced=0, jobs=0;
            for (const auto& event : r.trace) {
                traced+=event.job_work;
                if (event.event=="job_result") ++jobs;
                if (event.certified_lo && event.certified_hi)
                    require(!(vf::Z("1000000007")>=*event.certified_lo &&
                              vf::Z("1000000007")<=*event.certified_hi),"true factor excluded");
            }
            require(traced==r.counters.get("work_units"),"trace lost completed work");
            require(jobs==r.counters.get("jobs_started"),"started job lost its result");
        }
        options.algorithm.max_steps.reset();
        options.factor_lo=vf::Z(2000); options.factor_hi=vf::Z(2200);
        const vf::Z regional=vf::Z(1009)*vf::Z("1000000007");
        const auto empty=vf::factor_parallel(regional,options);
        require(empty.status=="region_exhausted","partial band promoted to global primality");
        require(!empty.lo && !empty.hi,"finished band still has a coverage hole");
        options.factor_lo.reset(); options.factor_hi.reset();
        options.algorithm.fermat_enabled=true; options.algorithm.vf_probes=true;
        options.algorithm.initial_fermat_tests=0;
        for (const auto& text : {"978508015703","244360414137386395849947182580145646989"}) {
            const vf::Z n(text);
            const auto found=vf::factor_parallel(n,options);
            require(found.status=="factor","parallel factor not found");
            require(found.factor>1 && found.factor<n && mpz_divisible_p(n.get_mpz_t(),found.factor.get_mpz_t()),"invalid parallel factor");
            require(found.outstanding_reservations==0,"factor cancellation leaked reservation");
        }
        std::cout << "{\"status\":\"passed\",\"checks\":" << checks
                  << ",\"work_caps\":8,\"factor_cases\":2,\"unlimited_deadline_rejected\":true}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
