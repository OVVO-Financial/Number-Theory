#include "gmp_product_worker.hpp"
#include <iostream>
#include <thread>

using vf::Z;
static std::uint64_t checks = 0, region_cases = 0, quota_cases = 0;

static void require(bool condition, const char* message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}

static bool wheel(unsigned d) {
    return d%2 && d%3 && d%5;
}

static Z root(const Z& n) {
    Z result;
    mpz_sqrt(result.get_mpz_t(),n.get_mpz_t());
    return result;
}

static bool reference_hit(unsigned n, unsigned lo, unsigned hi) {
    for (unsigned d=lo;d<=hi;++d)
        if (wheel(d) && mpz_gcd_ui(nullptr,Z(n).get_mpz_t(),d)>1) return true;
    return false;
}

static void verify_result(const Z& n, const Z& lo, const Z& hi,
                          const vf::ProductRegionResult& result,
                          std::uint64_t limit) {
    require(result.covered_lo==lo,"certificate lower endpoint");
    require(result.covered_hi>=lo-1 && result.covered_hi<=hi,"certificate upper endpoint");
    require(result.work_units<=limit,"work cap");
    require(result.work_units==result.counters.get("work_units"),"reported work agrees");
    require(result.work_units==result.counters.get("scalar_tests")+
                               result.counters.get("represented_candidates"),"exact work accounting");
    for (Z d=lo;d<=result.covered_hi;++d)
        require(!mpz_divisible_p(n.get_mpz_t(),d.get_mpz_t()),"certificate excludes a real divisor");
    if (result.factor!=0) {
        require(result.status=="factor","factor status");
        require(result.factor>1 && result.factor<n &&
                mpz_divisible_p(n.get_mpz_t(),result.factor.get_mpz_t()),"proper verified factor");
    } else require(result.status!="factor","factor status needs a factor");
    if (result.status=="completed") require(result.covered_hi==hi,"complete certificate");
}

template<class Function>
static void rejects(Function function, const char* message) {
    bool rejected=false;
    try { function(); } catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,message);
}

int main() {
    try {
        const auto started=vf::Clock::now();
        // Independent exhaustive short-region checks; reference_hit uses GCDs
        // only on wheel survivors, while certificates exclude actual divisors.
        for (unsigned n=49;n<=1201;n+=2) {
            if (!wheel(n)) continue;
            const unsigned upper=root(Z(n)).get_ui();
            if (upper<7) continue;
            vf::ProductRegionWorker worker{Z(n)};
            for (unsigned lo : {7U,std::max(7U,upper/2),upper}) {
                const unsigned hi=std::min(upper,lo+37);
                for (std::size_t degree : {1UL,7UL,16UL}) {
                    const auto r=worker.run(Z(lo),Z(hi),degree,{},-1,1000000);
                    verify_result(Z(n),Z(lo),Z(hi),r,1000000);
                    const bool hit=reference_hit(n,lo,hi);
                    require((r.factor!=0)==hit,"scalar witness reference");
                    require(r.status==(hit?"factor":"completed"),"unlimited region completes");
                    ++region_cases;
                }
            }
        }

        const Z nohit_n=Z(1009)*1013;
        vf::ProductRegionWorker limited(nohit_n);
        for (std::size_t degree : {1UL,2UL,7UL,8UL,16UL,32UL,64UL}) {
            for (std::uint64_t quota : {0ULL,1ULL,2ULL,7ULL,8ULL,9ULL,15ULL,16ULL,
                 31ULL,32ULL,63ULL,64ULL,65ULL,127ULL,128ULL,255ULL,256ULL,257ULL,512ULL}) {
                const auto r=limited.run(Z(7),Z(1008),degree,{},-1,quota);
                verify_result(nohit_n,Z(7),Z(1008),r,quota);
                require(r.factor==0,"no factor in quota-test region");
                if (quota) require(r.covered_hi>=7,"positive scalar quota makes progress");
                if (r.status=="budget") require(r.work_units==quota,"work budget fully consumed");
                else require(r.status=="completed","quota outcome");
                ++quota_cases;
            }
        }
        const auto gap=limited.run(Z(8),Z(10),16,{},-1,0);
        verify_result(nohit_n,Z(8),Z(10),gap,0);
        require(gap.status=="completed" && gap.work_units==0,"wheel-only gap needs no tests");

        // The cached engines refer to stable member addresses across fresh
        // job budgets and counters. Returned records must remain independent.
        const Z prime127=(Z(1)<<127)-1;
        for (bool shared : {true,false}) {
            vf::ProductRegionWorker cached(prime127,shared);
            std::vector<vf::Request> aux{{Z(5000),5,16,"vf"}};
            const auto first=cached.run(Z(120),Z(1079),16,aux,-1,10000);
            const auto second=cached.run(Z(1080),Z(2039),16,aux,-1,10000);
            verify_result(prime127,Z(120),Z(1079),first,10000);
            verify_result(prime127,Z(1080),Z(2039),second,10000);
            require(first.counters.get("polynomials_built")==unsigned(shared?1:2),"physical first polynomial builds");
            require(second.counters.get("polynomials_built")==0,"polynomials reused across jobs");
            require(second.counters.get("engine_cache_hits")>0,"cache-hit counter");
            require(first.counters.get("polynomials_built")==unsigned(shared?1:2),"previous result remains unchanged");
            require(second.counters.get("vf_candidates")==16,"auxiliary request is applied once");
        }

        // Aggregate GCD N from a trivial VF block must not hide a proper
        // factor in the later auxiliary request. The active strip is not
        // certified when a factor is returned.
        const Z mixed_n=Z(101)*103;
        vf::ProductRegionWorker mixed(mixed_n);
        std::vector<vf::Request> auxiliary{{mixed_n,5,1,"vf"},{Z(101),5,1,"vf"}};
        auto r=mixed.run(Z(7),Z(7),1,auxiliary,-1,3);
        verify_result(mixed_n,Z(7),Z(7),r,3);
        require(r.factor==101,"recovery continues past trivial VF hit");
        require(r.covered_hi==6,"factor-bearing batch does not certify its strip");
        require(r.counters.get("gcd_n_batches")==1,"GCD-N recovery was exercised");
        r=mixed.run(Z(7),Z(7),1,auxiliary,-1,2);
        verify_result(mixed_n,Z(7),Z(7),r,2);
        require(r.status=="completed" && r.factor==0,"unaffordable auxiliary probes skipped");
        require(r.counters.get("auxiliary_skipped_for_budget")==2,"skipped auxiliary count");
        r=mixed.run(Z(7),Z(7),1,{{Z(-101),-5,1,"vf"}},-1,2);
        require(r.factor==101,"signed auxiliary progression");

        // A genuine arbitrary-precision product job with a factor in its
        // assigned interval. All endpoints and candidates stay in GMP.
        const Z big_n=Z(65537)*((Z(1)<<257)-93);
        vf::ProductRegionWorker big(big_n);
        r=big.run(Z(60000),Z(70000),32,{},-1,10000);
        verify_result(big_n,Z(60000),Z(70000),r,10000);
        require(r.factor!=0,"large-modulus region finds a factor");

        std::atomic<bool> cancelled{true};
        r=limited.run(Z(7),Z(1008),16,{},1,10000,&cancelled);
        require(r.status=="cancelled" && r.covered_hi==6 && r.work_units==0,"initial cancellation");
        r=limited.run(Z(7),Z(1008),16,{},0,10000);
        require(r.status=="budget" && r.covered_hi==6 && r.work_units==0,"expired deadline certifies nothing");

        // Concurrent cancellation is polled between bounded polynomial stages.
        const Z prime521=(Z(1)<<521)-1;
        vf::ProductRegionWorker slow(prime521);
        cancelled.store(false,std::memory_order_relaxed);
        std::thread stopper([&]() {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            cancelled.store(true,std::memory_order_relaxed);
        });
        r=slow.run(Z(7),Z(500000),512,{},5,1000000,&cancelled);
        stopper.join();
        verify_result(prime521,Z(7),Z(500000),r,1000000);
        require(r.status=="cancelled" || r.status=="completed","cooperative concurrent cancellation");
        require(r.factor==0,"prime modulus returns no factor");

        rejects([](){vf::ProductRegionWorker bad(Z(15));},"reject nonunit wheel modulus");
        rejects([&](){limited.run(Z(6),Z(20),16,{},1,10);},"reject lower interval endpoint");
        rejects([&](){limited.run(Z(7),Z(2000),16,{},1,10);},"reject upper interval endpoint");
        rejects([&](){limited.run(Z(20),Z(7),16,{},1,10);},"reject reversed interval");
        rejects([&](){limited.run(Z(7),Z(20),0,{},1,10);},"reject zero degree");
        rejects([&](){limited.run(Z(7),Z(20),16,{},-.5,10);},"reject invalid negative timeout");
        rejects([&](){limited.run(Z(7),Z(20),16,{{Z(5),5,1,"vf"}},1,100);},"reject auxiliary degree mismatch");
        rejects([&](){mixed.run(Z(7),Z(20),1,{{Z(5),101,1,"vf"}},1,100);},"reject nonunit auxiliary step");

        std::cout << "{\"status\":\"passed\",\"checks\":" << checks
                  << ",\"region_cases\":" << region_cases << ",\"quota_cases\":" << quota_cases
                  << ",\"elapsed_seconds\":" << vf::seconds_since(started) << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "worker selftest failed: " << error.what() << '\n';
        return 1;
    }
}
