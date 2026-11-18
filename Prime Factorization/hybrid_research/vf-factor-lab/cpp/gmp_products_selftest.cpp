#include "gmp_products.hpp"
#include <iostream>
#include <sstream>

using vf::Z;
using vf::gmp_poly::Poly;
static std::uint64_t checks = 0;

static void require(bool condition, const std::string& message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}

static Poly scalar_multiply(const Poly& a, const Poly& b, const Z& n) {
    if (a.empty() || b.empty()) return {};
    Poly out(a.size() + b.size() - 1);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) out[i + j] += a[i] * b[j];
    for (auto& z : out) z = vf::gmp_poly::residue(z, n);
    vf::gmp_poly::trim(out);
    return out;
}

static Poly scalar_remainder(Poly a, const Poly& b, const Z& n) {
    if (b.empty() || b.back() != 1) throw std::runtime_error("invalid reference divisor");
    while (a.size() >= b.size()) {
        const std::size_t offset = a.size() - b.size();
        const Z multiplier = a.back();
        for (std::size_t j = 0; j < b.size(); ++j)
            a[offset + j] = vf::gmp_poly::residue(a[offset + j] - multiplier * b[j], n);
        vf::gmp_poly::trim(a);
    }
    return a;
}

static Z scalar_evaluate(const Poly& a, const Z& x, const Z& n) {
    Z out = 0, power = 1;
    for (const auto& coefficient : a) {
        out = vf::gmp_poly::residue(out + coefficient * power, n);
        power = vf::gmp_poly::residue(power * x, n);
    }
    return out;
}

static Z scalar_product(const Z& start, long step, std::size_t count, const Z& n) {
    Z out = 1, candidate = start;
    for (std::size_t j = 0; j < count; ++j, candidate += step)
        out = vf::gmp_poly::residue(out * candidate, n);
    return out;
}

static Poly random_poly(gmp_randclass& random, const Z& n, std::size_t length, bool monic = false) {
    Poly out(length);
    for (auto& z : out) z = random.get_z_range(n);
    if (!out.empty() && monic) out.back() = 1;
    vf::gmp_poly::trim(out);
    return out;
}

static void same_polynomial(const Poly& a, const Poly& b, const std::string& message) {
    require(a.size() == b.size(), message + " length");
    for (std::size_t i = 0; i < a.size(); ++i) require(a[i] == b[i], message + " coefficient");
}

int main() {
    try {
        const auto started = vf::Clock::now();
        gmp_randclass random(gmp_randinit_mt);
        random.seed(1729117UL);
        Z big128 = ((Z(1) << 63) - 25) * ((Z(1) << 65) - 49);
        Z big256 = ((Z(1) << 127) - 1) * ((Z(1) << 129) - 25);
        Z big512 = ((Z(1) << 255) - 19) * ((Z(1) << 257) - 93);
        std::vector<Z> moduli{Z(2), Z(4), Z(6), Z(15), Z(35), Z(49), Z(77),
                              Z(1009)*1013, big128, big256, big512};
        vf::Counters stats;
        for (const auto& n : moduli) {
            vf::gmp_poly::Ring ring(n, nullptr, &stats);
            for (std::size_t size : {1UL, 4UL, 5UL, 8UL, 15UL, 16UL, 17UL, 31UL,
                                    32UL, 33UL, 63UL, 64UL, 65UL, 127UL, 128UL,
                                    129UL, 255UL, 256UL, 257UL, 511UL, 512UL}) {
                const Poly a = random_poly(random, n, size, true);
                const Poly b = random_poly(random, n, size / 2 + 1, true);
                same_polynomial(ring.multiply(a,b), scalar_multiply(a,b,n), "random multiplication");
                const Poly worst_a(size, n - 1), worst_b(size + 1, n - 1);
                same_polynomial(ring.multiply(worst_a,worst_b), scalar_multiply(worst_a,worst_b,n), "carry boundary");
            }
            same_polynomial(ring.multiply({}, {Z(1)}), {}, "empty multiplication");
            for (std::size_t degree : {1UL,2UL,4UL,8UL,9UL,15UL,16UL,17UL,31UL,32UL,33UL,64UL,129UL,257UL}) {
                const Poly divisor = random_poly(random,n,degree+1,true);
                for (std::size_t qlength : {1UL,2UL,4UL,5UL,17UL,33UL,65UL}) {
                    Poly dividend = random_poly(random,n,degree+qlength,true);
                    same_polynomial(ring.remainder(dividend,divisor), scalar_remainder(dividend,divisor,n), "monic remainder");
                }
                const Poly short_a = random_poly(random,n,degree);
                same_polynomial(ring.remainder(short_a,divisor),short_a,"short remainder");
                same_polynomial(ring.remainder(divisor,divisor),{},"equal remainder");
            }
            for (std::size_t length : {1UL,2UL,3UL,5UL,8UL,9UL,16UL,17UL,31UL,32UL,33UL,64UL,65UL,127UL,128UL,129UL,257UL,512UL}) {
                Poly f = random_poly(random,n,length,true); f[0]=1;
                const Poly inverse = ring.inverse_series_one(f,length);
                Poly product = scalar_multiply(f,inverse,n); product.resize(length);
                require(product[0]==1,"inverse constant");
                for (std::size_t i=1;i<length;++i) require(product[i]==0,"inverse coefficient");
            }
            for (std::size_t degree : {1UL,7UL,16UL,17UL,64UL,65UL,128UL,257UL,512UL}) {
                Poly polynomial=ring.progression_polynomial(Z(1),degree);
                for (std::size_t count : {0UL,1UL,2UL,4UL,5UL,7UL,8UL,9UL,17UL,33UL,65UL,129UL}) {
                    std::vector<Z> points;
                    for (std::size_t i=0;i<count;++i)
                        points.push_back(i%3==0 ? Z(0) : i%3==1 ? Z(n-1) : Z(random.get_z_range(n)));
                    auto values=ring.multipoint(polynomial,points);
                    require(values.size()==points.size(),"multipoint result count");
                    for (std::size_t i=0;i<count;++i) {
                        require(values[i]==scalar_evaluate(polynomial,points[i],n),"multipoint Horner reference");
                        require(values[i]==scalar_product(points[i],1,degree,n),"multipoint progression reference");
                    }
                }
            }
        }
        for (const auto& n : std::vector<Z>{Z(77),Z(1009)*1013,big128,big256,big512}) {
            Z unit_gcd;
            mpz_gcd_ui(unit_gcd.get_mpz_t(),n.get_mpz_t(),30);
            if (unit_gcd != 1) continue;
            for (std::size_t degree : {1UL,7UL,16UL,17UL,31UL,65UL,128UL,257UL,512UL}) {
                vf::Counters shared_stats, unshared_stats;
                vf::NormalizedProductEngine shared(n,degree,nullptr,&shared_stats,true);
                vf::NormalizedProductEngine separate(n,degree,nullptr,&unshared_stats,false);
                require(shared.products({},{}).empty(),"empty normalized products");
                require(shared_stats.get("polynomials_built")==0,"lazy construction");
                std::vector<Z> starts;
                std::vector<long> steps;
                for (std::size_t i=0;i<21;++i) {
                    Z start=random.get_z_range(n*3);
                    if (i%2) start=-start;
                    starts.push_back(start);
                    steps.push_back(i%3==0?5:i%3==1?30:-5);
                }
                const auto values=shared.products(starts,steps);
                const auto unshared=separate.products(starts,steps);
                const auto repeated=shared.products(starts,steps);
                require(values==repeated,"repeat values");
                for (std::size_t i=0;i<starts.size();++i) {
                    const Z expected=scalar_product(starts[i],steps[i],degree,n);
                    Z step=vf::gmp_poly::residue(Z(steps[i]),n), power;
                    mpz_powm_ui(power.get_mpz_t(),step.get_mpz_t(),static_cast<unsigned long>(degree),n.get_mpz_t());
                    require(vf::gmp_poly::residue(values[i]*power,n)==expected,"normalized value identity");
                    require(unshared[i]==expected,"separate value identity");
                }
                require(shared_stats.get("polynomials_built")==1,"shared polynomial count");
                require(unshared_stats.get("polynomials_built")==3,"separate polynomial count");
                require(shared_stats.get("represented_candidates")==2*starts.size()*degree,"represented count");
            }
        }
        auto recover = [&](const Z& n,std::vector<vf::Request> requests) {
            std::vector<Z> values;
            for (const auto& r:requests) values.push_back(scalar_product(r.start,r.step,r.count,n));
            return vf::recover_mixed_batch(n,requests,values,nullptr,&stats);
        };
        require(recover(Z(77),{{Z(77),5,1,"vf"},{Z(7),30,1,"divisor"}})==7,"skip trivial VF hit");
        require(recover(Z(77),{{Z(7),4,3,"divisor"}})==7,"recover inside mixed-factor block");
        require(recover(Z(77),{{Z(77),5,1,"vf"},{Z(1),30,1,"divisor"}})==0,"trivial-only batch");
        require(recover(Z(49),{{Z(49),5,1,"vf"},{Z(7),30,1,"divisor"}})==7,"prime-power recovery");
        require(recover(Z(77),{{Z(7),30,1,"divisor"},{Z(11),30,1,"divisor"}})==7,"aggregate gcd N across blocks");
        require(recover(Z(77),{{Z(-77),5,1,"vf"},{Z(-7),30,1,"divisor"}})==7,"negative recovery");
        require(recover(Z(77),{})==0,"empty recovery");
        bool rejected=false;
        try { vf::NormalizedProductEngine e(Z(35),16); e.products({Z(1)},{5}); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"nonunit step rejection");
        rejected=false;
        try { vf::NormalizedProductEngine e(Z(77),0); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"zero degree rejection");
        rejected=false;
        try { vf::NormalizedProductEngine e(Z(77),16); e.products({Z(1)},{}); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"mismatched batch rejection");
        bool expired=false;
        try { vf::Budget b(0); vf::NormalizedProductEngine e(Z(77),16,&b); e.products({Z(1)},{5}); }
        catch (const vf::BudgetExceeded&) { expired=true; }
        require(expired,"expired deadline");
        require(stats.get("kronecker_multiplications")>0,"accelerated multiplication used");
        require(stats.get("newton_remainders")>0,"accelerated division used");
        std::cout << "{\"status\":\"passed\",\"checks\":" << checks
                  << ",\"modulus_bits_max\":512,\"maximum_degree\":512,\"kronecker_multiplications\":"
                  << stats.get("kronecker_multiplications") << ",\"newton_remainders\":"
                  << stats.get("newton_remainders") << ",\"elapsed_seconds\":"
                  << vf::seconds_since(started) << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "self-test failed after " << checks << " checks: " << error.what() << '\n';
        return 1;
    }
}
