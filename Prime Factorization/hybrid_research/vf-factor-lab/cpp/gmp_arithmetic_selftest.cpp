#include "gmp_arithmetic.hpp"

#include <iostream>
#include <random>
#include <sstream>

namespace {
std::uint64_t checks=0,root_cases=0,prime_cases=0,cursor_pairs=0,big_cursor_cases=0,bound_cases=0,big_resume_cases=0;
void require(bool condition,const char* expression,int line) {
    ++checks;
    if (!condition) {
        std::ostringstream message;
        message << "check failed at line " << line << ": " << expression;
        throw std::runtime_error(message.str());
    }
}
#define REQUIRE(x) require(static_cast<bool>(x),#x,__LINE__)

using vf::Z;

void compare(const vf::GmpFermatCursor& a,const vf::GmpFermatCursor& b) {
    REQUIRE(a.n==b.n);
    REQUIRE(a.next_midpoint==b.next_midpoint);
    REQUIRE(a.covered_through==b.covered_through);
    REQUIRE(a.factor==b.factor);
    REQUIRE(a.half_gap==b.half_gap);
    REQUIRE(a.candidates==b.candidates);
    REQUIRE(a.status==b.status);
    REQUIRE(a.exhausted==b.exhausted);
    REQUIRE(a.start_midpoint==b.start_midpoint);
    REQUIRE(a.first_midpoint==b.first_midpoint);
    REQUIRE(a.last_midpoint==b.last_midpoint);
    REQUIRE(a.root_ceiling==b.root_ceiling);
    REQUIRE(a.endpoint==b.endpoint);
    REQUIRE(a.residue_count==b.residue_count);
}

void root_tests() {
    gmp_randclass random(gmp_randinit_mt);
    random.seed(817629UL);
    for (unsigned bits : {1U,2U,31U,32U,63U,64U,65U,127U,128U,255U,256U,511U,512U,1024U}) {
        for (unsigned j=0;j<100;++j) {
            Z r=random.get_z_bits(bits);
            if (r==0) r=1;
            const Z sq=r*r;
            REQUIRE(vf::floor_sqrt(sq)==r);
            REQUIRE(vf::ceil_sqrt(sq)==r);
            REQUIRE(vf::floor_sqrt(Z(sq-1))==r-1);
            REQUIRE(vf::ceil_sqrt(Z(sq-1))==(r==1 ? Z(0):r));
            REQUIRE(vf::floor_sqrt(Z(sq+1))==r);
            REQUIRE(vf::ceil_sqrt(Z(sq+1))==r+1);
            root_cases+=3;
        }
    }
    REQUIRE(vf::floor_sqrt(Z(0))==0);
    REQUIRE(vf::ceil_sqrt(Z(0))==0);
    root_cases++;
    for (int a=-200;a<=200;++a) for (int b=1;b<=37;++b) {
        const Z q=vf::floor_div(Z(a),Z(b));
        REQUIRE(q*b<=a && a<(q+1)*b);
        const Z residue=vf::mod(Z(a),Z(b));
        REQUIRE(residue>=0 && residue<b);
        REQUIRE(q*b+residue==a);
    }
    for (std::uint64_t v : {std::uint64_t(0),std::uint64_t(1),std::uint64_t(1)<<63,
                           std::numeric_limits<std::uint64_t>::max()})
        REQUIRE(vf::z_to_u64(vf::z_from_u64(v))==v);
    bool caught=false;
    try { vf::floor_sqrt(Z(-1)); } catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
    caught=false;
    try { vf::z_to_u64(Z(1)<<64); } catch (const std::overflow_error&) { caught=true; }
    REQUIRE(caught);
}

void primality_tests() {
    const unsigned limit=200000;
    std::vector<bool> sieve(limit+1,true);
    sieve[0]=sieve[1]=false;
    for (unsigned p=2;p*p<=limit;++p) if (sieve[p])
        for (unsigned x=p*p;x<=limit;x+=p) sieve[x]=false;
    for (unsigned n=0;n<=limit;++n) {
        REQUIRE(vf::is_prime_u64(n)==sieve[n]);
        ++prime_cases;
    }
    for (const char* text : {"341550071728321","3825123056546413051",
                              "18446744073709551615","18446744030759878681"}) {
        REQUIRE(vf::prime_status(Z(text))=="composite"); ++prime_cases;
    }
    REQUIRE(vf::prime_status(Z("18446744073709551557"))=="prime"); ++prime_cases;
    for (unsigned exponent : {89U,107U,127U,521U}) {
        const Z prime=(Z(1)<<exponent)-1;
        REQUIRE(vf::prime_status(prime)=="probable_prime");
        REQUIRE(vf::prime_status(prime,false)=="unknown");
        REQUIRE(vf::prime_status(Z(prime*prime))=="composite");
        prime_cases+=3;
    }
}

void small_cursor_tests() {
    for (unsigned n=3;n<=2001;n+=2) {
        vf::GmpFermatCursor native(Z(n),true),gmp(Z(n),false);
        compare(native,gmp);
        REQUIRE(native.backend=="native64" && gmp.backend=="gmp");
        const Z initial=gmp.next_midpoint;
        REQUIRE(native.advance(0,-1,native.endpoint)==0);
        REQUIRE(gmp.advance(0,-1,gmp.endpoint)==0);
        compare(native,gmp);
        REQUIRE(native.advance(100,0,native.endpoint)==0);
        REQUIRE(gmp.advance(100,0,gmp.endpoint)==0);
        compare(native,gmp);
        REQUIRE(native.advance(100,-1,Z(initial-1))==0);
        REQUIRE(gmp.advance(100,-1,Z(initial-1))==0);
        compare(native,gmp);
        while (!gmp.exhausted && gmp.factor==0) {
            REQUIRE(native.advance(3,-1,native.endpoint)==gmp.advance(3,-1,gmp.endpoint));
            compare(native,gmp);
        }
        // Independent direct midpoint search with no wheel.
        Z expected=0;
        for (Z r=vf::ceil_sqrt(Z(n));r<=(n+1)/2;++r) {
            const Z diff=r*r-n;
            const Z b=vf::floor_sqrt(diff);
            if (b*b==diff && r-b>1 && r-b<n) { expected=r-b; break; }
        }
        REQUIRE(gmp.factor==expected);
        if (expected==0) {
            REQUIRE(gmp.exhausted);
            REQUIRE(gmp.covered_through==gmp.endpoint);
            REQUIRE(gmp.next_midpoint==gmp.endpoint+1);
        } else REQUIRE(gmp.covered_through<gmp.last_midpoint);
        REQUIRE(native.advance(100,-1,native.endpoint)==0);
        REQUIRE(gmp.advance(100,-1,gmp.endpoint)==0);
        compare(native,gmp);
        ++cursor_pairs;
    }
}

void random_cursor_tests() {
    std::mt19937_64 random(1730951);
    for (unsigned i=0;i<250;++i) {
        std::uint64_t value=random()|1;
        if (value<3) value=3;
        const Z n=vf::z_from_u64(value);
        vf::GmpFermatCursor native(n,true),gmp(n,false);
        compare(native,gmp);
        for (unsigned j=0;j<12;++j) {
            const Z cap=gmp.next_midpoint+static_cast<unsigned long>(random()%500);
            const auto maximum=random()%25;
            REQUIRE(native.advance(maximum,-1,cap)==gmp.advance(maximum,-1,cap));
            compare(native,gmp);
            if (gmp.exhausted || gmp.factor!=0) break;
        }
        ++cursor_pairs;
    }
    for (const char* text : {"18446743979220271189","18446744073709551557",
                            "18446744073709551615"}) {
        const Z n(text);
        vf::GmpFermatCursor native(n,true),gmp(n,false);
        REQUIRE(native.advance(4000,-1,native.endpoint)==gmp.advance(4000,-1,gmp.endpoint));
        compare(native,gmp);
        ++cursor_pairs;
    }
}

void big_cursor_tests() {
    for (unsigned bits : {65U,80U,96U,128U,256U,512U,1024U,2048U}) {
        const Z p=(Z(1)<<bits)+159;
        for (unsigned gap : {0U,2U,20U,2000U}) {
            const Z q=p+gap, n=p*q;
            vf::GmpFermatCursor cursor(n,true);
            REQUIRE(cursor.backend=="gmp");
            const Z original=cursor.next_midpoint;
            REQUIRE(cursor.advance(1,-1,Z(original-1))==0);
            REQUIRE(cursor.status=="range_limit");
            REQUIRE(!cursor.exhausted);
            REQUIRE(cursor.advance(2,-1,cursor.endpoint)==1);
            REQUIRE(cursor.factor==p);
            REQUIRE(n%cursor.factor==0);
            REQUIRE(cursor.covered_through<cursor.last_midpoint);
            REQUIRE(cursor.advance(100,-1,cursor.endpoint)==0);
            ++big_cursor_cases;
        }
    }
    bool caught=false;
    try { vf::GmpFermatCursor cursor(Z(4)); } catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
    vf::GmpFermatCursor cursor(Z(551),false);
    caught=false;
    try { cursor.advance(1,std::numeric_limits<double>::quiet_NaN(),Z(100)); }
    catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
    caught=false;
    try { cursor.advance(1,-1,Z(-1)); } catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
}

void big_resume_tests() {
    gmp_randclass random(gmp_randinit_mt);
    random.seed(738429UL);
    for (unsigned bits : {40U,64U,128U,256U,512U}) for (unsigned j=0;j<30;++j) {
        Z midpoint=random.get_z_bits(bits);
        mpz_setbit(midpoint.get_mpz_t(),bits-1);
        const Z d2=midpoint*static_cast<unsigned long>(60+4*j);
        Z half_gap=vf::floor_sqrt(d2);
        if (mpz_even_p(midpoint.get_mpz_t())==mpz_even_p(half_gap.get_mpz_t())) ++half_gap;
        const Z n=midpoint*midpoint-half_gap*half_gap;
        vf::GmpFermatCursor cursor(n,false);
        // Locate the first factor-bearing midpoint by direct integer arithmetic,
        // with no wheel and no incremental difference-of-squares recurrence.
        Z first=0,expected_factor=0;
        const Z start=vf::ceil_sqrt(n);
        for (Z r=start;r<=midpoint;++r) {
            const Z diff=r*r-n,b=vf::floor_sqrt(diff);
            if (b*b==diff && r-b>1 && r-b<n) { first=r;expected_factor=r-b;break; }
        }
        REQUIRE(first!=0);
        while (cursor.factor==0) {
            const Z next=cursor.next_midpoint,covered=cursor.covered_through;
            const auto before=cursor.candidates;
            REQUIRE(cursor.advance(1,-1,Z(next-1))==0);
            REQUIRE(!cursor.exhausted);
            REQUIRE(cursor.next_midpoint==next && cursor.covered_through==covered);
            const auto delta=cursor.advance(2,-1,midpoint);
            REQUIRE(delta>0 && delta<=2);
            REQUIRE(cursor.candidates==before+delta);
            REQUIRE(cursor.covered_through<first);
            if (cursor.factor==0) REQUIRE(cursor.next_midpoint>next);
        }
        REQUIRE(cursor.factor==expected_factor);
        REQUIRE(cursor.last_midpoint==first);
        REQUIRE(n%cursor.factor==0);
        ++big_resume_cases;
    }
}

void bound_tests() {
    gmp_randclass random(gmp_randinit_mt);
    random.seed(925037UL);
    for (unsigned bits : {8U,32U,64U,80U,128U,256U}) for (unsigned j=0;j<100;++j) {
        Z p=random.get_z_bits(bits),q=random.get_z_bits(bits+1);
        p|=1; q|=1;
        if (p<3) p=3;
        if (q<p) q+=p+1;
        if (mpz_even_p(q.get_mpz_t())) ++q;
        const Z n=p*q,midpoint=(p+q)/2;
        const Z root=vf::ceil_sqrt(n);
        const Z width=midpoint-root+1;
        const Z b=root+random.get_z_range(width);
        REQUIRE(vf::factor_upper_bound(n,Z(b-1))>=p);
        REQUIRE(vf::factor_upper_bound(n,midpoint)<p);
        const Z lo=1+random.get_z_range(p);
        REQUIRE(vf::fermat_midpoint_cap(n,lo)>=midpoint);
        REQUIRE(vf::fermat_midpoint_cap(n,p)==midpoint);
        ++bound_cases;
    }
}
}

int main() {
    try {
        root_tests(); primality_tests(); small_cursor_tests(); random_cursor_tests();
        big_cursor_tests(); big_resume_tests(); bound_tests();
        std::cout << "{\"status\":\"passed\",\"checks\":" << checks
                  << ",\"root_cases\":" << root_cases << ",\"primality_cases\":" << prime_cases
                  << ",\"native_gmp_cursor_pairs\":" << cursor_pairs
                  << ",\"big_cursor_cases\":" << big_cursor_cases
                  << ",\"big_resume_oracle_cases\":" << big_resume_cases
                  << ",\"bound_cases\":" << bound_cases << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
