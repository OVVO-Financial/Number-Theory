#include "gmp_regions.hpp"
#include <iostream>
#include <random>
#include <sstream>

namespace {
using vf::Z;
std::uint64_t checks=0,mapping_cases=0,ledger_cases=0,cursor_pairs=0,big_ranges=0;
void require(bool condition,const char* expression,int line) {
    ++checks;
    if (!condition) {
        std::ostringstream text;
        text << "check failed at line " << line << ": " << expression;
        throw std::runtime_error(text.str());
    }
}
#define REQUIRE(x) require(static_cast<bool>(x),#x,__LINE__)

void mapping_tests() {
    std::mt19937_64 random(8172207);
    for (unsigned value=3;value<=1001;value+=2) {
        const Z n(value), root=vf::floor_sqrt(n), first=vf::ceil_sqrt(n), endpoint=(n+1)/2;
        const auto span=Z(endpoint-first+1).get_ui();
        for (unsigned sample=0;sample<12;++sample) {
            Z a=first+static_cast<unsigned long>(random()%span);
            Z b=first+static_cast<unsigned long>(random()%span);
            if (a>b) std::swap(a,b);
            const auto interval=vf::excluded_factor_interval(n,a,b,Z(1));
            const auto integer_interval=vf::excluded_integer_midpoint_interval(n,a,b,Z(1));
            for (Z d=1;d<=root;++d) {
                const bool expected=2*a*d<=d*d+n && d*d+n<=2*b*d;
                const bool actual=interval && interval->lo<=d && d<=interval->hi;
                REQUIRE(expected==actual);
                if (n%d==0) {
                    const Z midpoint=(d+n/d)/2;
                    const bool expected_integer=a<=midpoint && midpoint<=b;
                    const bool actual_integer=integer_interval && integer_interval->lo<=d && d<=integer_interval->hi;
                    REQUIRE(expected_integer==actual_integer);
                }
            }
            if (a<b) {
                const Z split=(a+b)/2;
                const auto left=vf::excluded_integer_midpoint_interval(n,a,split,Z(1));
                const auto right=vf::excluded_integer_midpoint_interval(n,Z(split+1),b,Z(1));
                vf::CoverageLedger union_of_parts(Z(1),root);
                if (left) union_of_parts.add(*left);
                if (right) union_of_parts.add(*right);
                if (integer_interval) {
                    REQUIRE(union_of_parts.covers(integer_interval->lo,integer_interval->hi));
                    REQUIRE(union_of_parts.covered_length()==integer_interval->size());
                } else REQUIRE(union_of_parts.covered_length()==0);
            }
            ++mapping_cases;
        }
        for (unsigned sample=0;sample<6;++sample) {
            const auto high=root.get_ui();
            Z lo=1+static_cast<unsigned long>(random()%high);
            Z hi=1+static_cast<unsigned long>(random()%high);
            if (lo>hi) std::swap(lo,hi);
            const auto interval=vf::midpoint_interval_for_factors(n,lo,hi);
            for (Z r=first;r<=endpoint;++r) {
                // Compare the continuous inverse inequalities without reusing
                // roots or integer endpoint formulas from the implementation.
                const bool expected=hi*hi+n<=2*hi*r && 2*lo*r<=lo*lo+n;
                const bool actual=interval && interval->lo<=r && r<=interval->hi;
                REQUIRE(expected==actual);
            }
            ++mapping_cases;
        }
    }
    REQUIRE(!vf::excluded_factor_interval(Z(91),Z(100),Z(10)));
    REQUIRE(!vf::excluded_factor_interval(Z(91),Z(0),Z(9)));
    REQUIRE(!vf::midpoint_interval_for_factors(Z(91),Z(10),Z(8)));
    bool caught=false;
    try { vf::excluded_factor_interval(Z(90),Z(10),Z(20)); }
    catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
    for (unsigned p=3;p<=99;p+=2) for (unsigned q=p;q<=p+80;q+=2) {
        const Z n=Z(p)*q, midpoint=(p+q)/2;
        const auto at=vf::excluded_integer_midpoint_interval(n,midpoint,midpoint,Z(1));
        const auto before=vf::excluded_integer_midpoint_interval(n,Z(midpoint-20),Z(midpoint-1),Z(1));
        const auto after=vf::excluded_integer_midpoint_interval(n,Z(midpoint+1),Z(midpoint+20),Z(1));
        REQUIRE(at && at->lo<=p && p<=at->hi);
        REQUIRE(!before || !(before->lo<=p && p<=before->hi));
        REQUIRE(!after || !(after->lo<=p && p<=after->hi));
        // A prefix certificate recovers the existing global upper bound.
        const Z end=midpoint+3;
        const auto prefix=vf::excluded_integer_midpoint_interval(n,vf::ceil_sqrt(n),end,Z(1));
        REQUIRE(prefix && prefix->lo==vf::factor_upper_bound(n,end)+1);
    }
    caught=false;
    try { vf::midpoint_interval_for_factors(Z(90),Z(2),Z(8)); }
    catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
}

void ledger_tests() {
    vf::CoverageLedger ledger(Z(7),Z(100));
    REQUIRE(ledger.add(Z(40),Z(60)));
    REQUIRE(ledger.first_uncovered()==std::optional<Z>(Z(7)));
    REQUIRE(ledger.last_uncovered()==std::optional<Z>(Z(100)));
    REQUIRE(ledger.add(Z(7),Z(20)));
    REQUIRE(ledger.first_uncovered()==std::optional<Z>(Z(21)));
    REQUIRE(ledger.add(Z(80),Z(100)));
    REQUIRE(ledger.last_uncovered()==std::optional<Z>(Z(79)));
    REQUIRE(ledger.covered_length()==56);
    REQUIRE(ledger.add(Z(21),Z(39)));
    REQUIRE(ledger.first_uncovered()==std::optional<Z>(Z(61)));
    REQUIRE(ledger.intervals().size()==2);
    REQUIRE(!ledger.add(Z(10),Z(60)));
    REQUIRE(ledger.add(Z(61),Z(79)));
    REQUIRE(ledger.complete() && ledger.covered_length()==94);
    REQUIRE(!ledger.first_uncovered() && !ledger.last_uncovered());
    REQUIRE(ledger.gaps().empty());
    REQUIRE(!ledger.add(Z(-100),Z(200)));
    REQUIRE(!ledger.covers(Z(1),Z(100)));
    REQUIRE(ledger.covers(Z(7),Z(100)));

    std::mt19937 random(93179);
    for (unsigned trial=0;trial<60;++trial) {
        vf::CoverageLedger actual(Z(7),Z(170));
        std::vector<bool> covered(171,false);
        for (unsigned iteration=0;iteration<70;++iteration) {
            int lo=int(random()%220)-20,hi=lo+int(random()%18);
            bool changed=false;
            for (int x=std::max(7,lo);x<=std::min(170,hi);++x) {
                changed |= !covered[x]; covered[x]=true;
            }
            REQUIRE(actual.add(Z(lo),Z(hi))==changed);
            unsigned count=0;
            for (unsigned x=7;x<=170;++x) {
                REQUIRE(actual.covers(Z(x),Z(x))==covered[x]);
                count+=covered[x];
            }
            REQUIRE(actual.covered_length()==count);
            int first=7,last=170;
            while (first<=170 && covered[first]) ++first;
            while (last>=7 && covered[last]) --last;
            REQUIRE(actual.first_uncovered()==(first<=170?std::optional<Z>(Z(first)):std::nullopt));
            REQUIRE(actual.last_uncovered()==(last>=7?std::optional<Z>(Z(last)):std::nullopt));
            const int probe=int(random()%200)-10;
            int next=std::max(7,probe),previous=std::min(170,probe);
            while (next<=170 && covered[next]) ++next;
            while (previous>=7 && covered[previous]) --previous;
            REQUIRE(actual.next_uncovered(Z(probe))==(next<=170?std::optional<Z>(Z(next)):std::nullopt));
            REQUIRE(actual.previous_uncovered(Z(probe))==(previous>=7?std::optional<Z>(Z(previous)):std::nullopt));
            std::vector<bool> gap(171,false);
            for (const auto& interval:actual.gaps())
                for (unsigned x=interval.lo.get_ui();x<=interval.hi.get_ui();++x) gap[x]=true;
            for (unsigned x=7;x<=170;++x) REQUIRE(gap[x]!=covered[x]);
            ++ledger_cases;
        }
    }
    const Z huge=Z(1)<<256;
    vf::CoverageLedger large(huge,Z(huge+1000));
    REQUIRE(large.add(Z(huge+100),Z(huge+800)));
    REQUIRE(large.covered_length()==701);
    REQUIRE(large.first_uncovered()==std::optional<Z>(huge));
    REQUIRE(large.last_uncovered()==std::optional<Z>(Z(huge+1000)));
    REQUIRE(large.add(huge,Z(huge+99)));
    REQUIRE(large.first_uncovered()==std::optional<Z>(Z(huge+801)));
    vf::CoverageLedger empty(Z(9),Z(8));
    REQUIRE(empty.complete() && empty.gaps().empty() && empty.covered_length()==0);
}

void compare(const vf::RegionFermatCursor& a,const vf::RegionFermatCursor& b) {
    REQUIRE(a.start==b.start && a.end==b.end);
    REQUIRE(a.next_midpoint==b.next_midpoint);
    REQUIRE(a.covered_through==b.covered_through);
    REQUIRE(a.factor==b.factor && a.half_gap==b.half_gap);
    REQUIRE(a.candidates==b.candidates);
    REQUIRE(a.first_midpoint==b.first_midpoint && a.last_midpoint==b.last_midpoint);
    REQUIRE(a.completed==b.completed && a.status==b.status);
}

std::pair<Z,Z> direct_first_factor(const Z& n,const Z& start,const Z& end) {
    for (Z r=start;r<=end;++r) {
        const Z diff=r*r-n,root=vf::floor_sqrt(diff);
        if (root*root==diff && r-root>1 && r-root<n) return {Z(r-root),r};
    }
    return {Z(0),Z(0)};
}

void cursor_test(const Z& n,const Z& start,const Z& end,bool count_pair=true) {
    vf::RegionFermatCursor native(n,start,end,true),gmp(n,start,end,false);
    compare(native,gmp);
    REQUIRE(native.advance(0,-1)==gmp.advance(0,-1)); compare(native,gmp);
    REQUIRE(native.advance(100,0)==gmp.advance(100,0)); compare(native,gmp);
    const auto expected=direct_first_factor(n,gmp.start,gmp.end);
    while (!gmp.completed && gmp.factor==0) {
        const auto before=gmp.candidates;
        const auto delta=gmp.advance(3,-1);
        REQUIRE(native.advance(3,-1)==delta);
        REQUIRE(delta<=3 && gmp.candidates==before+delta);
        compare(native,gmp);
        const auto certified=gmp.certified_interval();
        if (certified) {
            REQUIRE(certified->lo==gmp.start && certified->hi<=gmp.end);
            const auto none=direct_first_factor(n,certified->lo,certified->hi);
            REQUIRE(none.first==0);
        }
        if (expected.first!=0) REQUIRE(gmp.covered_through<expected.second);
        if (expected.first!=0) {
            const auto excluded=gmp.excluded_factors(Z(1));
            REQUIRE(!excluded || !(excluded->lo<=expected.first && expected.first<=excluded->hi));
        }
    }
    REQUIRE(gmp.factor==expected.first);
    if (gmp.factor!=0) {
        REQUIRE(!gmp.completed && gmp.last_midpoint==expected.second);
        REQUIRE(n%gmp.factor==0);
    } else {
        REQUIRE(gmp.completed);
        if (gmp.start<=gmp.end) REQUIRE(gmp.covered_through==gmp.end);
    }
    REQUIRE(native.advance(100,-1)==0 && gmp.advance(100,-1)==0);
    compare(native,gmp);
    if (count_pair) ++cursor_pairs;
}

void cursor_tests() {
    std::mt19937_64 random(872631);
    for (unsigned value=3;value<=1001;value+=2) {
        const Z n(value),root=vf::ceil_sqrt(n),endpoint=(n+1)/2;
        const unsigned span=Z(endpoint-root+1).get_ui();
        for (unsigned sample=0;sample<4;++sample) {
            const Z begin=root+static_cast<unsigned long>(random()%span);
            const Z limit=std::min(endpoint,Z(begin+static_cast<unsigned long>(random()%50)));
            cursor_test(n,begin,limit);
        }
    }
    for (const char* text: {"18446744073709551615","18446744073709551557","18446743979220271189"}) {
        const Z n(text),root=vf::ceil_sqrt(n),endpoint=(n+1)/2;
        cursor_test(n,root,Z(root+300));
        cursor_test(n,Z(endpoint-300),endpoint);
        cursor_test(n,Z(endpoint+1),Z(endpoint+300));
    }
    // Seeking after an earlier factor must not claim the missing prefix.
    vf::RegionFermatCursor late(Z(77),Z(12),Z(15),true);
    late.advance(100,-1);
    const auto interval=late.certified_interval();
    REQUIRE(interval && interval->lo==12);
    const auto excluded=late.excluded_factors(Z(1));
    REQUIRE(!excluded || !(excluded->lo<=7 && 7<=excluded->hi));

    for (unsigned bits: {65U,80U,128U,256U,512U,1024U,2048U}) {
        for (unsigned offset=0;offset<4;++offset) {
            const Z midpoint=(Z(1)<<bits)+127+offset;
            const Z p=(Z(1)<<(bits/3))+13;
            const Z n=p*(2*midpoint-p);
            const Z begin=midpoint-20,end=midpoint+20;
            cursor_test(n,begin,end,false);
            vf::RegionFermatCursor cursor(n,begin,end,true);
            REQUIRE(cursor.backend=="gmp");
            REQUIRE(cursor.start>vf::ceil_sqrt(n));
            ++big_ranges;
        }
    }
    vf::RegionFermatCursor empty(Z(91),Z(100),Z(200),true);
    REQUIRE(empty.completed && !empty.certified_interval());
    REQUIRE(empty.advance(100,-1)==0);
    bool caught=false;
    try { vf::RegionFermatCursor bad(Z(90),Z(10),Z(100)); }
    catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
    caught=false;
    try { late.advance(1,std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { caught=true; }
    REQUIRE(caught);
}
}

int main() {
    try {
        const auto started=vf::Clock::now();
        mapping_tests(); ledger_tests(); cursor_tests();
        std::cout << "{\"status\":\"passed\",\"checks\":" << checks
                  << ",\"interval_mapping_cases\":" << mapping_cases
                  << ",\"ledger_updates\":" << ledger_cases
                  << ",\"native_gmp_range_pairs\":" << cursor_pairs
                  << ",\"large_seek_ranges\":" << big_ranges
                  << ",\"elapsed_seconds\":" << vf::seconds_since(started) << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
