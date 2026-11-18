#include <gmpxx.h>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// An exact modular Fermat sieve on the squared complex-plane Y = 2ab line.
// Every actual odd factor pair (p,q) with pq=N, gcd(N,2m)=1 induces
// H=(p^2+q^2)/2 and Y=(q^2-p^2)/2, and (Y,H) modulo m is in the wheel.
// This is a sieve for EXACT factor points, NOT a sieve for generic smooth
// nonzero CTM residuals.
using Z=mpz_class;
using Pair=std::pair<unsigned,unsigned>;
static unsigned mod_small(const Z& n,unsigned m){return static_cast<unsigned>(mpz_fdiv_ui(n.get_mpz_t(),m));}
static Z isqrt(const Z& x){Z s;mpz_sqrt(s.get_mpz_t(),x.get_mpz_t());return s;}
static Z ceil_sqrt(const Z& x){Z s=isqrt(x);if(s*s<x)++s;return s;}
static Z q_floor(const Z& N,const Z& Y){
    if(N<=0||Y<0)throw std::runtime_error("require N>0,Y>=0");
    return isqrt(Y+isqrt(N*N+Y*Y));
}
static Z q_ceil(const Z& N,const Z& Y){
    Z q=q_floor(N,Y);
    if(q*q*q*q-2*Y*q*q-N*N<0)++q;
    return q;
}
static Z count_ap(const Z& lo,const Z& hi,unsigned mod,unsigned r){
    if(hi<lo)return 0;
    Z aa,bb,top=hi-Z(r),bottom=lo-Z(1)-Z(r),d=mod;
    mpz_fdiv_q(aa.get_mpz_t(),top.get_mpz_t(),d.get_mpz_t());
    mpz_fdiv_q(bb.get_mpz_t(),bottom.get_mpz_t(),d.get_mpz_t());
    return aa-bb;
}
static std::vector<Pair> paired_wheel(const Z& N,unsigned m){
    if(m<2||m>20000)throw std::runtime_error("m must be in [2,20000]");
    if(mpz_gcd_ui(nullptr,N.get_mpz_t(),2ul*m)!=1)
        throw std::runtime_error("factor-2/prime-wheel precheck needed before wheel");
    const unsigned M=2*m;
    std::set<Pair> classes;
    const Z mod=M;
    const Z nmod=mod_small(N,M);
    for(unsigned r=1;r<M;r+=2){
        if(std::gcd(r,M)>1u)continue;
        Z inv,rv=r;
        if(mpz_invert(inv.get_mpz_t(),rv.get_mpz_t(),mod.get_mpz_t())==0)
            throw std::runtime_error("unexpected nonunit");
        Z s=(nmod*inv)%mod;
        // Exact integer numerators are even for the chosen odd lifts r,s.
        Z Y=(s*s-rv*rv)/2;
        Z H=(s*s+rv*rv)/2;
        classes.insert(Pair{mod_small(Y,m),mod_small(H,m)});
    }
    return {classes.begin(),classes.end()};
}
static void stats(const Z& N,unsigned m,const Z& A,const Z& B){
    if(N<=0||A<0||B<A)throw std::runtime_error("N>0 and 0<=A<=B needed");
    Z early=mpz_gcd_ui(nullptr,N.get_mpz_t(),2ul*m);
    if(early>1){
        std::cout<<"{\"status\":\""<<(early==N?"trivial_gcd":"small_factor_precheck")<<"\",\"factor\":\""<<early<<"\",\"modulus\":"<<m<<"}\n";
        return;
    }
    const auto C=paired_wheel(N,m);
    std::set<unsigned> yres;
    for(const auto& c:C)yres.insert(c.first);
    Z Y_count=0;
    for(const auto r:yres)Y_count+=count_ap(A,B,m,r);
    const Z qlo=std::max(ceil_sqrt(N),q_ceil(N,A));
    const Z qhi=std::min(isqrt(2*N),q_floor(N,B));
    Z q_candidate=0,q_units10=0;
    std::size_t qresidues=0;
    for(unsigned r=1;r<2*m;r+=2){
        if(std::gcd(r,2*m)>1u)continue;
        q_candidate+=count_ap(qlo,qhi,2*m,r);
        ++qresidues;
    }
    for(unsigned r:{1u,3u,7u,9u})q_units10+=count_ap(qlo,qhi,10,r);
    Z q_all=qhi>=qlo?qhi-qlo+1:Z(0);
    std::cout<<"{\"status\":\"ok\",\"N\":\""<<N<<"\",\"m\":"<<m
             <<",\"Ylo\":\""<<A<<"\",\"Yhi\":\""<<B
             <<"\",\"Ycount\":\""<<B-A+1
             <<"\",\"Ywheel\":\""<<Y_count
             <<"\",\"Yresidue_count\":"<<yres.size()
             <<",\"paired_residue_count\":"<<C.size()
             <<",\"qlo\":\""<<qlo<<"\",\"qhi\":\""<<qhi
             <<"\",\"qall\":\""<<q_all
             <<"\",\"qunit10\":\""<<q_units10
             <<"\",\"qunit2m\":\""<<q_candidate
             <<"\",\"qresidue_count\":"<<qresidues
             <<",\"full_y_period\":"<<m
             <<",\"note\":\"Y wheel is necessary only at exact factor points; every unit q residue remains admissible\"}\n";
}
int main(int argc,char**argv){
    try {
        if(argc<4)throw std::runtime_error("usage: y_fermat_wheel <N> <m> <Ylo> [Yhi] [--classes]");
        Z N(argv[1]);unsigned m=std::stoul(argv[2]);Z A(argv[3]);Z B=argc>=5&&argv[4][0]!='-'?Z(argv[4]):A;
        stats(N,m,A,B);
        if(argc>=6&&std::string(argv[5])=="--classes"&&mpz_gcd_ui(nullptr,N.get_mpz_t(),2ul*m)==1){
            const auto C=paired_wheel(N,m);
            std::cout<<"{\"pairs\":[";
            for(std::size_t i=0;i<C.size();++i){if(i)std::cout<<",";std::cout<<"["<<C[i].first<<","<<C[i].second<<"]";}
            std::cout<<"]}\n";
        }
        return 0;
    }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}
}
