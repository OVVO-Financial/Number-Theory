#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include "gmp_products.hpp"

using Z = mpz_class;
static Z isqrt(const Z& n) {
    if (n < 0) throw std::invalid_argument("negative square root");
    Z out;
    mpz_sqrt(out.get_mpz_t(), n.get_mpz_t());
    return out;
}
static Z ceil_sqrt(const Z& n) {
    Z out = isqrt(n);
    if (out*out < n) ++out;
    return out;
}
// floor(sqrt(Y + sqrt(N*N+Y*Y))) is computed exactly by nested isqrt.
static Z q_floor_from_y(const Z& N, const Z& y) {
    if (N <= 0 || y < 0) throw std::invalid_argument("require N>0, Y>=0");
    return isqrt(y + isqrt(N*N + y*y));
}
static Z quartic_value(const Z& N, const Z& q, const Z& y) {
    const Z q2=q*q;
    return q2*q2-2*y*q2-N*N;
}
static Z q_ceil_from_y(const Z& N, const Z& y) {
    Z out=q_floor_from_y(N,y);
    if (quartic_value(N,out,y)<0) ++out;
    return out;
}
static Z q_to_y_floor(const Z& N,const Z& q) {
    if(q <= 0) throw std::invalid_argument("require q>0");
    const Z q2=q*q;
    Z res;
    const Z num=q2*q2-N*N;
    const Z den=2*q2;
    mpz_fdiv_q(res.get_mpz_t(),num.get_mpz_t(),den.get_mpz_t());
    return res;
}
static Z q_to_y_ceil(const Z& N,const Z& q) {
    if(q <= 0) throw std::invalid_argument("require q>0");
    const Z q2=q*q;
    Z res;
    const Z num=q2*q2-N*N;
    const Z den=2*q2;
    mpz_cdiv_q(res.get_mpz_t(),num.get_mpz_t(),den.get_mpz_t());
    return res;
}
static bool unit10(unsigned r) { return r==1 || r==3 || r==7 || r==9; }
static Z ap_count(const Z& L,const Z& U,unsigned modulus,unsigned residue) {
    if (U<L) return 0;
    Z r=residue, m=modulus, up=U-r,low=L-1-r;
    Z a,b;
    mpz_fdiv_q(a.get_mpz_t(),up.get_mpz_t(),m.get_mpz_t());
    mpz_fdiv_q(b.get_mpz_t(),low.get_mpz_t(),m.get_mpz_t());
    return a-b;
}
static Z unit10_count(const Z& L,const Z& U) {
    Z total=0;
    for(unsigned r:{1u,3u,7u,9u})total+=ap_count(L,U,10,r);
    return total;
}
static Z q_min(const Z& n) {return ceil_sqrt(n);}
static Z q_max_balanced(const Z& n) {return isqrt(2*n);}
static Z clip_max(const Z& a,const Z& b){return a<b?a:b;}
static Z clip_min(const Z& a,const Z& b){return a>b?a:b;}
static void print_sparse(const Z& n,const Z& ya,const Z& yb){
    if (ya<0 || yb<ya)throw std::invalid_argument("Y interval 0<=lo<=hi required");
    Z lo=q_ceil_from_y(n,ya), hi=q_floor_from_y(n,yb);
    lo=clip_min(lo,q_min(n));hi=clip_max(hi,q_max_balanced(n));
    Z admissible=unit10_count(lo,hi);
    Z all=hi>=lo?hi-lo+1:Z(0);
    std::cout<<"{\"n\":\""<<n<<"\",\"Ylo\":\""<<ya
             <<"\",\"Yhi\":\""<<yb<<"\",\"qlo\":\""<<lo
             <<"\",\"qhi\":\""<<hi<<"\",\"integer_q\":\""<<all
             <<"\",\"unit_q\":\""<<admissible
             <<"\",\"empty\":"<<(admissible==0?"true":"false")
             <<",\"scope\":\"balanced_only\"}\n";
}
struct Group {Z start;std::size_t count;};
struct BatchResult {
    Z factor=0;
    Z witness_lo=0;
    Z witness_hi=0;
    std::size_t inspected_blocks=0;
    bool budget_hit=false;
};
static BatchResult evaluate_batches(const Z& n,
                        const std::vector<Group>& groups,
                        std::size_t degree, std::size_t step) {
    BatchResult result;
    if(groups.empty())return result;
    vf::NormalizedProductEngine engine(n,degree);
    std::vector<Z> starts;
    std::vector<long> steps;
    starts.reserve(groups.size());steps.reserve(groups.size());
    for(const auto& g:groups){starts.push_back(g.start);steps.push_back(static_cast<long>(step));}
    std::vector<Z> values=engine.products(starts,steps);
    for(std::size_t i=0;i<groups.size();++i){
        ++result.inspected_blocks;
        Z factor;
        mpz_gcd(factor.get_mpz_t(),values[i].get_mpz_t(),n.get_mpz_t());
        if(factor==1)continue;
        result.witness_lo=groups[i].start;
        result.witness_hi=groups[i].start+Z(step)*(groups[i].count-1);
        if (factor>1 && factor<n){result.factor=factor;return result;}
        // Product gcd N can occur when different factors are present; recover.
        Z q=groups[i].start;
        for(std::size_t j=0;j<groups[i].count;++j,q+=Z(step)){
            mpz_gcd(factor.get_mpz_t(),q.get_mpz_t(),n.get_mpz_t());
            if(factor>1 && factor<n){result.factor=factor;return result;}
        }
    }
    return result;
}
static void emit_result(const Z& n, const std::string& method,
                        const std::string& status,const Z& f,
                        const Z& witness_lo,const Z& witness_hi,
                        const Z& qlo,const Z& qhi,
                        const Z& total_candidates,std::size_t blocks,
                        std::size_t degree,double elapsed,const Z& covered_qhi) {
    Z p=f>1?n/f:Z(0);
    Z q=f>1?f:Z(0);
    if(f>1 && p>q)std::swap(p,q);
    Z y=0,h=0;
    if(f>1){y=(q*q-p*p)/2;h=(q*q+p*p)/2;}
    std::cout<<"{\"N\":\""<<n<<"\",\"method\":\""<<method
             <<"\",\"status\":\""<<status<<"\",\"factor\":\""<<f
             <<"\",\"p\":\""<<p<<"\",\"q\":\""<<q
             <<"\",\"Y\":\""<<y<<"\",\"H\":\""<<h
             <<"\",\"witness_qlo\":\""<<witness_lo<<"\",\"witness_qhi\":\""<<witness_hi
             <<"\",\"candidate_qlo\":\""<<qlo<<"\",\"candidate_qhi\":\""<<qhi
             <<"\",\"covered_qhi\":\""<<covered_qhi
             <<"\",\"candidate_count\":\""<<total_candidates
             <<"\",\"batch_count\":"<<blocks<<",\"degree\":"<<degree
             <<",\"seconds\":"<<elapsed
             <<",\"scope\":\"balanced_candidate_q_interval\"}\n";
}
static void search(const Z& n, const std::string& method,
                   std::size_t max_degree,std::size_t max_blocks,
                   std::size_t specified_degree) {
    if(n<3)throw std::invalid_argument("N>=3 required");
    Z lo=q_min(n),hi=q_max_balanced(n), f=0, wlo=0, whi=0;
    const bool sieve840=(method=="batch840"||method=="wheel840");
    const unsigned modulus=sieve840?840u:10u;
    std::vector<unsigned> residues;
    for(unsigned r=1;r<modulus;++r)
        if(sieve840?std::gcd(r,modulus)==1u:unit10(r))residues.push_back(r);
    std::vector<std::uint8_t> admissible(modulus,0);
    Z total=0;for(unsigned r:residues){admissible[r]=1;total+=ap_count(lo,hi,modulus,r);}
    std::string status="no_factor_in_balanced_interval";
    const auto t0=std::chrono::steady_clock::now();
    std::size_t blocks=0,degree=0;
    Z covered=lo-1;
    if(n==2 || n==3 || n==5 || n==7){status="small_prime";}
    else if(n%2==0){f=2;status="factor_found";}
    else if(n%5==0){f=5;status="factor_found";}
    else if(sieve840 && n%3==0){f=3;status="factor_found";}
    else if(sieve840 && n%7==0){f=7;status="factor_found";}
    else if(method=="wheel"||method=="wheel840") {
        Z q=lo;
        while(q<=hi){
            unsigned r=mpz_fdiv_ui(q.get_mpz_t(),modulus);
            if(admissible[r]){
                ++blocks;
                if(blocks>max_blocks){status="budget_exhausted";break;}
                if(mpz_divisible_p(n.get_mpz_t(),q.get_mpz_t()) && q>1 && q<n){
                    f=q;wlo=q;whi=q;status="factor_found";break;
                }
            }
            covered=q;++q;
        }
    } else if(method=="batch"||method=="batch840"){
        if(specified_degree)degree=specified_degree;
        else {
            // Choose approximately sqrt(number of candidate values).
            if(!mpz_fits_ulong_p(total.get_mpz_t())) {
                status="budget_exhausted";
                emit_result(n,method,status,f,wlo,whi,lo,hi,total,blocks,degree,
                            std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count(),covered);
                return;
            }
            unsigned long t=mpz_get_ui(total.get_mpz_t());
            Z root=isqrt(Z(t));
            degree=std::max<std::size_t>(1,root.get_ui());
        }
        if(degree>max_degree||degree==0){status="budget_exhausted";}
        else {
            std::vector<Group> groups;
            std::vector<Z> tail;
            groups.reserve(std::min<std::size_t>(max_blocks,32768));
            for(unsigned r:residues) {
                Z q=lo+((r+modulus-mpz_fdiv_ui(lo.get_mpz_t(),modulus))%modulus);
                while(q<=hi){
                    Z available=(hi-q)/modulus+1;
                    if(available>=degree) {
                        if(groups.size()>=max_blocks){status="budget_exhausted";break;}
                        groups.push_back(Group{q,degree});
                        q+=Z(modulus)*degree;
                    } else {
                        while(q<=hi){tail.push_back(q);q+=modulus;}
                    }
                }
                if(status=="budget_exhausted")break;
            }
            if(status!="budget_exhausted"){
                BatchResult b=evaluate_batches(n,groups,degree,modulus);
                blocks=b.inspected_blocks;f=b.factor;wlo=b.witness_lo;whi=b.witness_hi;
                // A factor can stop GCD checks before all blocks and tails are
                // certified. Do NOT claim full interval coverage on early exit.
                if(f>1)status="factor_found";
                else for(const Z& q:tail){
                    Z g;
                    mpz_gcd(g.get_mpz_t(),q.get_mpz_t(),n.get_mpz_t());
                    if(g>1 && g<n){f=g;wlo=q;whi=q;status="factor_found";break;}
                }
                if(f==0)covered=hi;
            }
        }
    } else throw std::invalid_argument("method must be batch, wheel, batch840 or wheel840");
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    emit_result(n,method,status,f,wlo,whi,lo,hi,total,blocks,degree,elapsed,covered);
}
int main(int argc,char**argv) {
    try {
        if(argc<3)throw std::runtime_error("usage: global_vertical <sparse N Ylo Yhi | map N q | batch|batch840|wheel|wheel840 N [maxdegree maxblocks degree]>");
        const std::string mode=argv[1];
        const Z n(argv[2]);
        if(n<=0)throw std::runtime_error("N must be positive");
        if(mode=="sparse"){
            if(argc!=5)throw std::runtime_error("sparse N Ylo Yhi");
            print_sparse(n,Z(argv[3]),Z(argv[4]));
        }else if(mode=="map"){
            if(argc!=4)throw std::runtime_error("map N q");
            const Z q(argv[3]);
            if(q<=0)throw std::runtime_error("q must be positive");
            std::cout<<"{\"N\":\""<<n<<"\",\"q\":\""<<q
                <<"\",\"Yfloor\":\""<<q_to_y_floor(n,q)
                <<"\",\"Yceil\":\""<<q_to_y_ceil(n,q)
                <<"\",\"qback_floor\":\""<<q_floor_from_y(n,q_to_y_floor(n,q))<<"\"}\n";
        }else if(mode=="batch"||mode=="wheel"||mode=="batch840"||mode=="wheel840"){
            const std::size_t maxdegree=argc>3?std::stoull(argv[3]):2048;
            const std::size_t maxblocks=argc>4?std::stoull(argv[4]):20000;
            const std::size_t degree=argc>5?std::stoull(argv[5]):0;
            search(n,mode,maxdegree,maxblocks,degree);
        }else throw std::runtime_error("unknown action");
        return 0;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}
}
