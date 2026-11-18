#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

using Big = mpz_class;

Big fs(Big n, Big m, Big a, Big b) {
    if (m <= 0 || n < 0) throw std::runtime_error("invalid floor sum");
    Big aa, ab, bb, rb;
    mpz_fdiv_qr(aa.get_mpz_t(), ab.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t());
    mpz_fdiv_qr(bb.get_mpz_t(), rb.get_mpz_t(), b.get_mpz_t(), m.get_mpz_t());
    Big ans = aa * (n * (n-1) / 2) + bb*n;
    a=ab; b=rb;
    while (true) {
        if (a>=m) { Big d=a/m; ans += (n*(n-1)/2)*d; a%=m; }
        if (b>=m) { Big d=b/m; ans += n*d; b%=m; }
        Big ymax=a*n+b;
        if (ymax<m) break;
        n = ymax/m;
        b = ymax%m;
        std::swap(m,a);
    }
    return ans;
}

unsigned inverse10(unsigned r) {
    switch(r) { case 1:return 1; case 3:return 7; case 7:return 3; case 9:return 9; default:throw std::runtime_error("bad residue"); }
}

Big corridor_count_for_lane(const Big& N, const Big& Q, const Big& B,
                            unsigned rp, unsigned rq) {
    Big q0=Q+ ((rq+10-mpz_fdiv_ui(Q.get_mpz_t(),10))%10);
    if (q0>B) return 0;
    Big n=(B-q0)/10+1;
    Big H=B-Q, j0=q0-Q;
    const Big d1 = 10*Q*B;
    const Big d2 = 10*Q*Q;
    const Big a1 = -10*N;
    const Big b1 = N*(B-j0) - rp*Q*B;
    const Big a2 = 10*N;
    const Big b2 = rp*Q*Q - N*(Q-j0);
    return n + fs(n,d1,a1,b1) + fs(n,d2,a2,b2);
}

Big corridor_count(const Big& N,const Big& Q,const Big& B) {
    Big count=0;
    const unsigned n10=mpz_fdiv_ui(N.get_mpz_t(),10);
    for (unsigned rq : {1u,3u,7u,9u}) {
        const unsigned rp=(n10*inverse10(rq))%10;
        count += corridor_count_for_lane(N,Q,B,rp,rq);
    }
    return count;
}

struct Search {
    const Big& N;
    uint64_t nodes=0,pruned=0,positive=0,leaves=0,limit;
    explicit Search(const Big& N,uint64_t lim):N(N),limit(lim) {}
    bool find(const Big& L,const Big& U,Big& factor) {
        if (nodes>=limit) throw std::runtime_error("node limit reached");
        ++nodes;
        Big c=corridor_count(N,L,U);
        if (c==0) { ++pruned; return false; }
        ++positive;
        if (L==U) {
            ++leaves;
            if (L>1 && L<N && mpz_divisible_p(N.get_mpz_t(),L.get_mpz_t())) {
                factor=L; return true;
            }
            throw std::runtime_error("singleton corridor count inconsistent");
        }
        // The remainder of the search uses only exact monotonic geometry.
        Big mid=(L+U)/2;
        if (find(L,mid,factor)) return true;
        return find(mid+1,U,factor);
    }
};

int main(int argc,char** argv) {
    if(argc<2) {std::cerr<<"Usage: "<<argv[0]<<" N [max_nodes] [trial|corridor]\n";return 2;}
    try {
        Big N(argv[1]);
        if(N<=1) { std::cerr<<"N must be >1\n"; return 2; }
        uint64_t max_nodes=(argc>2)?std::stoull(argv[2]):10000000ull;
        std::string mode=(argc>3)?argv[3]:"corridor";
        if (mode=="count") {
            if(argc<6) throw std::runtime_error("count mode requires Q and B");
            Big Q(argv[4]), B(argv[5]);
            if(Q<=0 || B<Q || N%2==0 || N%5==0)
                throw std::runtime_error("count requires 0<Q<=B, odd N coprime to 5");
            std::cout<<"{\"count\":\""<<corridor_count(N,Q,B)<<"\"}"<<std::endl;
            return 0;
        }
        bool budget_exhausted=false;
        const auto t0=std::chrono::steady_clock::now();
        Big f=0;
        uint64_t nodes=0,pruned=0,positive=0;
        if(N==2 || N==5) { /* prime control */ }
        else if(N%2==0) f=2;
        else if(N%5==0) f=5;
        else {
            Big L,U,tmp=2*N;
            mpz_sqrt(L.get_mpz_t(),N.get_mpz_t());
            if(L*L<N) ++L;
            mpz_sqrt(U.get_mpz_t(),tmp.get_mpz_t());
            if(mode=="trial") {
                Big q=L;
                while(q<=U) {
                    unsigned r=mpz_fdiv_ui(q.get_mpz_t(),10);
                    if(r==1||r==3||r==7||r==9)break;
                    ++q;
                }
                while(q<=U && nodes<max_nodes) {
                    ++nodes;
                    if(q>1 && q<N && mpz_divisible_p(N.get_mpz_t(),q.get_mpz_t())) {
                        f=q;break;
                    }
                    unsigned r=mpz_fdiv_ui(q.get_mpz_t(),10);
                    q+=(r==3)?4:2; // 1->3->7->9->1, wheel 10
                }
                budget_exhausted=(f==0 && q<=U);
            } else if(mode=="corridor") {
                Search s(N,max_nodes);
                try{s.find(L,U,f);}catch(const std::runtime_error &e) {
                    if(std::string(e.what())!="node limit reached")throw;
                    budget_exhausted=true;
                }
                nodes=s.nodes;pruned=s.pruned;positive=s.positive;
            } else throw std::runtime_error("mode not recognized");
        }
        auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        const std::string status = (f>1 && f<N) ? "factor_found"
            : budget_exhausted ? "budget_exhausted"
            : (N==2 || N==5) ? "small_prime"
            : "no_factor_in_balanced_interval";
        std::cout << "{\"n\":\""<<N<<"\",\"method\":\""<<mode
                  <<"\",\"factor\":\""<<f<<"\",\"status\":\""<<status
                  <<"\",\"nodes\":"<<nodes<<",\"empty_regions\":"<<pruned
                  <<",\"positive_regions\":"<<positive<<",\"seconds\":"<<elapsed
                  <<",\"scope\":\"balanced q in [ceil_sqrt(N),floor_sqrt(2N)]\"}"
                  <<std::endl;
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl;return 2; }
}
