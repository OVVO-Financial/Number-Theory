// Exact dyadic iterated-average residual search, C++17/GMP.
// This is an arithmetic-progression / batch-GCD method; it does NOT
// establish a new factoring complexity bound or GNFS competitiveness.
#include <gmpxx.h>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Candidate {
    mpz_class value;
    std::uint64_t offset;
    char direction;
};

struct Counters {
    std::uint64_t candidates=0, products=0, gcd_calls=0, batches=0;
};

static bool is_factor(const mpz_class &d, const mpz_class &n) {
    return d > 1 && d < n && n % d == 0;
}

static bool inspect_tree(const std::vector<mpz_class> &tree,
                         const std::vector<Candidate> &values,
                         std::size_t index, std::size_t left,
                         std::size_t right, const mpz_class &n,
                         Counters &counter, mpz_class &answer,
                         std::uint64_t &offset, char &direction) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), tree[index].get_mpz_t(), n.get_mpz_t());
    ++counter.gcd_calls;
    if (g == 1) return false;
    if (right-left == 1) {
        if (left < values.size() && is_factor(g,n)) {
            answer=g; offset=values[left].offset; direction=values[left].direction;
            return true;
        }
        return false;
    }
    const std::size_t mid=(left+right)/2;
    return inspect_tree(tree,values,2*index,left,mid,n,counter,answer,offset,direction)
        || inspect_tree(tree,values,2*index+1,mid,right,n,counter,answer,offset,direction);
}

static bool inspect_batch(const std::vector<Candidate> &values,
                          const mpz_class &n, Counters &counter,
                          mpz_class &answer, std::uint64_t &offset,
                          char &direction) {
    if (values.empty()) return false;
    std::size_t base=1;
    while (base < values.size()) base*=2;
    std::vector<mpz_class> tree(base*2,1);
    for (std::size_t i=0;i<values.size();++i) {
        mpz_mod(tree[base+i].get_mpz_t(),values[i].value.get_mpz_t(),n.get_mpz_t());
    }
    for (std::size_t i=base-1;i>0;--i) {
        tree[i]=(tree[2*i]*tree[2*i+1])%n;
        ++counter.products;
    }
    ++counter.batches;
    return inspect_tree(tree,values,1,0,base,n,counter,answer,offset,direction);
}

int main(int argc, char **argv) {
    try {
        if (argc!=4 && argc!=5) {
            std::cerr << "Usage: dyadic_scan N level_k max_offset [batch_size=128]\n";
            return 2;
        }
        const mpz_class N(argv[1]);
        if (N <= 3 || N%2==0) throw std::runtime_error("N must be odd and > 3");
        const unsigned k=(unsigned)std::stoul(argv[2]);
        if (k<1 || k>256) throw std::runtime_error("level_k must be 1..256");
        const std::uint64_t max_offset=std::stoull(argv[3]);
        const std::size_t batch_size = argc==5 ? std::stoull(argv[4]):128;
        if (batch_size<1 || batch_size>65536) throw std::runtime_error("batch_size must be 1..65536");
        mpz_class r; mpz_sqrt(r.get_mpz_t(), N.get_mpz_t());
        if (r*r<N) ++r;
        const mpz_class w=mpz_class(1)<<k;
        const mpz_class C=(w-1)*N+r;
        mpz_class rounded, frac;
        mpz_fdiv_qr(rounded.get_mpz_t(),frac.get_mpz_t(),C.get_mpz_t(),w.get_mpz_t());
        const mpz_class origin=r-frac;
        Counters counter; mpz_class answer; std::uint64_t offset=0; char direction='?';
        bool success=false;
        for (std::uint64_t j=0;;) {
            std::vector<Candidate> candidates;
            candidates.reserve(2*batch_size);
            for (std::size_t b=0;b<batch_size && j<=max_offset;++b,++j) {
                const mpz_class d=origin-w*j;
                candidates.push_back({d,j,'-'});
                if (j!=0) candidates.push_back({origin+w*j,j,'+'});
                counter.candidates += (j==0 ? 1:2);
            }
            if (inspect_batch(candidates,N,counter,answer,offset,direction)) {
                success=true;break;
            }
            if (j>max_offset) break;
        }
        std::cout << "N=" << N << " k=" << k << " denominator=" << w
                  << " floor_average=" << rounded << " fractional_numerator=" << frac
                  << " batch_size=" << batch_size << " checked_terms=" << counter.candidates
                  << " product_multiplications=" << counter.products
                  << " gcd_calls=" << counter.gcd_calls
                  << " batches=" << counter.batches << '\n';
        if (success) {
            std::cout << "FACTOR " << answer << " * " << (N/answer)
                      << " offset=" << offset << " direction=" << direction << '\n';
            return 0;
        }
        std::cout << "UNRESOLVED (finite budget; not a primality result)\n";
        return 1;
    } catch(const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }
}
