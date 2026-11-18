// Certified Gaussian norm-one conditional factoring path; N-only execution.
// For z=a+bi with norm c=a*a+b*b a unit mod N, xi=z/conj(z)=z*z/c.
// Its imaginary numerator is exactly 2ab, the vertical green-line coordinate.
// Over a prime factor p, xi has multiplicative order dividing p-Legendre(-1,p).
// It is a special case of known Lucas/p+1 group-order factorization, not a new exponent bound.
// g++ -std=c++17 -O3 -Wall -Wextra gaussian_norm_one.cpp -lgmpxx -lgmp -o gaussian_norm_one
#include <gmpxx.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <stdexcept>
using Big=mpz_class;
static Big norm_mod(Big a,const Big&n){a%=n;if(a<0)a+=n;return a;}
static Big gcd_big(const Big&a,const Big&b){Big z;mpz_gcd(z.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return z;}
static bool proper(const Big&a,const Big&n){return a>1&&a<n;}
static std::vector<unsigned> primes(unsigned B){std::vector<unsigned char> s(B+1,1);s[0]=0;if(B>=1)s[1]=0;for(unsigned i=2;i<=B/i;++i)if(s[i])for(unsigned j=i*i;j<=B;j+=i)s[j]=0;std::vector<unsigned>v;for(unsigned i=2;i<=B;++i)if(s[i])v.push_back(i);return v;}
struct G {Big a,b;};
static G gmul(const G&u,const G&v,const Big&n){return {norm_mod(u.a*v.a-u.b*v.b,n),norm_mod(u.a*v.b+u.b*v.a,n)};}
static G gpow(G x,uint64_t e,const Big&n,uint64_t &cost){G r={1,0};for(;e;e>>=1){if(e&1){r=gmul(r,x,n);++cost;}if(e>1){x=gmul(x,x,n);++cost;}}return r;}
static uint64_t prime_power(unsigned p,unsigned B){uint64_t x=p;while(x<=B/p)x*=p;return x;}
struct GaussianResult{Big factor=0;std::string method="unresolved";uint64_t gauss_mults=0,gcd_calls=0;unsigned seed_a=0,seed_b=0;};
static GaussianResult solve(const Big&n,unsigned B,unsigned max_seeds,bool near_line){
 GaussianResult out;
 for(unsigned s: {2U,3U,5U,7U,11U})if(n%s==0){out.factor=s;out.method="small_factor";return out;}
 auto bases=primes(B);
 const std::vector<std::pair<unsigned,unsigned>> seeds={{2,1},{3,1},{3,2},{4,1},{4,3},{5,1},{5,2},{5,3},{5,4},{6,1},{7,2},{7,3}};
 Big root;mpz_sqrt(root.get_mpz_t(),n.get_mpz_t());if(root*root<n)++root;
 for(size_t si=0;si<std::min<size_t>(seeds.size(),max_seeds);++si){
  Big a,b;
  if(near_line){
   a=root+static_cast<unsigned long>(si+1);
   Big discr=a*a-n;mpz_sqrt(b.get_mpz_t(),discr.get_mpz_t());
   if(b==0)b=1;
   if((a-b)*(a+b)==n){out.factor=a-b;out.method="exact_strip_hit";return out;}
  }else{a=seeds[si].first;b=seeds[si].second;}
  out.seed_a=near_line?0:seeds[si].first;out.seed_b=near_line?0:seeds[si].second;
  Big c=a*a+b*b;
  Big found=gcd_big(c,n);++out.gcd_calls;
  if(proper(found,n)){out.factor=found;out.method="gaussian_norm_gcd";return out;}
  if(found!=1)continue;
  Big inv;mpz_invert(inv.get_mpz_t(),c.get_mpz_t(),n.get_mpz_t());
  G t={norm_mod((a*a-b*b)*inv,n),norm_mod(2*a*b*inv,n)};
  if(proper(found,n)){out.factor=found;out.method="gaussian_norm_gcd";return out;}
  if(found!=1)continue;
  for(unsigned j=0;j<bases.size();++j){
   t=gpow(t,prime_power(bases[j],B),n,out.gauss_mults);
   if((j%4)==3 || j+1==bases.size()){
    Big v=gcd_big(t.a-1,n);++out.gcd_calls;
    Big g=gcd_big(v,t.b);++out.gcd_calls;
    if(proper(g,n)){out.factor=g;out.method="gaussian_order";return out;}
   }
  }
 }
 return out;
}

extern "C" void dx_gaussian(mpz_srcptr N,mpz_ptr factor,unsigned B,unsigned seeds,bool near_line){
 Big n;mpz_set(n.get_mpz_t(),N);auto r=solve(n,B,seeds,near_line);mpz_set(factor,r.factor.get_mpz_t());
}

#ifndef DYADIC_LIB
int main(int argc,char**argv){
 if(argc<2){std::cerr<<"usage: gaussian_norm_one N [B] [seeds]\n";return 2;}
 try{
 Big n(argv[1]);if(n<3||n%2==0){std::cerr<<"odd N required\n";return 2;}
 unsigned B=argc>2?static_cast<unsigned>(std::stoul(argv[2])):512;
 unsigned ns=argc>3?static_cast<unsigned>(std::stoul(argv[3])):4;
 bool near_line=argc>4 ? std::string(argv[4])=="near" : false;
 auto st=std::chrono::steady_clock::now();auto out=solve(n,B,ns,near_line);
 auto us=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()-st).count();
 if(out.factor!=0&&(!proper(out.factor,n)||n%out.factor!=0))return 3;
 std::cout<<"{\"N\":\""<<n<<"\",\"factor\":\""<<out.factor<<"\",\"method\":\""<<out.method<<"\",\"us\":"<<us<<",\"gaussian_mults\":"<<out.gauss_mults<<",\"gcd_calls\":"<<out.gcd_calls<<",\"seed\":["<<out.seed_a<<","<<out.seed_b<<"]}\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}

#endif // DYADIC_LIB
