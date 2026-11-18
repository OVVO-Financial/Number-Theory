// N-only conditional smooth-order fast path: Pollard p-1, Lucas p+1, Brent rho control.
// Distinct from the dyadic residual sieve: exploits GROUP ORDERS, not just binary digits.
// Compile: g++ -std=c++17 -O3 -Wall -Wextra smooth_order.cpp -lgmpxx -lgmp -o smooth_order
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <limits>

using Big=mpz_class;
static Big mod_pos(Big x,const Big& N){x%=N;if(x<0)x+=N;return x;}
static Big gcd_big(const Big& a, const Big& b){Big d;mpz_gcd(d.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return d;}
static bool proper(const Big&d,const Big&N){return d>1&&d<N;}
static std::vector<unsigned> primes_up_to(unsigned B){
 std::vector<bool> prime(B+1,true);prime[0]=false;if(B>=1)prime[1]=false;
 for(unsigned j=2;j<=B/j;++j)if(prime[j])for(unsigned z=j*j;z<=B;z+=j)prime[z]=false;
 std::vector<unsigned> v;for(unsigned i=2;i<=B;++i)if(prime[i])v.push_back(i);return v;
}
static uint64_t highest_prime_power(unsigned p,unsigned B){uint64_t a=p;while(a<=B/p)a*=p;return a;}
static Big powm_ui(const Big& a,unsigned long e,const Big&N){Big out;mpz_powm_ui(out.get_mpz_t(),a.get_mpz_t(),e,N.get_mpz_t());return out;}
// Compute V_e(P,1) modulo N. V_0=2, V_1=P, V_(n+1)=P*V_n - V_(n-1).
static Big lucas_V(const Big &Pin, uint64_t e, const Big &N){
 Big P=mod_pos(Pin,N),v=Big(2)%N,v1=P;
 if(e==0)return v;
 int top=63;while(top>0 && ((e>>top)&1ULL)==0)--top;
 for(int b=top;b>=0;--b){
  Big v0=mod_pos(v*v-2,N);
  Big vo=mod_pos(v*v1-P,N);
  Big v2=mod_pos(v1*v1-2,N);
  if((e>>b)&1ULL){v=std::move(vo);v1=std::move(v2);}else{v=std::move(v0);v1=std::move(vo);}
 }
 return v;
}
struct SmoothResult{Big factor=0;std::string which="unresolved";uint64_t modular_steps=0,gcd_calls=0;};
static SmoothResult smooth_orders(const Big&N,unsigned B,bool do_minus,bool do_lucas){
 SmoothResult out;
 for(unsigned p: {2U,3U,5U,7U,11U}){if(N%p==0){out.factor=p;out.which="small_factor";return out;}}
 auto primes=primes_up_to(B);
 auto test=[&](const Big& value,const std::string&name){++out.gcd_calls;Big g=gcd_big(value,N);if(proper(g,N)){out.factor=g;out.which=name;return true;}return false;};
 if(do_minus){
  for(unsigned base: {2U,3U,5U}){
   Big t=base;
   ++out.gcd_calls;Big gd=gcd_big(t,N);if(proper(gd,N)){out.factor=gd;out.which="pminus_gcd_base";return out;}
   for(unsigned j=0;j<primes.size();++j){
    uint64_t pe=highest_prime_power(primes[j],B);
    t=powm_ui(t,static_cast<unsigned long>(pe),N);++out.modular_steps;
    if((j%8)==7 || j+1==primes.size())if(test(t-1,"pminus_stage1"))return out;
   }
  }
 }
 if(do_lucas){
  unsigned selected=0;
  for(unsigned P=3;P<65;++P){
   Big D=Big(P)*P-4;
   ++out.gcd_calls;Big gd=gcd_big(D,N);
   if(proper(gd,N)){out.factor=gd;out.which="lucas_discriminant_gcd";return out;}
   if(gd==1 && mpz_jacobi(D.get_mpz_t(),N.get_mpz_t())==-1){selected=P;break;}
  }
  // Jacobi -1 guarantees opposite quadratic characters across two prime factors.
  // If none in the bounded deterministic scan, use P=3 (unconditional correctness of found GCD).
  if(selected==0)selected=3;
  for(unsigned shift=0;shift<3;++shift){
   Big t=selected+2*shift;
   for(unsigned j=0;j<primes.size();++j){
    uint64_t pe=highest_prime_power(primes[j],B);
    t=lucas_V(t,pe,N);++out.modular_steps;
    if((j%8)==7 || j+1==primes.size())if(test(t-2,"lucas_stage1"))return out;
   }
  }
 }
 return out;
}
static SmoothResult brent_rho(const Big&N,uint64_t budget){
 SmoothResult out;
 for(unsigned p: {2U,3U,5U,7U,11U})if(N%p==0){out.factor=p;out.which="small_factor";return out;}
 Big y=2,c=1,x=0,ys=0,g=1,q=1;
 uint64_t r=1,k=0;
 auto f=[&](const Big&v)->Big {++out.modular_steps;return mod_pos(v*v+c,N);};
 while(g==1 && out.modular_steps<budget){
  x=y;
  for(uint64_t i=0;i<r&&out.modular_steps<budget;++i)y=f(y);
  k=0;
  while(k<r&&g==1&&out.modular_steps<budget){
   ys=y;
   uint64_t block=std::min(uint64_t(64),r-k);
   for(uint64_t j=0;j<block&&out.modular_steps<budget;++j){y=f(y);q=mod_pos(q*abs(x-y),N);}
   g=gcd_big(q,N);++out.gcd_calls;
   k+=block;
  }
  if(r>budget/2)r=budget;else r*=2;
 }
 if(g==N){
  do {ys=f(ys);g=gcd_big(abs(x-ys),N);++out.gcd_calls;}
  while(g==1 && out.modular_steps<budget);
 }
 if(proper(g,N)){out.factor=g;out.which="brent_rho";}
 return out;
}

extern "C" void dx_smooth(mpz_srcptr N,mpz_ptr factor,unsigned B,int which){
 Big n;mpz_set(n.get_mpz_t(),N);auto r=smooth_orders(n,B,which==0,which==1);mpz_set(factor,r.factor.get_mpz_t());
}
extern "C" void dx_rho(mpz_srcptr N,mpz_ptr factor,uint64_t budget){
 Big n;mpz_set(n.get_mpz_t(),N);auto r=brent_rho(n,budget);mpz_set(factor,r.factor.get_mpz_t());
}

#ifndef DYADIC_LIB
int main(int argc,char **argv){
 if(argc<2){std::cerr<<"usage: ./smooth_order N [B] [mode:hybrid|pminus|lucas|rho] [rho_budget]\n";return 2;}
 try{
  Big N(argv[1]);if(N<3 || N%2==0){std::cerr<<"N must be odd >=3\n";return 2;}
  unsigned B=argc>2 ? static_cast<unsigned>(std::stoul(argv[2])):512;
  std::string mode=argc>3 ? argv[3]:"hybrid";
  uint64_t budget=argc>4 ? static_cast<uint64_t>(std::stoull(argv[4])):1000000;
  auto start=std::chrono::steady_clock::now();
  SmoothResult out;
  if(mode=="rho")out=brent_rho(N,budget);
  else if(mode=="pminus")out=smooth_orders(N,B,true,false);
  else if(mode=="lucas")out=smooth_orders(N,B,false,true);
  else if(mode=="hybrid")out=smooth_orders(N,B,true,true);
  else {std::cerr<<"unknown mode\n";return 2;}
  auto us=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()-start).count();
  if(out.factor!=0 && (N%out.factor!=0 || !proper(out.factor,N)))return 3;
  std::cout<<"{\"N\":\""<<N<<"\",\"factor\":\""<<out.factor<<"\",\"method\":\""<<out.which<<"\",\"us\":"<<us<<",\"mod_steps\":"<<out.modular_steps<<",\"gcd_calls\":"<<out.gcd_calls<<"}\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}

#endif // DYADIC_LIB
