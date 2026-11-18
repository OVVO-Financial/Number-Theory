// Newton on S^2=1 over Z/NZ is conjugate to repeated squaring:
// t=(S-1)/(S+1), t'=t^2. For S=2ab/(a*a+b*b) we have
// t= -(a-b)^2/(a+b)^2 (mod N), provided denominators are units.
// This is a conditional 2-power-order Pollard p-1 search, NOT modular contraction.
#include <gmpxx.h>
#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
using Big=mpz_class;
static Big pos(Big x,const Big&N){x%=N;if(x<0)x+=N;return x;}
static Big gcd_big(const Big&a,const Big&b){Big r;mpz_gcd(r.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return r;}
static bool proper(const Big&d,const Big&N){return d>1&&d<N;}
struct NewtonResult{Big factor=0;unsigned k=0,a=0,b=0;uint64_t gcds=0;std::string method="unresolved";};
static NewtonResult solve(const Big&N,unsigned kmax,unsigned max_seeds){
 NewtonResult out;
 for(unsigned p: {2U,3U,5U,7U,11U})if(N%p==0){out.factor=p;out.method="small_factor";return out;}
 const std::vector<std::pair<unsigned,unsigned>> seeds={{2,1},{3,1},{3,2},{4,1},{4,3},{5,1},{5,2},{5,3},{5,4},{6,1},{7,2},{7,3}};
 for(size_t j=0;j<std::min<size_t>(seeds.size(),max_seeds);++j){
  unsigned a=seeds[j].first,b=seeds[j].second;out.a=a;out.b=b;
  Big am=a-b,ap=a+b;
  Big d=gcd_big(ap,N);++out.gcds;
  if(proper(d,N)){out.factor=d;out.method="linear_gcd";return out;}
  if(d!=1)continue;
  Big inv;mpz_invert(inv.get_mpz_t(),ap.get_mpz_t(),N.get_mpz_t());
  Big t=pos(-am*am*inv*inv,N);
  for(unsigned k=0;k<=kmax;++k){
   if(k>0)t=pos(t*t,N);
   ++out.gcds;Big g=gcd_big(t-1,N);
   if(proper(g,N)){out.factor=g;out.k=k;out.method="newton_power2_minus";return out;}
   ++out.gcds;g=gcd_big(t+1,N);
   if(proper(g,N)){out.factor=g;out.k=k;out.method="newton_power2_plus";return out;}
  }
 }
 return out;
}

extern "C" void dx_newton(mpz_srcptr N,mpz_ptr factor,unsigned kmax,unsigned seeds){
 Big n;mpz_set(n.get_mpz_t(),N);auto r=solve(n,kmax,seeds);mpz_set(factor,r.factor.get_mpz_t());
}

#ifndef DYADIC_LIB
int main(int argc,char**argv){
 if(argc<2){std::cerr<<"usage: newton2_primary N [kmax] [seeds]\n";return 2;}
 try{
  Big N(argv[1]);unsigned k=argc>2?static_cast<unsigned>(std::stoul(argv[2])):48;
  unsigned seeds=argc>3?static_cast<unsigned>(std::stoul(argv[3])):4;
  auto start=std::chrono::steady_clock::now();auto r=solve(N,k,seeds);
  auto us=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()-start).count();
  if(r.factor!=0 && (!proper(r.factor,N)||N%r.factor!=0))return 3;
  std::cout<<"{\"N\":\""<<N<<"\",\"factor\":\""<<r.factor<<"\",\"method\":\""<<r.method<<"\",\"k\":"<<r.k<<",\"a\":"<<r.a<<",\"b\":"<<r.b<<",\"gcd_calls\":"<<r.gcds<<",\"us\":"<<us<<"}\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}

#endif // DYADIC_LIB
