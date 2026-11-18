// Two independent GCD channels on the same near-strip CTM states:
// E=a^2-b^2-N and H=a^2+b^2. One batch can test product(E*H) mod N.
// This is a constant-factor additional target family, not a new factoring exponent.
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
using Big=mpz_class;
static Big gcd_big(const Big&a,const Big&b){Big r;mpz_gcd(r.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return r;}
static bool proper(const Big&d,const Big&n){return d>1&&d<n;}
struct Item{Big E,H,p,q;uint64_t j;};
int main(int argc,char **argv){
 if(argc<2){std::cerr<<"usage: dual_residual N [states] [batch] [E|H|both]\n";return 2;}
 try{
  Big n(argv[1]);unsigned long width=argc>2?std::stoul(argv[2]):12000;
  unsigned long batch=argc>3?std::stoul(argv[3]):256;
  std::string mode=argc>4?argv[4]:"both";
  if(n<3||n%2==0||batch==0||!(mode=="E"||mode=="H"||mode=="both"))return 2;
  Big a;mpz_sqrt(a.get_mpz_t(),n.get_mpz_t());if(a*a<n)++a;
  Big product=1,found=0;
  std::vector<Item> items;items.reserve(batch);
  uint64_t gcd_calls=0,states=0,hit_j=0;std::string hit_source="unresolved";
  auto begin=std::chrono::steady_clock::now();
  auto flush=[&](){
   if(items.empty())return false;
   Big d=gcd_big(product,n);++gcd_calls;
   if(d==1){items.clear();product=1;return false;}
   // Localize: even when product gcd=N, either channel may carry a proper factor.
   for(const auto&it:items){
    if(it.E==0 && proper(it.p,n)&&n%it.p==0){found=it.p;hit_source="exact_strip_hit";hit_j=it.j;break;}
    if(mode!="H"){
     Big g=gcd_big(it.E,n);++gcd_calls;
     if(proper(g,n)){found=g;hit_source="E";hit_j=it.j;break;}
    }
    if(mode!="E"){
     Big g=gcd_big(it.H,n);++gcd_calls;
     if(proper(g,n)){found=g;hit_source="H";hit_j=it.j;break;}
    }
   }
   items.clear();product=1;
   return proper(found,n);
  };
  for(unsigned long j=0;j<width;++j,++a){
   Big delta=a*a-n;
   Big b;mpz_sqrt(b.get_mpz_t(),delta.get_mpz_t());
   Big E=delta-b*b,H=a*a+b*b;
   items.push_back({E,H,a-b,a+b,j});++states;
   Big v=(mode=="E")?E:((mode=="H")?H:E*H);
   v%=n;if(v<0)v+=n;
   product*=v;product%=n;
   if(items.size()>=batch && flush())break;
  }
  if(!proper(found,n))flush();
  auto us=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()-begin).count();
  if(proper(found,n)&&n%found!=0)return 3;
  std::cout<<"{\"N\":\""<<n<<"\",\"factor\":\""<<found<<"\",\"source\":\""<<hit_source<<"\",\"hit_j\":"<<hit_j<<",\"states\":"<<states<<",\"gcd_calls\":"<<gcd_calls<<",\"us\":"<<us<<"}\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}
