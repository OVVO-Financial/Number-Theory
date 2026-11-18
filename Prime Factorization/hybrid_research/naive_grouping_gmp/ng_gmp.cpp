// Exact GMP naive-grouping seed chains: an independently verifiable research worker.
// A failed bounded search is never called a primality proof.
#include <gmpxx.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
using Big=mpz_class;
static Big isqrt(const Big &n){Big r;mpz_sqrt(r.get_mpz_t(),n.get_mpz_t());return r;}
static Big gcd_(const Big &a,const Big &b){Big g;mpz_gcd(g.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return g;}
struct Answer {Big factor=0, seed=0, group=0; unsigned long long starts=0,states=0;};
static bool proper(const Big &n,const Big &d){return d>1&&d<n&&n%d==0;}
static bool test_group(const Big &N,const Big &g,bool gcd_mode,Big &factor){
 if(g<=1||g>=N)return false;
 Big a=N/g,b=N%g,c=g-b;
 if(b==0 && proper(N,g)){factor=g;return true;}
 if(b>1 && a%b==0 && proper(N,b)){factor=b;return true;}
 if(c>1 && (a+1)%c==0 && proper(N,c)){factor=c;return true;}
 if(gcd_mode){
  Big d=gcd_(a,b);if(proper(N,d)){factor=d;return true;}
  d=gcd_(a+1,c);if(proper(N,d)){factor=d;return true;}
 }
 return false;
}
// On doubling: floor(N/(2g))=floor(a/2); N mod (2g)=b+(a mod 2)*g.
static bool chain(const Big &N,const Big &seed,const Big &maxg,unsigned long long maxstates,
                  bool gcd_mode,Answer &ans){
 if(seed<=0||seed>=N)return false;
 Big g=seed,a=N/g,b=N%g;
 ++ans.starts;
 while(g<=maxg && ans.states<maxstates){
  ++ans.states;
  Big c=g-b,d=0;
  if(b==0&&proper(N,g))d=g;
  else if(b>1 && a%b==0 && proper(N,b))d=b;
  else if(c>1 && (a+1)%c==0 && proper(N,c))d=c;
  if(d==0 && gcd_mode){d=gcd_(a,b);if(!proper(N,d))d=gcd_(a+1,c);if(!proper(N,d))d=0;}
  if(d>1){ans.factor=d;ans.seed=seed;ans.group=g;return true;}
  if(g>maxg/2)break;
  b+=(a%2)*g;a/=2;g*=2;
 }
 return false;
}
int main(int argc,char **argv){
 try{
  if(argc<3)throw std::runtime_error("usage: ng_gmp N seed|seeds|groups [seed_or_maxstarts] [maxstates] [--gcd]");
  Big N(argv[1]);if(N<=1)throw std::runtime_error("N must be >1");
  std::string mode=argv[2]; bool gcd_mode=false;
  for(int j=3;j<argc;++j)if(std::string(argv[j])=="--gcd")gcd_mode=true;
  Big factor=0;Answer ans;
  if(N%2==0)factor=2;
  if(factor==0 && N>5 && N%5==0)factor=5;
  const Big sqrtN=isqrt(N), half=sqrtN/2, qhi=isqrt(2*N)+2;
  const unsigned long long maxstates=(argc>4&&argv[4][0]!='-')?std::stoull(argv[4]):1000000ull;
  if(factor==0 && mode=="seed"){
   if(argc<4)throw std::runtime_error("seed mode requires seed");
   Big seed(argv[3]);chain(N,seed,std::max(qhi,Big(2*sqrtN)),maxstates,gcd_mode,ans);
  }else if(factor==0 && mode=="seeds"){
   unsigned long long maxstarts=argc>3&&argv[3][0]!='-'?std::stoull(argv[3]):100000;
   for(Big s=1;s<half && ans.starts<maxstarts && ans.states<maxstates;s+=2)
    if(chain(N,s,std::max(qhi,Big(2*sqrtN)),maxstates,gcd_mode,ans))break;
  }else if(factor==0 && mode=="groups"){
   // CONDITIONAL: complete for the nontrivial regrouping success of p<q<2p.
   // Group q-1 or q+1 divisible by four is inside the range below.
   Big lo=sqrtN>2?sqrtN-2:Big(4),hi=qhi;
   Big g=lo+((4-mpz_fdiv_ui(lo.get_mpz_t(),4))%4);
   unsigned long long maxgroups=argc>3&&argv[3][0]!='-'?std::stoull(argv[3]):100000;
   while(g<=hi && ans.states<maxstates && ans.starts<maxgroups){
    ++ans.starts;++ans.states;
    if(test_group(N,g,gcd_mode,ans.factor)){ans.group=g;break;}
    g+=4;
   }
  }else if(factor==0)throw std::runtime_error("mode must be seed, seeds, or groups");
  if(factor==0)factor=ans.factor;
  if(factor!=0 && !proper(N,factor))throw std::runtime_error("nontrivial factor validation failed");
  std::cout<<"{\"status\":\""<<(factor>1?"factor_found":"budget_or_region_exhausted")
   <<"\",\"factor\":\""<<factor<<"\",\"cofactor\":\""<<(factor>1?N/factor:Big(0))
   <<"\",\"seed\":\""<<ans.seed<<"\",\"group\":\""<<ans.group
   <<"\",\"starts\":"<<ans.starts<<",\"states\":"<<ans.states<<"}\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}
