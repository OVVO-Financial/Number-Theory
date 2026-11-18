// Compare mathematically valid square-congruence relation *candidates*, not just small residuals.
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using Big=mpz_class;
struct Fact {bool smooth=false,negative=false;std::vector<unsigned short> exp;};
struct Relation {std::vector<unsigned short> left,right;bool neg=false;Big x;};
static Big root(const Big& n){Big x;mpz_sqrt(x.get_mpz_t(),n.get_mpz_t());return x;}
static Big ceilroot(const Big& n){Big x=root(n);if(x*x<n)++x;return x;}
static std::vector<unsigned> primes_up_to(unsigned B){
 std::vector<bool> sieve(B+1,true);sieve[0]=false;if(B>=1)sieve[1]=false;
 for(unsigned p=2;p*p<=B;++p)if(sieve[p])for(unsigned q=p*p;q<=B;q+=p)sieve[q]=false;
 std::vector<unsigned> v;for(unsigned p=2;p<=B;++p)if(sieve[p])v.push_back(p);return v;
}
static Fact smooth(const Big& value,const std::vector<unsigned>& primes){
 Fact r;r.negative=value<0;r.exp.assign(primes.size(),0);
 if(value==0)return r;
 Big v=abs(value);
 for(size_t j=0;j<primes.size();++j){
  unsigned p=primes[j];
  while(mpz_divisible_ui_p(v.get_mpz_t(),p)!=0){
   mpz_divexact_ui(v.get_mpz_t(),v.get_mpz_t(),p);
   if(r.exp[j]==std::numeric_limits<unsigned short>::max())throw std::runtime_error("exponent overflow");
   ++r.exp[j];
  }
  if(v==1)break;
 }
 r.smooth=(v==1);return r;
}
static void xor_words(std::vector<uint64_t>&a,const std::vector<uint64_t>&b){for(size_t i=0;i<a.size();++i)a[i]^=b[i];}
static bool bit(const std::vector<uint64_t>&a,size_t i){return (a[i/64]>>(i%64))&1;}
struct Engine{
 bool one_sided;
 Big N;std::vector<unsigned> primes;std::string mode;
 unsigned long long valid=0, raw_smooth=0,direct=0, rank=0, dep=0, attempts=0, factor_deps=0, duplicate_residuals=0;
 std::set<Big> residuals;
 Big extracted=0;
 size_t cols=0,words=0,witness_words=0;
 std::vector<std::vector<uint64_t>> basis,basis_witness;
 std::vector<Relation> rel;
 Engine(Big n,std::vector<unsigned> f,std::string mode_,size_t max_candidates):one_sided(mode_=="QS"||mode_=="CTMWalkOneSide"||mode_=="GroupingBandOneSide"||mode_=="VerticalOneSide"),N(n),primes(std::move(f)),mode(std::move(mode_)){
  cols=one_sided?primes.size()+1:2*primes.size()+1;
  words=(cols+63)/64;witness_words=(max_candidates+63)/64;
  basis.resize(cols);basis_witness.resize(cols);
 }
 std::vector<uint64_t> parity(const Relation&r){
  std::vector<uint64_t> bits(words,0);
  size_t off=one_sided?0:primes.size();
  if(!one_sided)for(size_t i=0;i<primes.size();++i)if(r.left[i]&1)bits[i/64]^=uint64_t(1)<<(i%64);
  for(size_t i=0;i<primes.size();++i)if(r.right[i]&1){size_t col=off+i;bits[col/64]^=uint64_t(1)<<(col%64);}
  if(r.neg){size_t col=cols-1;bits[col/64]^=uint64_t(1)<<(col%64);}
  return bits;
 }
 void dependency(const std::vector<uint64_t>&w){
  if(attempts>=64)return;
  ++attempts;
  std::vector<uint64_t> lexp(primes.size(),0),rexp(primes.size(),0);
  Big X=1,Y=1;unsigned negatives=0;
  for(size_t j=0;j<rel.size();++j){
   if(!bit(w,j))continue;
   const auto& r=rel[j];negatives+=r.neg;
   if(one_sided)X=(X*r.x)%N;
   else for(size_t i=0;i<primes.size();++i)lexp[i]+=r.left[i];
   for(size_t i=0;i<primes.size();++i)rexp[i]+=r.right[i];
  }
  if(negatives%2)throw std::runtime_error("invalid negative sign dependency");
  for(size_t i=0;i<primes.size();++i){
   if(rexp[i]%2||lexp[i]%2)throw std::runtime_error("invalid parity dependency");
   Big base=primes[i],pp;
   if(!one_sided&&lexp[i]){mpz_powm_ui(pp.get_mpz_t(),base.get_mpz_t(),lexp[i]/2,N.get_mpz_t());X=(X*pp)%N;}
   if(rexp[i]){mpz_powm_ui(pp.get_mpz_t(),base.get_mpz_t(),rexp[i]/2,N.get_mpz_t());Y=(Y*pp)%N;}
  }
  if((X*X-Y*Y)%N!=0)throw std::runtime_error("relation dependency invalid modulo N");
  Big d;
  Big diff=X-Y; mpz_gcd(d.get_mpz_t(),diff.get_mpz_t(),N.get_mpz_t());
  if(d>1&&d<N){++factor_deps;extracted=d;return;}
  Big summ=X+Y; mpz_gcd(d.get_mpz_t(),summ.get_mpz_t(),N.get_mpz_t());
  if(d>1&&d<N){++factor_deps;extracted=d;}
 }
 void offer(const Big&left,const Big&right,const Big& x=0){
  if(right==0){
   ++direct;
   if(x>1 && x<N && N%x==0)extracted=x;
   return;
  }
  if(!one_sided && !residuals.insert(right).second){++duplicate_residuals;return;}
  Fact rf=smooth(right,primes);
  if(!rf.smooth)return;
  ++raw_smooth;
  Fact lf;
  if(one_sided){
   lf.smooth=true;lf.exp.assign(primes.size(),0);
  }else{
   if(left<=0)return;
   lf=smooth(left,primes);if(!lf.smooth)return;
  }
  if((left-right)%N!=0)throw std::runtime_error("invalid modular relation");
  Relation r{std::move(lf.exp),std::move(rf.exp),rf.negative,x};
  ++valid;rel.push_back(std::move(r));
  auto bits=parity(rel.back());
  std::vector<uint64_t>w(witness_words,0);
  size_t id=rel.size()-1;w[id/64]|=uint64_t(1)<<(id%64);
  while(true){
   size_t col=0;for(;col<cols;++col)if(bit(bits,col))break;
   if(col==cols){++dep;dependency(w);break;}
   if(basis[col].empty()){
    basis[col]=std::move(bits);basis_witness[col]=std::move(w);++rank;break;
   }
   xor_words(bits,basis[col]);xor_words(w,basis_witness[col]);
  }
 }
};
int main(int argc,char**argv){
 try{
  if(argc!=5 && argc!=6)throw std::runtime_error("usage: relation_bench N candidate_count smooth_bound QS|CTM|CTMWalk|Grouping|GroupingBand|SquareY|NearPlus|NearMinus|CTMWalkOneSide|GroupingBandOneSide|VerticalOneSide [--stop-on-factor]");
  bool stop=(argc==6&&std::string(argv[5])=="--stop-on-factor");
  Big N(argv[1]);
  if(N<=5||N%2==0)throw std::runtime_error("N must be odd >5");
  unsigned long K=std::stoul(argv[2]),B=std::stoul(argv[3]);std::string mode=argv[4];
  if(!K||K>250000||B<2||B>10000)throw std::runtime_error("K or B outside bounds");
  if(mode!="QS"&&mode!="CTM"&&mode!="CTMWalk"&&mode!="Grouping"&&mode!="GroupingBand"&&mode!="SquareY"&&mode!="NearPlus"&&mode!="NearMinus"&&mode!="CTMWalkOneSide"&&mode!="GroupingBandOneSide"&&mode!="VerticalOneSide")throw std::runtime_error("unknown mode");
  auto ps=primes_up_to(B);Engine eng(N,ps,mode,K);
  Big a0=ceilroot(N), low=a0-1, hi=ceilroot(2*N)+1;
  struct Walk {Big a,b,e;};
  std::vector<Walk> walks;
  if(mode=="CTMWalk"||mode=="CTMWalkOneSide"){
   unsigned nm=mpz_fdiv_ui(N.get_mpz_t(),20);
   unsigned am=mpz_fdiv_ui(a0.get_mpz_t(),10);
   for(int ar=0;ar<10;++ar)for(int br=0;br<10;++br){
    int dif=((ar*ar-br*br-(int)nm)%20+20)%20;
    if(dif!=0)continue;
    Big a=a0+((ar+10-(int)am)%10),b=br;
    if(a<=b+1)continue;
    Big e=a*a-b*b-N;
    walks.push_back({a,b,e});
   }
   if(walks.empty())throw std::runtime_error("No CTM lanes found");
  }
  auto t0=std::chrono::steady_clock::now();
  unsigned long evaluated=0;std::set<Big> seen_x;
  for(unsigned long j=0;j<K;++j){
   if(mode=="QS"||mode=="CTM"){
    ++evaluated;
    Big a=a0+j, f=a*a-N;
    if(mode=="QS"){
     if(mpz_perfect_square_p(f.get_mpz_t())){
      Big bb=root(f),d=a-bb;
      if(d>1&&d<N&&N%d==0){++eng.direct;eng.extracted=d;}
     }else eng.offer(a*a,f,a);
    }
    else{
     Big b=root(f), E=f-b*b;
     // In CTM coordinates L=a-b,U=a+b and L*U=N+E.
     eng.offer((a-b)*(a+b),E,a-b);
    }
   }else if(mode=="CTMWalk"||mode=="CTMWalkOneSide"){
    auto &w=walks[j%walks.size()];
    if(w.a-w.b>1){
     ++evaluated;
     Big left=(w.a-w.b)*(w.a+w.b);
     if(left-w.e!=N)throw std::runtime_error("CTM staircase residual drift");
     if(mode=="CTMWalkOneSide"){
      Big x=w.a;
      // Distinct x values, not repeated staircase states, define independent QS trials.
      if(seen_x.insert(x).second){
       Big f=x*x-N;
       if(mpz_perfect_square_p(f.get_mpz_t())){Big b=root(f),d=x-b;
        if(d>1&&d<N&&N%d==0){++eng.direct;eng.extracted=d;}
       }else eng.offer(x*x,f,x);
      }
     }else eng.offer(left,w.e,w.a-w.b);
     if(w.e<0){Big change=20*w.a+100;w.a+=10;w.e+=change;}
     else {Big change=20*w.b+100;w.b+=10;w.e-=change;}
    }
   }else if(mode=="Grouping"){
    Big g=2*j+3; ++evaluated; Big quo=N/g,r=N%g;
    eng.offer(g*quo,-r,g);
   }else if(mode=="GroupingBand"||mode=="GroupingBandOneSide"){
    Big g=low+((4-mpz_fdiv_ui(low.get_mpz_t(),4))%4)+4*j;
    if(g>hi)break;
    ++evaluated;
    if(mode=="GroupingBandOneSide"){
     for(int direction : {-1,1}){
      Big q=g+direction;if(q<=1)continue;
      Big x=(q+N/q+1)/2;if(x<a0)x=a0;
      if(!seen_x.insert(x).second)continue;
      Big f=x*x-N;
      if(mpz_perfect_square_p(f.get_mpz_t())){Big b=root(f),d=x-b;
       if(d>1&&d<N&&N%d==0){++eng.direct;eng.extracted=d;}
      }else eng.offer(x*x,f,x);
     }
    }else{Big quo=N/g,r=N%g;eng.offer(g*quo,-r,g);}
   }else if(mode=="VerticalOneSide"){
    ++evaluated;
    Big Y=Big(3)*N*(j+1)/(Big(4)*(K+1));
    // N-only inverse of the monotone vertical-line map to the real Fermat axis.
    Big h=ceilroot(N*N+Y*Y);
    Big x=ceilroot((N+h+1)/2);
    if(x<a0)x=a0;
    if(seen_x.insert(x).second){
     Big f=x*x-N;
     if(mpz_perfect_square_p(f.get_mpz_t())){Big b=root(f),d=x-b;
      if(d>1&&d<N&&N%d==0){++eng.direct;eng.extracted=d;}
     }else eng.offer(x*x,f,x);
    }
   }else if(mode=="NearPlus" || mode=="NearMinus"){
    ++evaluated;
    Big E=j+1;if(mode=="NearMinus")E=-E;
    eng.offer(N+E,E);
   }else if(mode=="SquareY"){
    Big Y=Big(3)*N*(j+1)/(Big(4)*(K+1));Y-=Y%4;
    if(Y==0)continue;
    ++evaluated;
    Big h=root(N*N+Y*Y),l=(h-Y)*(h+Y),E=l-N*N;
    eng.offer(l,E);
   }
   if(stop&&eng.extracted>1)break;
  }
  auto micros=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now()-t0).count();
  std::cout<<"{\"method\":\""<<mode<<"\",\"N\":\""<<N<<"\",\"smooth_bound\":"<<B<<",\"candidates_requested\":"<<K<<",\"candidates_evaluated\":"<<evaluated<<",\"factor_base_primes\":"<<ps.size()<<",\"unique_residuals\":"<<eng.residuals.size()<<",\"unique_one_sided_x\":"<<seen_x.size()<<",\"duplicate_residuals\":"<<eng.duplicate_residuals<<",\"residual_smooth\":"<<eng.raw_smooth<<",\"valid_relations\":"<<eng.valid<<",\"rank\":"<<eng.rank<<",\"dependencies\":"<<eng.dep<<",\"direct_exact_hits\":"<<eng.direct<<",\"factor_from_dependency\":"<<(eng.factor_deps>0?"true":"false")<<",\"dependency_attempts\":"<<eng.attempts<<",\"factor_recovered\":\""<<eng.extracted<<"\",\"time_us\":"<<micros<<"}\n";
 }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}
}
