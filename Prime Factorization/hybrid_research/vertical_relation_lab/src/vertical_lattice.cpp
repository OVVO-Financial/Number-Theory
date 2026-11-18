#include <gmpxx.h>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using Big=mpz_class;
static Big isqrt(const Big& x){Big r; mpz_sqrt(r.get_mpz_t(),x.get_mpz_t());return r;}
static Big ceil_sqrt(const Big& x){Big r=isqrt(x);if(r*r<x)++r;return r;}
static Big floor_sum(Big n,Big m,Big a,Big b){
 if(m<=0||n<0)throw std::runtime_error("bad floor_sum");
 Big aa,ab,bb,rb;
 mpz_fdiv_qr(aa.get_mpz_t(),ab.get_mpz_t(),a.get_mpz_t(),m.get_mpz_t());
 mpz_fdiv_qr(bb.get_mpz_t(),rb.get_mpz_t(),b.get_mpz_t(),m.get_mpz_t());
 Big ans=aa*n*(n-1)/2+bb*n;a=ab;b=rb;
 while(true){
  if(a>=m){Big d=a/m;ans+=n*(n-1)/2*d;a%=m;}
  if(b>=m){Big d=b/m;ans+=n*d;b%=m;}
  Big y=a*n+b;if(y<m)break;
  n=y/m;b=y%m;std::swap(m,a);
 }
 return ans;
}
using Pair=std::pair<unsigned,unsigned>;
static std::set<Pair> residue_pairs(const Big& N,unsigned m){
 if(m==0||m>3000)throw std::runtime_error("m must be 1..3000");
 if(mpz_gcd_ui(nullptr,N.get_mpz_t(),2ul*m)!=1)throw std::runtime_error("gcd(N,2m)>1; precheck small factors");
 std::set<Pair> s;Big mod=2ul*m;
 for(unsigned p=1;p<2*m;p+=2){
  if(std::gcd(p,2u*m)!=1u)continue;
  Big inv,pp=p;
  if(mpz_invert(inv.get_mpz_t(),pp.get_mpz_t(),mod.get_mpz_t())==0)continue;
  Big q=(N*inv)%mod;
  Big y=(q*q-pp*pp)/2, h=(q*q+pp*pp)/2;
  y%=m;h%=m;if(y<0)y+=m;if(h<0)h+=m;
  s.insert({y.get_ui(),h.get_ui()});
 }
 return s;
}
static void sqrt_upper(const Big& x, int steps, Big& u, Big& v){
 u=ceil_sqrt(x);v=1;
 for(int t=0;t<steps;++t){
  Big unew=u*u+x*v*v;
  Big vnew=2*u*v;
  u=unew;v=vnew;
  // Division by gcd does not change the rational enclosure.
  Big g;mpz_gcd(g.get_mpz_t(),u.get_mpz_t(),v.get_mpz_t());
  if(g>1){u/=g;v/=g;}
 }
}
static Big residue_count(const Big& N,unsigned m,const Big& A,const Big& B,
                         const std::set<Pair>& pairs,int newton_steps=0){
 if(A<0||B<A||newton_steps<0||newton_steps>3)throw std::runtime_error("invalid interval or Newton count");
 Big N2=N*N, u0,v0,u1,v1;
 sqrt_upper(N2+A*A,newton_steps,u0,v0);
 sqrt_upper(N2+B*B,newton_steps,u1,v1);
 Big count=0,span=B-A;
 for(auto [ry,rh]:pairs){
  Big y0=A+((ry+ m -mpz_fdiv_ui(A.get_mpz_t(),m))%m);
  if(y0>B)continue;
  Big n=(B-y0)/m+1;
  if(span==0){
   // L=(N^2+A*Y)*v0/u0, U=u0/v0: exactly rigorous for any Newton depth.
   Big low_n=(N2+A*y0)*v0-Big(rh)*u0;
   Big up_n=u0-Big(rh)*v0;
   Big lo,lden=Big(m)*u0;
   mpz_cdiv_q(lo.get_mpz_t(),low_n.get_mpz_t(),lden.get_mpz_t());
   Big hi;
   Big den=m*v0;
   mpz_fdiv_q(hi.get_mpz_t(),up_n.get_mpz_t(),den.get_mpz_t());
   if(hi>=lo)count+=hi-lo+1;
   continue;
  }
  // Upper chord between exact Newton rational upper endpoints.
  // Lower left tangent: (N^2+A*Y)/sqrt(N^2+A^2) >= (N^2+A*Y)*v0/u0.
  const Big chord_d=span*v0*v1;
  const Big d_up=Big(m)*chord_d;
  const Big a_up=Big(m)*(u1*v0-u0*v1);
  const Big b_up=u0*v1*(B-y0)+u1*v0*(y0-A)-Big(rh)*chord_d;
  const Big d_low=Big(m)*u0;
  const Big a_low=-Big(m)*A*v0;
  const Big b_low=Big(rh)*u0-(N2+A*y0)*v0;
  count+=n+floor_sum(n,d_up,a_up,b_up)+floor_sum(n,d_low,a_low,b_low);
 }
 return count;
}
static bool factor_at_y(const Big& N,const Big& y,Big& p,Big& q){
 Big h2=N*N+y*y;
 if(!mpz_perfect_square_p(h2.get_mpz_t()))return false;
 Big h=isqrt(h2),p2=h-y,q2=h+y;
 if(p2<=0||!mpz_perfect_square_p(p2.get_mpz_t())||!mpz_perfect_square_p(q2.get_mpz_t()))return false;
 p=isqrt(p2);q=isqrt(q2);return p>1&&q>1&&p*q==N;
}
static Big visited=0,pruned=0, leaves=0, exhausted=0;
static bool recurse(const Big& N,unsigned m,const std::set<Pair>& pairs,const Big&A,const Big&B,
                    const Big& budget,Big& p,Big&q,int newton_steps){
 if(visited>=budget){exhausted=1;return false;}
 ++visited;
 if(residue_count(N,m,A,B,pairs,newton_steps)==0){++pruned;return false;}
 if(B-A < Big(m)*2){
  ++leaves;
  for(auto [ry,rh]:pairs){
   Big y=A+((ry+m-mpz_fdiv_ui(A.get_mpz_t(),m))%m);
   for(;y<=B;y+=m){
    Big h2=N*N+y*y;
    if(!mpz_perfect_square_p(h2.get_mpz_t()))continue;
    Big h=isqrt(h2);
    if(h%m!=rh)continue;
    if(factor_at_y(N,y,p,q))return true;
   }
  }
  return false;
 }
 Big mid=(A+B)/2;
 if(recurse(N,m,pairs,A,mid,budget,p,q,newton_steps))return true;
 if(exhausted!=0)return false;
 return recurse(N,m,pairs,mid+1,B,budget,p,q,newton_steps);
}
int main(int argc,char**argv){
 try{
  if(argc<3)throw std::runtime_error("usage: vertical_lattice <sieve|count|search|verify> N [arguments]");
  std::string action=argv[1];Big N(argv[2]);
  if(N<=0)throw std::runtime_error("N>0 required");
  if(action=="verify"){
   if(argc!=5)throw std::runtime_error("verify N p q");
   Big p(argv[3]),q(argv[4]);Big a=(p+q)/2,b=(q-p)/2;
   Big Y=2*a*b,H=a*a+b*b,E=a*a-b*b-N;
   std::cout<<"{\"valid\":"<<(p>1&&q>1&&p*q==N&&((p+q)%2==0)?"true":"false")<<",\"Y\":\""<<Y<<"\",\"H\":\""<<H<<"\",\"E\":\""<<E<<"\",\"squared_residual\":\""<<H*H-N*N-Y*Y<<"\"}\n";
   return 0;
  }
  // Factor 5 (and other wheel primes) before constructing a semiprime wheel.
  // No search for semiprimes ending in 5 is needed: N%5==0 gives a factor.
  unsigned m=argc>3?std::stoul(argv[3]):420;
  if(m==0||m>3000)throw std::runtime_error("m must be 1..3000");
  Big g; mpz_gcd_ui(g.get_mpz_t(),N.get_mpz_t(),2ul*m);
  if(N>5 && N%5==0)g=5; // trivial last-digit-5 semiprime exclusion
  if(g==N){
   for(unsigned long d=2;d<=2ul*m;++d)
    if(Big(d)<N && N%d==0){g=d;break;}
  }
  if(g>1){
   if(g<N){
    std::cout<<"{\"status\":\"precheck_factor\",\"factor\":\""<<g
             <<"\",\"cofactor\":\""<<N/g<<"\",\"wheel_not_applicable\":true}\n";
   }else{
    std::cout<<"{\"status\":\"wheel_not_applicable\",\"gcd_N_2m\":\""<<g<<"\"}\n";
   }
   return 0;
  }
  auto pairs=residue_pairs(N,m);
  if(action=="sieve"){
   std::set<unsigned> ys;for(auto [y,h]:pairs)ys.insert(y);
   std::cout<<"{\"m\":"<<m<<",\"pair_count\":"<<pairs.size()<<",\"y_residue_count\":"<<ys.size()<<",\"pairs\":[";
   bool first=true;for(auto [y,h]:pairs){if(!first)std::cout<<",";first=false;std::cout<<"["<<y<<","<<h<<"]";}
   std::cout<<"]}\n";return 0;
  }
  if(action=="count"){
   if(argc!=6&&argc!=7)throw std::runtime_error("count N m Ylo Yhi [newton_steps]");
   Big lo(argv[4]),hi(argv[5]);int steps=argc==7?std::stoi(argv[6]):0;
   std::cout<<"{\"lo\":\""<<lo<<"\",\"hi\":\""<<hi<<"\",\"envelope_lattice_count\":\""<<residue_count(N,m,lo,hi,pairs,steps)<<"\",\"newton_steps\":"<<steps<<",\"m\":"<<m<<"}\n";return 0;
  }
  if(action=="search"){
   if(argc<5||argc>8)throw std::runtime_error("search N m maxnodes [Yhi [Ylo [newton_steps]]]");
   Big budget(argv[4]), lo=argc>6?Big(argv[6]):Big(0), hi=argc>5?Big(argv[5]):3*N/4;
   int steps=argc>7?std::stoi(argv[7]):0;
   Big p=0,q=0;bool found=recurse(N,m,pairs,lo,hi,budget,p,q,steps);
   std::cout<<"{\"status\":\""<<(found?"factor_found":exhausted!=0?"budget_exhausted":"region_exhausted")<<"\",\"p\":\""<<p<<"\",\"q\":\""<<q<<"\",\"regions\":\""<<visited<<"\",\"empty\":\""<<pruned<<"\",\"leaves\":\""<<leaves<<"\"}\n";return 0;
  }
  throw std::runtime_error("unknown action");
 }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}
}
