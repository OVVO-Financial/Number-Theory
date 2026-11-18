// In-process exact GMP hybrid. No subprocess overhead; optional bounded rho fallback.
#include <gmpxx.h>
#include <iostream>
#include <chrono>
#include <string>
#include <cstdint>
#include <stdexcept>
using Big=mpz_class;
extern "C" void dx_newton(mpz_srcptr N,mpz_ptr factor,unsigned kmax,unsigned seeds);
extern "C" void dx_gaussian(mpz_srcptr N,mpz_ptr factor,unsigned B,unsigned seeds,bool near_line);
extern "C" void dx_smooth(mpz_srcptr N,mpz_ptr factor,unsigned B,int which);
extern "C" void dx_rho(mpz_srcptr N,mpz_ptr factor,uint64_t budget);
int main(int argc,char**argv){
 if(argc<2){std::cerr<<"usage: hybrid_native N [mode:hybrid|rho|gauss|newton|pminus|lucas] [B] [rho_budget]\n";return 2;}
 try{
  Big N(argv[1]);if(N<3){std::cerr<<"N must be >=3\n";return 2;}
  std::string mode=argc>2?argv[2]:"hybrid";
  unsigned B=argc>3?static_cast<unsigned>(std::stoul(argv[3])):256;
  uint64_t budget=argc>4?static_cast<uint64_t>(std::stoull(argv[4])):400000;
  Big found=0;
  auto now=std::chrono::steady_clock::now;
  auto start=now();
  std::string method="unresolved";
  std::string traces;
  auto stage=[&](const std::string&name,auto fn){
   if(found>1 && found<N)return;
   auto st=now();fn();
   auto dt=std::chrono::duration_cast<std::chrono::microseconds>(now()-st).count();
   if(!traces.empty())traces+=",";
   traces+="{\"stage\":\""+name+"\",\"us\":"+std::to_string(dt)+"}";
   if(found>1 && found<N){if(N%found!=0)throw std::runtime_error("invalid factor");method=name;}
  };
  if(mode=="hybrid"||mode=="newton")stage("newton2",[&](){dx_newton(N.get_mpz_t(),found.get_mpz_t(),48,4);});
  if(mode=="hybrid"||mode=="gauss")stage("gaussian",[&](){dx_gaussian(N.get_mpz_t(),found.get_mpz_t(),B*2,4,false);});
  if(mode=="hybrid"||mode=="pminus")stage("pminus",[&](){dx_smooth(N.get_mpz_t(),found.get_mpz_t(),B,0);});
  if(mode=="hybrid"||mode=="lucas")stage("lucas",[&](){dx_smooth(N.get_mpz_t(),found.get_mpz_t(),B,1);});
  if(mode=="hybrid"||mode=="rho")stage("rho",[&](){dx_rho(N.get_mpz_t(),found.get_mpz_t(),budget);});
  auto dt=std::chrono::duration_cast<std::chrono::microseconds>(now()-start).count();
  std::cout<<"{\"N\":\""<<N<<"\",\"factor\":\""<<found<<"\",\"method\":\""<<method<<"\",\"us\":"<<dt<<",\"stages\":["<<traces<<"]}"<<std::endl;
  return 0;
 }catch(const std::exception&e){std::cerr<<"Error: "<<e.what()<<"\n";return 2;}
}
