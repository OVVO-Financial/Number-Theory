// ==== fermat_dual_gpu.cu ====
// CPU-only build (works without CUDA):
//   g++ -x c++ -O3 -std=gnu++17 -fopenmp fermat_dual_gpu.cu -lgmpxx -lgmp -o fermat_dual
//
// GPU build (CUDA sieve enabled; requires nvcc):
//   nvcc -O3 -std=c++17 -DUSE_CUDA fermat_dual_gpu.cu -Xcompiler "-O3 -fopenmp" -lgmpxx -lgmp -o fermat_dual_gpu
//
// Run examples:
//   ./fermat_dual --stream=dual --vspace=1
//   ./fermat_dual --stream=a    --vspace=2
//   ./fermat_dual --stream=b
//
// --stream=dual|a|b   : choose both streams (race), or a-only, or b-only
// --vspace=1|2        : a-stream spacing: 1=uniform, 2=exponential bands

#include <gmpxx.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <cstdint>
#include <atomic>
#include <cstring>

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef USE_CUDA
  #include <cuda_runtime.h>
#endif

// ===== Helpers =====
static inline bool is_square(const mpz_class& x){ return mpz_perfect_square_p(x.get_mpz_t()) != 0; }
static inline mpz_class isqrt(const mpz_class& x){ mpz_class r; mpz_sqrt(r.get_mpz_t(), x.get_mpz_t()); return r; }
static inline long lcm_long(long a,long b){ return std::lcm(a,b); }
static inline unsigned long umod(const mpz_class& n, unsigned long m){ return mpz_tdiv_ui(n.get_mpz_t(), m); }
static std::vector<bool> quad_residues_mod(long m){ std::vector<bool>qr(m,false); for(long r=0;r<m;++r) qr[(1LL*r*r)%m]=true; return qr; }
static inline mpz_class mpz_from_u64(uint64_t x){ mpz_class z; mpz_import(z.get_mpz_t(),1,-1,sizeof(x),0,0,&x); return z; }

// SAFE 64-bit export for mpz_class (no truncation on Windows)
static inline uint64_t mpz_to_u64(const mpz_class& z){
  uint64_t v = 0; size_t count = 0;
  mpz_export(&v, &count, -1, sizeof(uint64_t), 0, 0, z.get_mpz_t());
  return v;
}

// Ceil/Floor division by small unsigned long, returning uint64_t safely
static inline uint64_t ceil_div_mpz_ui(const mpz_class& num, unsigned long den){
  mpz_class q = (num + (den - 1)) / den;   // ceil(num/den) for nonnegative num
  return mpz_to_u64(q);
}
static inline uint64_t floor_div_mpz_ui(const mpz_class& num, unsigned long den){
  mpz_class q = num / den;                 // floor(num/den) for nonnegative num
  return mpz_to_u64(q);
}

// ===== Tiny TD (wheel-30) =====
struct StageRes{ bool found=false; mpz_class p,q; size_t work=0; double secs=0.0; };
static StageRes tiny_td(const mpz_class& n, long limit=500000){
  StageRes r; auto t0=std::chrono::high_resolution_clock::now();
  if(mpz_divisible_ui_p(n.get_mpz_t(),2)){ r.found=true; r.p=2; r.q=n/2; }
  else if(mpz_divisible_ui_p(n.get_mpz_t(),3)){ r.found=true; r.p=3; r.q=n/3; }
  else if(mpz_divisible_ui_p(n.get_mpz_t(),5)){ r.found=true; r.p=5; r.q=n/5; }
  size_t tested=0;
  if(!r.found){
    static const int OFF[8]={1,-1,7,-7,11,-11,13,-13};
    for(long k=1;30L*k-1<=limit;++k){
      long base=30L*k;
      for(int i=0;i<8;++i){
        long p=base+OFF[i]; if(p>limit) break;
        ++tested;
        if(mpz_divisible_ui_p(n.get_mpz_t(), (unsigned long)p)){ r.found=true; r.p=p; r.q=n/r.p; break; }
      }
      if(r.found) break;
    }
  }
  r.work=tested; r.secs=std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t0).count();
  return r;
}

// ===== Iterated Averages (IA) =====
static StageRes ia_probe(const mpz_class& n){
  StageRes r; auto t0=std::chrono::high_resolution_clock::now();
  if(n<=3){ r.secs=0.0; return r; }
  mpz_class a0; mpz_sqrt(a0.get_mpz_t(), n.get_mpz_t()); if(a0*a0<n) a0 += 1;

  int n_last = mpz_class(n % 10).get_si();
  std::vector<long> offs = {0,-1,1};
  if(n_last==1 || n_last==9){ offs.push_back(-2); offs.push_back(2); }

  size_t checks=0;
  auto try_g = [&](const mpz_class& x)->mpz_class{
    if(x<=0) return 0; ++checks; mpz_class g; mpz_gcd(g.get_mpz_t(), n.get_mpz_t(), x.get_mpz_t());
    return (g>1 && g<n) ? g : 0;
  };

  const int MAX=32;
  for(int i=1;i<=MAX;++i){
    mpz_class num; mpz_mul_2exp(num.get_mpz_t(), n.get_mpz_t(), i); num += (a0 - n);
    mpz_class f;   mpz_tdiv_q_2exp(f.get_mpz_t(), num.get_mpz_t(), i);
    mpz_class c = num + (mpz_class(1) << i) - 1; mpz_tdiv_q_2exp(c.get_mpz_t(), c.get_mpz_t(), i);
    for(long o:offs){
      if(auto g=try_g(f+o)){ r.found=true; r.p=g; r.q=n/g; break; }
      if(r.found) break;
      if(auto g=try_g(c+o)){ r.found=true; r.p=g; r.q=n/g; break; }
    }
    if(r.found) break;
  }
  r.work=checks; r.secs=std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t0).count();
  return r;
}

// ===== On-the-fly QR sieve =====
struct OnFlyQR{
  std::vector<int> primes;
  std::vector<std::vector<char>> qr;
  void init(const std::vector<int>& ps){
    primes=ps; qr.clear(); qr.reserve(primes.size());
    for(int p:primes){ std::vector<char> v(p,0); for(int r=0;r<p;++r) v[(1LL*r*r)%p]=1; qr.emplace_back(std::move(v)); }
  }
  inline bool pass_b(const mpz_class& n, const mpz_class& b) const{
    for(size_t i=0;i<primes.size();++i){
      int p=primes[i]; unsigned long nm=umod(n,p), br=umod(b,p);
      unsigned long a2=(nm + ( (1ULL*br*br)%p )) % p;
      if(!qr[i][a2]) return false;
    }
    return true;
  }
  inline bool pass_a(const mpz_class& n, const mpz_class& a) const{
    for(size_t i=0;i<primes.size();++i){
      int p=primes[i]; unsigned long nm=umod(n,p), ar=umod(a,p);
      unsigned long d = ( (1ULL*ar*ar)%p + p - nm ) % p;
      if(!qr[i][d]) return false;
    }
    return true;
  }
};

// ===== Deterministic winner =====
struct Winner{ bool valid=false; mpz_class b,a,p,q; long r=0; int stream=0; /*0=b,1=a*/ };
static inline bool winner_less(const Winner& x, const Winner& y){
  if(!y.valid) return true;
  int cmp = mpz_cmp(x.b.get_mpz_t(), y.b.get_mpz_t());
  if(cmp<0) return true; if(cmp>0) return false; return x.r < y.r;
}

// ===== Stats =====
struct FermatStats{
  bool found=false; mpz_class a,b,p,q,k2b;
  size_t candidates_b=0, candidates_a=0;
  size_t allow_b=0, allow_a=0; long M=0;
  int stream=-1; // 0=b, 1=a
  double secs=0.0;
};

// ===== GPU sieve (optional) =====
#ifdef USE_CUDA
__constant__ uint32_t d_P[128];
__constant__ uint32_t d_P_off[128];

__global__ void sieve_b_kernel(
  uint64_t k0, uint64_t k1, uint64_t M, uint64_t bmax,
  int npr, int nallow,
  const uint32_t* __restrict__ n_mod_p,
  const uint32_t* __restrict__ M_mod_p,
  const uint32_t* __restrict__ r_mod_p,
  const uint32_t* __restrict__ p_off,
  const uint8_t*  __restrict__ QR_flat,
  const uint64_t* __restrict__ allow_b,
  uint64_t* __restrict__ out_b,
  unsigned int* __restrict__ out_count)
{
  const uint64_t i = blockIdx.y * blockDim.y + threadIdx.y; // residue index
  const uint64_t k = k0 + blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= (uint64_t)nallow || k > k1) return;

  const uint64_t r = allow_b[i];
  const uint64_t b = k * M + r;
  if (b > bmax) return;

  bool ok = true;
  #pragma unroll 4
  for (int j = 0; j < npr && ok; ++j){
    const uint32_t p   = d_P[j];
    const uint32_t km  = (uint32_t)(k % p);
    const uint32_t Mp  = M_mod_p[j];
    const uint32_t rm  = r_mod_p[p_off[j] + i];
    const uint32_t bp  = ( (uint64_t)km * Mp + rm ) % p;
    const uint32_t a2  = ( (uint64_t)bp * bp + n_mod_p[j] ) % p;
    const uint32_t off = d_P_off[j];
    if (QR_flat[off + a2] == 0) ok = false;
  }
  if (!ok) return;

  unsigned int idx = atomicAdd(out_count, 1u);
  out_b[idx] = b;
}

struct GPUBandResult { std::vector<uint64_t> b_hits; };

static GPUBandResult gpu_sieve_b_band(
  uint64_t k0, uint64_t k1, uint64_t M, uint64_t bmax,
  const std::vector<int>& primes,
  const std::vector<uint8_t>& QR_flat,
  const std::vector<uint32_t>& P_off,
  const std::vector<uint32_t>& n_mod_p,
  const std::vector<uint32_t>& M_mod_p,
  const std::vector<uint64_t>& allow_b,
  const std::vector<uint32_t>& r_mod_p_flat,
  const std::vector<uint32_t>& p_off_flat,
  int threads_x = 256, int threads_y = 4)
{
  GPUBandResult R;
  const int npr = (int)primes.size();
  const int nallow = (int)allow_b.size();
  if(nallow==0 || npr==0 || k0>k1) return R;

  // device buffers
  uint32_t *d_nmod=nullptr, *d_Mmod=nullptr, *d_rmod=nullptr, *d_poff=nullptr;
  uint64_t *d_allow=nullptr, *d_outb=nullptr;
  uint8_t  *d_QR=nullptr;
  unsigned int *d_count=nullptr;

  cudaMalloc(&d_nmod,  npr * sizeof(uint32_t));
  cudaMalloc(&d_Mmod,  npr * sizeof(uint32_t));
  cudaMalloc(&d_rmod,  r_mod_p_flat.size() * sizeof(uint32_t));
  cudaMalloc(&d_poff,  p_off_flat.size() * sizeof(uint32_t));
  cudaMalloc(&d_allow, nallow * sizeof(uint64_t));
  cudaMalloc(&d_QR,    QR_flat.size() * sizeof(uint8_t));

  const uint64_t kspan = k1 - k0 + 1;
  const size_t cap = (size_t)std::min<uint64_t>((uint64_t)nallow * kspan, 50'000'000ULL);
  cudaMalloc(&d_outb,  cap * sizeof(uint64_t));
  cudaMalloc(&d_count, sizeof(unsigned int));
  cudaMemset(d_count, 0, sizeof(unsigned int));

  cudaMemcpy(d_nmod,  n_mod_p.data(), npr*sizeof(uint32_t), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Mmod,  M_mod_p.data(), npr*sizeof(uint32_t), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rmod,  r_mod_p_flat.data(), r_mod_p_flat.size()*sizeof(uint32_t), cudaMemcpyHostToDevice);
  cudaMemcpy(d_poff,  p_off_flat.data(), p_off_flat.size()*sizeof(uint32_t), cudaMemcpyHostToDevice);
  cudaMemcpy(d_allow, allow_b.data(), nallow*sizeof(uint64_t), cudaMemcpyHostToDevice);
  cudaMemcpy(d_QR,    QR_flat.data(), QR_flat.size()*sizeof(uint8_t), cudaMemcpyHostToDevice);

  std::vector<uint32_t> P_u(primes.begin(), primes.end());
  cudaMemcpyToSymbol(d_P,     P_u.data(), npr*sizeof(uint32_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_P_off, P_off.data(), npr*sizeof(uint32_t), 0, cudaMemcpyHostToDevice);

  dim3 block(threads_x, threads_y);
  const uint64_t nx = (kspan + threads_x - 1) / threads_x;
  const uint64_t ny = (nallow + threads_y - 1) / threads_y;
  dim3 grid((unsigned)nx, (unsigned)ny);

  sieve_b_kernel<<<grid, block>>>(
    k0, k1, M, bmax, npr, nallow,
    d_nmod, d_Mmod, d_rmod, d_poff, d_QR,
    d_allow, d_outb, d_count
  );
  cudaDeviceSynchronize();

  unsigned int count=0;
  cudaMemcpy(&count, d_count, sizeof(unsigned int), cudaMemcpyDeviceToHost);
  R.b_hits.resize(count);
  if (count>0) cudaMemcpy(R.b_hits.data(), d_outb, count*sizeof(uint64_t), cudaMemcpyDeviceToHost);

  cudaFree(d_nmod); cudaFree(d_Mmod); cudaFree(d_rmod); cudaFree(d_poff);
  cudaFree(d_allow); cudaFree(d_QR); cudaFree(d_outb); cudaFree(d_count);
  return R;
}
#endif // USE_CUDA

// ===== CRT filters & QR primes =====
struct BuildCRT {
  long M=1;
  std::vector<long> mods;
  struct Filter{ long m; std::vector<char> ok_b; std::vector<char> ok_a; };
  std::vector<Filter> filt;
  std::vector<long> allow_b, allow_a;
};

static BuildCRT build_crt_filters(const mpz_class& n){
  BuildCRT B;
  const long MAX_M = 50'000'000L; // cap stride
  std::vector<long> base_mods = {10, 16, 9, 5, 7, 11, 13};
  auto try_add = [&](long m){
    long newM = lcm_long(B.M, m);
    if(newM>0 && newM<=MAX_M){ B.M=newM; B.mods.push_back(m); return true; }
    return false;
  };
  for(long m: base_mods) (void)try_add(m);

  int log2_n = mpz_sizeinbase(n.get_mpz_t(), 2);
  int max_crt_prime = std::min(47, log2_n / 2);
  const long extras[] = {17,19,23,29,31,37,41,43,47};
  for(long p: extras) if(p<=max_crt_prime) (void)try_add(p);

  // build filters
  B.filt.reserve(B.mods.size());
  for(long m: B.mods){
    auto qr = quad_residues_mod(m);
    std::vector<char> okb(m,0), oka(m,0);
    long nm = (long)umod(n,m);
    for(long r=0;r<m;++r){
      long a2 = (nm + (1LL*r*r)%m) % m; if(qr[a2]) okb[r]=1;
      long d  = ((1LL*r*r)%m - nm) % m; if(d<0) d+=m; if(qr[d]) oka[r]=1;
    }
    B.filt.push_back({m, std::move(okb), std::move(oka)});
  }

  // Parity: for odd n, b and a have opposite parity.
  int n4 = (int)umod(n,4);
  for(long r=0;r<B.M;++r){
    const bool parity_b_ok = ( (n4==1 && !(r&1)) || (n4==3 &&  (r&1)) ); // n≡1: b even; n≡3: b odd
    const bool parity_a_ok = !parity_b_ok;                               // a opposite
    if(parity_b_ok){
      bool ok=true; for(auto& f:B.filt){ if(!f.ok_b[r%f.m]){ ok=false; break; } }
      if(ok) B.allow_b.push_back(r);
    }
    if(parity_a_ok){
      bool ok=true; for(auto& f:B.filt){ if(!f.ok_a[r%f.m]){ ok=false; break; } }
      if(ok) B.allow_a.push_back(r);
    }
  }
  return B;
}

struct QRPrimes {
  std::vector<int> primes;
  OnFlyQR extra;
};

static QRPrimes build_qr_onfly(const mpz_class& n, const BuildCRT& B){
  QRPrimes Q;
  std::vector<int> qr_primes = {
    17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
    191, 193, 197, 199, 211
  };
  int log2_n = mpz_sizeinbase(n.get_mpz_t(), 2);
  int max_qr_prime = std::min(211, log2_n);
  for(int p: qr_primes){
    if(p<=max_qr_prime && std::find(B.mods.begin(), B.mods.end(), p)==B.mods.end())
      Q.primes.push_back(p);
  }
  Q.extra.init(Q.primes);
  return Q;
}

// Vertical spacing and stream selectors
enum VSpaceMode { VSPACE_UNIFORM=1, VSPACE_EXPONENTIAL=2 };
enum StreamMode { STREAM_DUAL=0, STREAM_A=1, STREAM_B=2 };

static FermatStats fermat_dual_scan(
  const mpz_class& n,
  VSpaceMode vspace_mode,
  StreamMode stream_mode,
  bool use_gpu_sieve // affects only b-stream when CUDA build is enabled
){
  FermatStats S; auto t0=std::chrono::high_resolution_clock::now();

  BuildCRT B = build_crt_filters(n);
  S.M = B.M; S.allow_b = B.allow_b.size(); S.allow_a = B.allow_a.size();

  // Bounds
  const mpz_class bmax_mpz = isqrt(n);
  const uint64_t  bmax     = mpz_to_u64(bmax_mpz);
  mpz_class a_base = isqrt(n); if(a_base*a_base < n) a_base += 1; // ceil(sqrt(n))

  // Ensure a_base has the required parity (opposite of b)
  {
    int n4 = (int)umod(n,4);
    bool a_should_be_odd = (n4 == 1); // n≡1 ⇒ a odd; n≡3 ⇒ a even
    bool a_is_odd = mpz_odd_p(a_base.get_mpz_t()) != 0;
    if (a_should_be_odd != a_is_odd) a_base += 1;
  }

  const mpz_class amax_mpz = a_base + bmax_mpz; // vertical window up to sqrt(n)+bmax
  const uint64_t  amax_u   = mpz_to_u64(amax_mpz);
  const uint64_t  M_u      = (uint64_t)B.M;

  QRPrimes Q = build_qr_onfly(n, B);

  const uint64_t STATUS_EVERY = 10'000'000ULL;
  const uint64_t TARGET_CANDIDATES_PER_BAND = 50'000'000ULL;

  const uint64_t Ab = (uint64_t)B.allow_b.size();
  const uint64_t Aa = (uint64_t)B.allow_a.size();

  auto band_size_uniform = [&](uint64_t A)->uint64_t{
    if(A==0) return 0;
    uint64_t K = (TARGET_CANDIDATES_PER_BAND + A - 1) / A;
    return (K==0)?1:K;
  };

  std::atomic<bool> any_found{false};
  Winner best; uint64_t printed_b_at=0, printed_a_at=0;

  const bool run_b = (stream_mode != STREAM_A) && (Ab>0);
  const bool run_a = (stream_mode != STREAM_B) && (Aa>0);

#ifdef USE_CUDA
  // Prepare GPU tables once (only if we'll run b-stream with GPU)
  std::vector<uint8_t>  QR_flat;
  std::vector<uint32_t> P_off;
  std::vector<uint32_t> n_mod_p, M_mod_p, r_mod_p_flat, p_off_flat;

  auto build_gpu_tables_for_b = [&](void){
    const int npr = (int)Q.primes.size();
    P_off.resize(npr);
    n_mod_p.resize(npr);
    M_mod_p.resize(npr);
    p_off_flat.resize(npr+1);
    p_off_flat[0]=0;
    QR_flat.clear();
    size_t sum_off = 0;

    for(int j=0;j<npr;++j){
      int p = Q.primes[j];
      auto qr = quad_residues_mod(p);
      P_off[j] = (uint32_t)sum_off;
      for(int x=0;x<p;++x) QR_flat.push_back(qr[x]?1:0);
      sum_off += p;
      n_mod_p[j] = (uint32_t)umod(n, (unsigned long)p);
      M_mod_p[j] = (uint32_t)(B.M % p);
      p_off_flat[j+1] = p_off_flat[j] + (uint32_t)B.allow_b.size();
    }
    // r_mod_p_flat: for each prime, all (r % p) for r in allow_b
    r_mod_p_flat.resize(p_off_flat.back());
    for(size_t j=0;j<Q.primes.size();++j){
      int p = Q.primes[j];
      uint32_t base = p_off_flat[j];
      for(size_t i=0;i<B.allow_b.size();++i){
        r_mod_p_flat[base + (uint32_t)i] = (uint32_t)(B.allow_b[i] % p);
      }
    }
  };
  if(use_gpu_sieve && run_b && !Q.primes.empty()) build_gpu_tables_for_b();
#endif

#ifdef _OPENMP
  #pragma omp parallel
  #pragma omp single nowait
#endif
  {
    // ===== b-stream task =====
    if(run_b){
#ifdef _OPENMP
      #pragma omp task shared(any_found, best, printed_b_at)
#endif
      {
        const uint64_t kmax_b = mpz_to_u64(bmax_mpz / B.M);
        uint64_t k0 = 0;
        while(k0 <= kmax_b && !any_found.load(std::memory_order_relaxed)){
          uint64_t Kpb_b = band_size_uniform(Ab);
          uint64_t k1 = k0 + Kpb_b - 1; if(k1 > kmax_b) k1 = kmax_b;

#ifdef USE_CUDA
          if(use_gpu_sieve && !Q.primes.empty()){
            auto GR = gpu_sieve_b_band(
              k0, k1, M_u, bmax,
              Q.primes, QR_flat, P_off, n_mod_p, M_mod_p,
              B.allow_b, r_mod_p_flat, p_off_flat
            );
            for(uint64_t b_val : GR.b_hits){
              if(any_found.load(std::memory_order_relaxed)) break;
              S.candidates_b += 1; // count verified attempts
              mpz_class b = mpz_from_u64(b_val);
              if(!Q.extra.pass_b(n, b)) continue; // usually redundant
              mpz_class a2 = n + b*b;
              if(!is_square(a2)) continue;
              mpz_class a = isqrt(a2);
              mpz_class p = a - b, q = a + b;
              if(p<=0 || q<=0) continue;
              Winner w; w.valid=true; w.b=b; w.a=a; w.p=p; w.q=q; w.r=(long)(b_val%M_u); w.stream=0;
#ifdef _OPENMP
              #pragma omp critical
#endif
              {
                if(!any_found.load(std::memory_order_relaxed) || winner_less(w, best)){
                  best=w; any_found.store(true, std::memory_order_release);
                }
              }
            }
          } else
#endif
          {
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
              Winner local; bool have=false;
#ifdef _OPENMP
              #pragma omp for schedule(static)
#endif
              for(size_t i=0;i<B.allow_b.size();++i){
                if(any_found.load(std::memory_order_relaxed)) continue;
                long r = B.allow_b[i];
                for(uint64_t k=k0; k<=k1 && !any_found.load(std::memory_order_relaxed); ++k){
                  mpz_class b = mpz_from_u64(k) * B.M + r;
                  if(b > bmax_mpz) break;
#ifdef _OPENMP
                  #pragma omp atomic update
#endif
                  S.candidates_b += 1;

                  if((S.candidates_b - printed_b_at) >= STATUS_EVERY){
#ifdef _OPENMP
                    #pragma omp critical
#endif
                    {
                      if((S.candidates_b - printed_b_at) >= STATUS_EVERY){
                        printed_b_at = S.candidates_b;
                        double secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t0).count();
                        std::cerr << "[SF-b] tested="<< S.candidates_b <<" allow="<< S.allow_b
                                  <<" stride(M)="<< S.M <<" k="<< k <<" elapsed="<< secs <<"s\n";
                      }
                    }
                  }

                  if(!Q.extra.pass_b(n, b)) continue;
                  mpz_class a2 = n + b*b;
                  if(!is_square(a2)) continue;

                  mpz_class a = isqrt(a2);
                  mpz_class p = a - b, q = a + b;
                  if(p <= 0 || q <= 0) continue;

                  Winner w; w.valid=true; w.b=b; w.a=a; w.p=p; w.q=q; w.r=r; w.stream=0;
                  if(!have || winner_less(w, local)){ local=w; have=true; }
                  break; // smallest b for this (i,k)
                }
              }
              if(have){
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                  if(!any_found.load(std::memory_order_relaxed) || winner_less(local, best)){
                    best = local; any_found.store(true, std::memory_order_release);
                  }
                }
              }
            } // end parallel
          }

          if(any_found.load(std::memory_order_acquire)) break;
          k0 = (k1 == UINT64_MAX ? k1 : k1 + 1);
        } // while bands
      }
    } // run_b

    // ===== a-stream task (vertical) =====
    if(run_a){
#ifdef _OPENMP
      #pragma omp task shared(any_found, best, printed_a_at)
#endif
      {
        auto do_band = [&](uint64_t k0, uint64_t k1){
#ifdef _OPENMP
          #pragma omp parallel
#endif
          {
            Winner local; bool have=false;
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for(size_t i=0;i<B.allow_a.size();++i){
              if(any_found.load(std::memory_order_relaxed)) continue;
              long r = B.allow_a[i];

              mpz_class r_mpz = (unsigned long)r;  // r ∈ [0,M)

              // Clamp k so that a = r + kM ∈ [a_base, amax]
              uint64_t k_start = 0;
              if (a_base > r_mpz) {
                mpz_class diff = a_base - r_mpz;
                k_start = ceil_div_mpz_ui(diff, (unsigned long)B.M);
              }
              uint64_t k_stop = 0;
              if (amax_mpz > r_mpz) {
                mpz_class diff2 = amax_mpz - r_mpz;
                k_stop = floor_div_mpz_ui(diff2, (unsigned long)B.M);
              } else {
                continue; // cannot reach a >= a_base for this residue
              }

              uint64_t kb = (k0 > k_start ? k0 : k_start);
              uint64_t ke = (k1 < k_stop  ? k1 : k_stop);
              if (kb > ke) continue;

              for(uint64_t k=kb; k<=ke && !any_found.load(std::memory_order_relaxed); ++k){
                mpz_class a = mpz_from_u64(k) * B.M + r;

#ifdef _OPENMP
                #pragma omp atomic update
#endif
                S.candidates_a += 1;

                if((S.candidates_a - printed_a_at) >= STATUS_EVERY){
#ifdef _OPENMP
                  #pragma omp critical
#endif
                  {
                    if((S.candidates_a - printed_a_at) >= STATUS_EVERY){
                      printed_a_at = S.candidates_a;
                      double secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t0).count();
                      std::cerr << "[SF-a] tested="<< S.candidates_a <<" allow="<< S.allow_a
                                <<" stride(M)="<< S.M <<" k="<< k <<" elapsed="<< secs <<"s\n";
                    }
                  }
                }

                if(!Q.extra.pass_a(n, a)) continue;
                mpz_class delta = a*a - n;
                if(!is_square(delta)) continue;

                mpz_class b = isqrt(delta);
                mpz_class p = a - b, q = a + b;
                if(p <= 0 || q <= 0) continue;

                Winner w; w.valid=true; w.a=a; w.b=b; w.p=p; w.q=q; w.r=r; w.stream=1;
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                  if(!any_found.load(std::memory_order_relaxed) || winner_less(w, best)){
                    best = w; any_found.store(true, std::memory_order_release);
                  }
                }
                break; // smallest b for this (i,k) via smallest a
              }
            }
          }; // end parallel
        };

        const uint64_t kmax_a = (amax_u / (uint64_t)B.M);
        if(vspace_mode == VSPACE_UNIFORM){
          uint64_t Kpb_a = band_size_uniform(Aa);
          for(uint64_t k0=0; k0<=kmax_a && !any_found.load(std::memory_order_relaxed); ){
            uint64_t k1 = k0 + Kpb_a - 1; if(k1 > kmax_a) k1 = kmax_a;
            do_band(k0, k1);
            if(any_found.load(std::memory_order_acquire)) break;
            k0 = (k1 == UINT64_MAX ? k1 : k1 + 1);
          }
        }else{ // VSPACE_EXPONENTIAL
          uint64_t k0=0, span=1;
          while(k0<=kmax_a && !any_found.load(std::memory_order_relaxed)){
            uint64_t k1 = k0 + span - 1; if(k1 > kmax_a) k1 = kmax_a;
            do_band(k0, k1);
            if(any_found.load(std::memory_order_acquire)) break;
            k0 = (k1 == UINT64_MAX ? k1 : k1 + 1);
            if(span < (1u<<28)) span <<= 1; // grow cautiously
          }
        }
      }
    } // run_a
  } // end omp single/region

  if(any_found.load(std::memory_order_acquire)){
    S.found = true; S.a=best.a; S.b=best.b; S.p=best.p; S.q=best.q; S.k2b = 2*best.b; S.stream = best.stream;
  }
  S.secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t0).count();
  return S;
}

// ===== Main =====
int main(int argc, char** argv){
  bool use_gpu = false;
  VSpaceMode vspace = VSPACE_UNIFORM;
  StreamMode stream = STREAM_DUAL;

  for(int i=1;i<argc;++i){
    if(std::strcmp(argv[i],"--use-gpu")==0) use_gpu=true;
    else if(std::strncmp(argv[i],"--vspace=",9)==0){
      int m = std::atoi(argv[i]+9);
      vspace = (m==2)?VSPACE_EXPONENTIAL:VSPACE_UNIFORM;
    }else if(std::strncmp(argv[i],"--stream=",9)==0){
      const char* v = argv[i]+9;
      if(std::strcmp(v,"a")==0) stream=STREAM_A;
      else if(std::strcmp(v,"b")==0) stream=STREAM_B;
      else stream=STREAM_DUAL;
    }
  }
#ifndef USE_CUDA
  if(use_gpu){
    std::cerr << "[warn] --use-gpu requested but binary not built with -DUSE_CUDA; using CPU path.\n";
    use_gpu = false;
  }
#endif

  std::cout << "Enter n (integer): ";
  std::string s; if(!(std::cin>>s)){ std::cerr<<"Input error\n"; return 1; }
  mpz_class n; if(mpz_set_str(n.get_mpz_t(), s.c_str(), 10)!=0){ std::cerr<<"Invalid n\n"; return 1; }

  if(mpz_even_p(n.get_mpz_t())){ mpz_class half=n/2;
    std::cout << "\n=== Result (even n) ===\n"
              << "p = 2\nq = " << half << "\n"; return 0; }

  // Stage 0: TD
  auto td = tiny_td(n);
  if(td.found){ if(td.p>td.q) std::swap(td.p,td.q);
    std::cout << "\n=== Result (tiny TD) ===\n"
              << "p = " << td.p << "\nq = " << td.q << "\nVerified: " << td.p*td.q
              << "\nstats: candidates="<<td.work<<"  elapsed="<<td.secs<<"s\n"; return 0; }
  std::cout << "\n[tiny TD] no hit  (candidates="<<td.work<<", elapsed="<<td.secs<<"s)\n";

  // Stage 1: IA
  auto ia = ia_probe(n);
  if(ia.found){ if(ia.p>ia.q) std::swap(ia.p,ia.q);
    std::cout << "\n=== Result (IA) ===\n"
              << "p = " << ia.p << "\nq = " << ia.q << "\nVerified: " << ia.p*ia.q
              << "\nstats: gcd_checks="<<ia.work<<"  elapsed="<<ia.secs<<"s\n"; return 0; }
  std::cout << "[IA] no hit  (gcd_checks="<<ia.work<<", elapsed="<<ia.secs<<"s)\n";

  // Stage 2: Dual Sieved Fermat (configurable streams + vertical spacing)
  auto fm = fermat_dual_scan(n, vspace, stream, use_gpu);
  std::cout << "\n=== Result (sieved Fermat) ===\n";
  if(fm.found){
    mpz_class p=fm.p, q=fm.q; if(p>q) std::swap(p,q);
    std::cout << "(q - p)=2b="<< fm.k2b << "\n"
              << "a="<< fm.a <<"  b="<< fm.b << "  stream=" << (fm.stream==0 ? "b" : "a") << "\n"
              << "p="<< p <<"  q="<< q << "\n"
              << "Verified: " << p*q << "\n"
              << "stats: cand_b="<< fm.candidates_b
              << "  cand_a="<< fm.candidates_a
              << "  allow_b="<< fm.allow_b
              << "  allow_a="<< fm.allow_a
              << "  stride(M)="<< fm.M
              << "  elapsed="<< fm.secs <<"s\n";
  }else{
    std::cout << "No hit. cand_b="<< fm.candidates_b
              << "  cand_a="<< fm.candidates_a
              << "  allow_b="<< fm.allow_b
              << "  allow_a="<< fm.allow_a
              << "  stride(M)="<< fm.M
              << "  elapsed="<< fm.secs <<"s\n";
  }
  return 0;
}
