#!/usr/bin/env python3
"""Reproducible, exact-integer tests of real-pair/Hadamard bit encoding for factor screening.

Source: Viole (2026), Real-Pair Representations of Complex Numbers and Quaternions.
Experiment: extend T=[[1,-1],[1,1]] to T^(tensor k), then test invertibility,
transported integer multiplication with carries, XOR-vs-positional convolution,
and candidate factor screening against directly comparable modular arithmetic.

No claims of asymptotic improvement are assumed. Run: python3 real_pair_factor_test.py
"""
import csv
import json
import math
import random
import statistics
import time
from pathlib import Path

SEED = 20261118
RNG = random.Random(SEED)
HERE = Path(__file__).resolve().parent


def bits(x: int, width: int) -> list[int]:
    assert 0 <= x < 1 << width
    return [(x >> j) & 1 for j in range(width)]


def bitvalue(a) -> int:
    return sum(int(v) << j for j, v in enumerate(a))


def realpair_hadamard(vector) -> tuple[int, ...]:
    """Exact T^(tensor k) where T(x,y)=(x-y,x+y) in butterfly order."""
    a = list(vector)
    n = len(a)
    assert n > 0 and n & (n - 1) == 0
    half = 1
    while half < n:
        stride = 2 * half
        for lo in range(0, n, stride):
            for j in range(lo, lo + half):
                x, y = a[j], a[j + half]
                a[j] = x - y
                a[j + half] = x + y
        half *= 2
    return tuple(a)


def realpair_inverse(vector) -> tuple[int, ...]:
    """Exact inverse of realpair_hadamard for vectors in the integer image."""
    a = list(vector)
    n = len(a)
    assert n > 0 and n & (n - 1) == 0
    half = n // 2
    while half:
        stride = 2 * half
        for lo in range(0, n, stride):
            for j in range(lo, lo + half):
                x, y = a[j], a[j + half]
                assert (x + y) % 2 == 0 and (y - x) % 2 == 0
                a[j] = (x + y) // 2
                a[j + half] = (y - x) // 2
        half //= 2
    return tuple(a)


def coeff_conv_carry(a, b):
    """Polynomial convolution of bit coefficients followed by base-2 carries."""
    n, m = len(a), len(b)
    out = [0] * (n + m)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i+j] += ai * bj
    carry = 0
    for k in range(len(out)):
        v = out[k] + carry
        out[k] = v & 1
        carry = v >> 1
    assert carry == 0
    return out


def realpair_mult_complete(a: int, b: int, width: int) -> int:
    """Full encode -> exact transported product -> decode, no native multiplication."""
    A = realpair_hadamard(bits(a, width))
    B = realpair_hadamard(bits(b, width))
    product_bits = coeff_conv_carry(realpair_inverse(A), realpair_inverse(B))
    result = realpair_hadamard(product_bits)
    return bitvalue(realpair_inverse(result))


def realpair_mult_core(A, B):
    """Product of already encoded operands; returns an encoded product."""
    product_bits = coeff_conv_carry(realpair_inverse(A), realpair_inverse(B))
    return realpair_hadamard(product_bits)


def schoolbook(a, b, width):
    return bitvalue(coeff_conv_carry(bits(a, width), bits(b, width)))


def canonical_fwht(a):
    """Standard Walsh-Hadamard H(x,y)=(x+y,x-y), diagonalizing XOR convolution."""
    out = list(a)
    step = 1
    while step < len(out):
        for off in range(0, len(out), 2*step):
            for j in range(off, off+step):
                x, y = out[j], out[j+step]
                out[j] = x+y
                out[j+step] = x-y
        step *= 2
    return out


def xor_convolution(a, b, width):
    A = canonical_fwht(bits(a, width))
    B = canonical_fwht(bits(b, width))
    C = canonical_fwht([x*y for x, y in zip(A, B)])
    assert all(z % width == 0 for z in C)
    return [z // width for z in C]


def is_prime(n):
    """Deterministic Miller-Rabin for n < 2^64."""
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        s += 1; d //= 2
    for a in (2, 325, 9375, 28178, 450775, 9780504, 179526502):
        a %= n
        if a == 0: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s-1):
            x = x*x % n
            if x == n - 1: break
        else:
            return False
    return True


def next_prime(n):
    k = n | 1
    while not is_prime(k): k += 2
    return k


def time_call(fn, repeat=1, reps=5):
    times = []
    for _ in range(reps):
        t0 = time.perf_counter_ns()
        for _ in range(repeat): fn()
        times.append((time.perf_counter_ns()-t0)/repeat)
    return statistics.median(times)


def quaternion_realpair(q):
    a,b,c,d=q
    return (a-b,a+b,c-d,c+d)


def quaternion_realpair_inverse(v):
    x,y,u,w=v
    assert (x+y)%2==0 and (y-x)%2==0 and (u+w)%2==0 and (w-u)%2==0
    return ((x+y)//2, (y-x)//2, (u+w)//2, (w-u)//2)


def quaternion_multiply(q,r):
    a,b,c,d=q
    e,f,g,h=r
    return (a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)


def complex_realpair_product(v, w):
    # Equation (3) of the source paper; all intermediate halves are integral
    # on the parity lattice that is the image of the Gaussian integers.
    x,y=v
    u,z=w
    a=x*z+y*u-y*z+x*u
    b=x*z+y*u+y*z-x*u
    assert a%2==0 and b%2==0
    return (a//2,b//2)


def quaternion_realpair_product(v,w):
    # Transport equation (14): (z1*w1-z2*conj(w2), z1*w2+z2*conj(w1)).
    z1,z2=v[:2],v[2:]
    w1,w2=w[:2],w[2:]
    cm=lambda z: (z[1],z[0])
    a=complex_realpair_product(z1,w1)
    b=complex_realpair_product(z2,cm(w2))
    c=complex_realpair_product(z1,w2)
    d=complex_realpair_product(z2,cm(w1))
    return (a[0]-b[0],a[1]-b[1],c[0]+d[0],c[1]+d[1])


def spectral_mod3(d,width):
    # T tensor row (width-2) equals sum_j (-1)^j b_j.
    # Because 2 = -1 (mod 3), this is exactly d modulo 3.
    return realpair_hadamard(bits(d,width))[width-2] % 3


def check_encoding_and_multiply():
    rows=[]
    for width in [2, 4, 8, 16, 32, 64, 128, 256]:
        count = 160 if width <= 32 else 48
        for _ in range(count):
            a = RNG.randrange(1 << width)
            b = RNG.randrange(1 << width)
            A = realpair_hadamard(bits(a, width))
            B = realpair_hadamard(bits(b, width))
            assert realpair_inverse(A) == tuple(bits(a, width))
            assert realpair_inverse(B) == tuple(bits(b, width))
            assert sum(x*x for x in A) == width * a.bit_count()
            assert realpair_mult_complete(a, b, width) == a*b
        rows.append({'operand_bits':width,'random_pairs':count,'inverse_fails':0,'product_fails':0,'norm_fails':0})
    # The quaternion coordinate product is exact but not binary integer multiplication.
    assert quaternion_multiply((0,1,0,0), (0,1,0,0)) == (-1,0,0,0)
    assert 2*2 == 4
    for _ in range(1000):
        qa=tuple(RNG.randrange(-3,4) for _ in range(4))
        qb=tuple(RNG.randrange(-3,4) for _ in range(4))
        pa=quaternion_realpair(qa)
        pb=quaternion_realpair(qb)
        assert quaternion_realpair_product(pa,pb) == quaternion_realpair(quaternion_multiply(qa,qb))
        assert quaternion_realpair_inverse(pa)==qa
        assert quaternion_realpair_inverse(pb)==qb
    return rows


def check_xor_is_not_integer_mult():
    rows=[]
    for width, count in ((4,256), (8,4096), (16,1000)):
        matches = 0
        for idx in range(count):
            if width==4:
                a,b=divmod(idx,16)
            else:
                a=RNG.randrange(1,1<<width);b=RNG.randrange(1,1<<width)
            xor=bitvalue(xor_convolution(a,b,width))
            matches += xor==a*b
        rows.append({'bits':width, 'pairs':count, 'matched_integer_product':matches,
                     'failed_integer_product':count-matches})
    assert bitvalue(xor_convolution(2,2,4)) == 1
    assert 2*2 == 4
    return rows


def factor_cases():
    # All semiprimes are generated from close primes deliberately so a tight factor
    # strip can be supplied *externally*, making screening power testable.
    cases=[]
    for prime_bits in (12,16,20,24,28,31):
        for rep in range(6):
            p = next_prime((1 << (prime_bits - 1)) + RNG.randrange(200, min(1 << (prime_bits-2), 10000)))
            q = next_prime(p + RNG.randrange(50, 600))
            n = p*q
            assert p < q and is_prime(p) and is_prime(q) and n < (1<<64)
            assert p.bit_length() == prime_bits and q.bit_length() == prime_bits
            root = math.isqrt(n)
            assert root - 512 <= p <= root and root+1 <= q <= root+512
            cases.append({'prime_bits':prime_bits,'p':p,'q':q,'n':n,'root':root})
    return cases


def screen_cases(cases):
    rows=[]
    total_candidate=0
    for c in cases:
        n, root, p, q = c['n'], c['root'], c['p'], c['q']
        ds=list(range(max(3,root-512)|1,root+1,2))
        assert p in ds and n%p==0
        total_candidate+=len(ds)
        # The alternating-bit Walsh coefficient equals the residue modulo 3.
        m_width = 1 << ((n.bit_length()-1).bit_length())
        assert spectral_mod3(n,m_width) == n%3
        # Real factors also have the correct coefficient modulo 3.
        assert spectral_mod3(p,m_width) == p%3
        qmin,qmax=root+1,root+512
        for t in (8, 12, 16):
            modulus=1<<t
            width=1 << ((t-1).bit_length())  # 8 or 16 dimensional bit vector
            residues = {z%modulus for z in range(qmin,qmax+1)}
            spectral = {realpair_hadamard(bits(z, width)) for z in residues}
            def get_res(d):
                return (n*pow(d,-1,modulus)) % modulus
            qres=[get_res(d) for d in ds]
            direct_filter=lambda: [d for d, z in zip(ds,qres) if z in residues]
            spectral_filter=lambda: [d for d, z in zip(ds,qres)
                                     if realpair_hadamard(bits(z,width)) in spectral]
            direct_survive=direct_filter()
            spectral_survive=spectral_filter()
            assert direct_survive==spectral_survive and p in direct_survive
            # A no-cofactor-bound 2-adic filter cannot reject ANY odd candidate.
            assert all(((d*z-n)%modulus)==0 for d,z in zip(ds,qres))
            hits=[d for d in ds if n%d==0]
            assert hits==[p]
            rows.append({**c,'mod_bits':t,'spectral_dimensions':width,
                         'window_factor_low':root-512,'window_factor_high':root,
                         'window_cofactor_low':qmin,'window_cofactor_high':qmax,
                         'candidates':len(ds),'divisors_found':len(hits),
                         'mod_survivors':len(direct_survive),
                         'spectral_survivors':len(spectral_survive),
                         'false_negative':0})
    return rows


def bench_factor(cases, t):
    modulus=1<<t
    width=1 << ((t-1).bit_length())
    prepared=[]
    for c in cases:
        n, root = c['n'], c['root']
        ds = list(range(max(3,root-512)|1,root+1,2))
        rset={z%modulus for z in range(root+1,root+513)}
        hset={realpair_hadamard(bits(z,width)) for z in rset}
        prepared.append((n, ds,rset,hset))
    allN=sum(len(a[1]) for a in prepared)
    # The q-residue is computed inside each filter, fairly; both filters use exactly
    # the same modular inverse, one using transformed coordinates for membership.
    def trial():
        hit=[]
        for n, ds,_,_ in prepared:
            for d in ds:
                if n%d==0: hit.append(d)
        return hit
    def mod_screen():
        hit=[]
        for n, ds,rset,_ in prepared:
            for d in ds:
                z=(n*pow(d,-1,modulus))%modulus
                if z in rset and n%d==0: hit.append(d)
        return hit
    def spec_screen():
        hit=[]
        for n, ds,_,hset in prepared:
            for d in ds:
                z=(n*pow(d,-1,modulus))%modulus
                if realpair_hadamard(bits(z,width)) in hset and n%d==0: hit.append(d)
        return hit
    def smallprime_direct():
        hit=[]
        for n,ds,_,_ in prepared:
            for d in ds:
                if (n%3!=0 and d%3==0): continue
                if n%d==0:hit.append(d)
        return hit
    def smallprime_spectral():
        hit=[]
        for n,ds,_,_ in prepared:
            w=1 << ((max(n.bit_length(),3)-1).bit_length())
            n3=spectral_mod3(n,w)
            for d in ds:
                if n3 and spectral_mod3(d,w)==0:continue
                if n%d==0:hit.append(d)
        return hit
    expected=trial()
    assert mod_screen()==expected and spec_screen()==expected
    assert smallprime_direct()==expected and smallprime_spectral()==expected
    trials=time_call(trial,reps=7)
    mod3=time_call(smallprime_direct,reps=7)
    spec3=time_call(smallprime_spectral,reps=7)
    mod=time_call(mod_screen,reps=7)
    spec=time_call(spec_screen,reps=7)
    return {'mod_bits':t,'cases':len(cases),'candidates':allN,
            'trial_ms':round(trials/1e6,4),
            'mod3_ms':round(mod3/1e6,4),'spectral_mod3_ms':round(spec3/1e6,4),
            'modular_ms':round(mod/1e6,4),
            'spectral_ms':round(spec/1e6,4),
            'modular_vs_trial':round(mod/trials,2),
            'mod3_vs_trial':round(mod3/trials,2),'spectral_mod3_vs_trial':round(spec3/trials,2),
            'spectral_vs_trial':round(spec/trials,2),
            'spectral_vs_modular':round(spec/mod,2),
            'set_construction_excluded':True}


def bench_multiplication():
    rows=[]
    for width in [8,16,32,64,128,256]:
        pairs=[]
        for _ in range(12):
            a=(1<<(width-1))|RNG.getrandbits(width-1)
            b=(1<<(width-1))|RNG.getrandbits(width-1)
            pairs.append((a,b))
        encoded=[(realpair_hadamard(bits(a,width)),realpair_hadamard(bits(b,width))) for a,b in pairs]
        def native():
            return [a*b for a,b in pairs]
        def textbook():
            return [schoolbook(a,b,width) for a,b in pairs]
        def transported():
            return [realpair_mult_complete(a,b,width) for a,b in pairs]
        def core():
            return [realpair_mult_core(a,b) for a,b in encoded]
        assert textbook()==native()==transported()
        assert [bitvalue(realpair_inverse(v)) for v in core()]==native()
        # Amortize timer overhead for native work but avoid large slow cases.
        rep_native=1000 if width<=32 else 400
        nat=time_call(native,repeat=rep_native,reps=7)
        sb=time_call(textbook,repeat=3,reps=5)
        tr=time_call(transported,repeat=2,reps=5)
        co=time_call(core,repeat=2,reps=5)
        rows.append({'operand_bits':width,'pairs_per_batch':len(pairs),
                     'native_microseconds_per_pair':round(nat/len(pairs)/1000,4),
                     'schoolbook_microseconds_per_pair':round(sb/len(pairs)/1000,4),
                     'transported_full_microseconds_per_pair':round(tr/len(pairs)/1000,4),
                     'transported_core_microseconds_per_pair':round(co/len(pairs)/1000,4),
                     'transported_full_vs_native':round(tr/nat,1),
                     'transported_full_vs_schoolbook':round(tr/sb,2)})
    return rows


def fermat_iterations(cases):
    vals=[]
    for c in cases:
        n=c['n'];a=math.isqrt(n)
        if a*a<n:a+=1
        t=0
        while True:
            b2=a*a-n
            b=math.isqrt(b2)
            t+=1
            if b*b==b2:
                assert (a-b,a+b)==(c['p'],c['q'])
                break
            a+=1
        vals.append(t)
    return {'min':min(vals),'median':statistics.median(vals),'max':max(vals),
            'total':sum(vals)}


def write_csv(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)


def main():
    t0=time.perf_counter()
    checks=check_encoding_and_multiply()
    xor=check_xor_is_not_integer_mult()
    cases=factor_cases()
    screens=screen_cases(cases)
    fbench=[bench_factor(cases,t) for t in (8,12,16)]
    multiplication=bench_multiplication()
    fermat=fermat_iterations(cases)
    write_csv('real_pair_factor_cases.csv', screens)
    write_csv('real_pair_multiplication_benchmarks.csv', multiplication)
    aggregates=[]
    for t in (8,12,16):
        a=[r for r in screens if r['mod_bits']==t]
        candidates=sum(v['candidates'] for v in a)
        survivors=sum(v['mod_survivors'] for v in a)
        aggregates.append({'mod_bits':t,'cases':len(a),'candidates':candidates,
                           'survivors':survivors,'rejected_pct':round((1-survivors/candidates)*100,3),
                           'missed_factors':sum(v['false_negative'] for v in a)})
    result={'seed':SEED,'source_paper':'Real-Pair Representations of Complex Numbers and Quaternions, 2026-11-04',
            'transform':'T=[[1,-1],[1,1]], T tensor power, exact integer coefficients',
            'checks':checks,'xor':xor,'semiprimes':len(cases),
            'prime_bits':sorted({c['prime_bits'] for c in cases}),
            'experimental_geometry':'Near-square semiprimes; a provided +/-512 cofactor/factor strip (external prior); odd candidates only',
            'factor_summary':aggregates,'factor_benchmarks':fbench,
            'fermat_iterations':fermat,'multiplication_benchmarks':multiplication,
            'duration_s':round(time.perf_counter()-t0,2)}
    (HERE/'real_pair_factor_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()