#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>

// ---------- Precision helpers ----------
static void set_working_prec_bits(const mpz_class& n) {
    // Bits ~ 4x number of bits of n, min 512, cap at something sensible
    size_t bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    size_t prec = std::max<size_t>(512, 4 * bits);
    mpf_set_default_prec(prec);
}

// ---------- Integer perfect square checks ----------
static bool is_perfect_square(const mpz_class& x, mpz_class* root = nullptr) {
    if (x < 0) return false;
    mpz_class r;
    mpz_sqrt(r.get_mpz_t(), x.get_mpz_t());
    if (r * r == x) { if (root) *root = r; return true; }
    return false;
}

// Newton integer sqrt for mpz (thin wrapper around mpz_sqrt)
static mpz_class isqrt(const mpz_class& n) {
    mpz_class r;
    mpz_sqrt(r.get_mpz_t(), n.get_mpz_t());
    return r;
}

// ---------- Gaussian-square test for z = x + i y ----------
// Returns true and fills a,b if (a+bi)^2 == x + i y  (a,b integers)
static bool gaussian_integer_sqrt(const mpz_class& x, const mpz_class& y, mpz_class& a, mpz_class& b) {
    // 1) m^2 = x^2 + y^2 must be a perfect square
    mpz_class m2 = x*x + y*y;
    mpz_class m;
    if (!is_perfect_square(m2, &m)) return false;

    // 2) U=(m+x)/2 and V=(m-x)/2 must be nonneg squares
    mpz_class mx = m + x;
    mpz_class mx2 = m - x;
    if ((mx & 1) != 0 || (mx2 & 1) != 0) return false; // parity
    mpz_class U = mx >> 1;
    mpz_class V = mx2 >> 1;
    mpz_class ra, rb;
    if (!is_perfect_square(U, &ra)) return false;
    if (!is_perfect_square(V, &rb)) return false;

    // Signs so that 2ab == y
    // Try (ra, rb) and flip one sign if needed
    mpz_class twoab = 2 * ra * rb;
    if (twoab == y) { a = ra; b = rb; return true; }
    if (-twoab == y) { a = ra; b = -rb; return true; }

    // If y==0 and either ra==0 or rb==0, both sign choices work
    if (y == 0) {
        a = ra;
        b = mpz_class(0);
        return true;
    }

    return false;
}

// ---------- Big complex type over mpf_class ----------
struct BigC {
    mpf_class re, im;

    BigC() : re(0), im(0) {}
    BigC(const mpf_class& r, const mpf_class& i) : re(r), im(i) {}
};

static inline BigC add(const BigC& a, const BigC& b) { return BigC(a.re + b.re, a.im + b.im); }
static inline BigC subc(const BigC& a, const BigC& b) { return BigC(a.re - b.re, a.im - b.im); }
static inline BigC mul(const BigC& a, const BigC& b) {
    return BigC(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}
static inline BigC scale(const BigC& a, const mpf_class& s) { return BigC(a.re*s, a.im*s); }
static inline mpf_class norm2(const BigC& a) { return a.re*a.re + a.im*a.im; }
static inline BigC conjc(const BigC& a) { return BigC(a.re, -a.im); }

// a / b
static inline BigC divc(const BigC& a, const BigC& b) {
    mpf_class d = norm2(b);
    BigC num = mul(a, conjc(b));
    return BigC(num.re / d, num.im / d);
}

// |a| (Euclidean)
static inline mpf_class abs_c(const BigC& a) {
    mpf_class n2 = norm2(a);
    // sqrt for mpf_class
    mpf_class r; mpf_sqrt(r.get_mpf_t(), n2.get_mpf_t());
    return r;
}

// ---------- Core: Newton along vertical line ----------
// Iteration:
//  w_{k+1} = 0.5 * ( w_k + (n + i*y_k) / w_k )
//  y_{k+1} = Im( w_{k+1}^2 )
//
// Returns true if it converges to integers (a,b) such that (a+bi)^2 == n + i*(2ab).
static bool vertical_line_newton(const mpz_class& n,
                                 BigC w0,
                                 mpf_class y0,
                                 mpf_class tol,
                                 int max_iter,
                                 mpz_class& a_out,
                                 mpz_class& b_out) {
    BigC w = w0;
    mpf_class y = y0;

    // Constant z's real part:
    mpf_class n_f(n); // mpf copy of n

    for (int k = 0; k < max_iter; ++k) {
        // Compute next Newton step for sqrt(n + i*y)
        BigC z(n_f, y);
        if (w.re == 0 && w.im == 0) {
            // avoid division by zero: nudge
            w.re = n_f > 0 ? n_f : 1;
            w.im = 1;
        }
        BigC w_inv = divc(z, w);
        BigC w_next = scale(add(w, w_inv), mpf_class(0.5));

        // Enforce the vertical-line constraint by setting y_{k+1} = Im(w_next^2)
        BigC w2 = mul(w_next, w_next);
        y = w2.im;

        // Check residual |w_next^2 - (n + i*y)| = |(re - n) + i(im - y)| = |(w2.re - n, 0)|
        mpf_class re_err = w2.re - n_f;
        if (re_err < 0) re_err = -re_err;
        // Im error is zero by construction
        if (re_err < tol) {
            // Round to nearest integers
            mpf_class a_f = w_next.re;
            mpf_class b_f = w_next.im;
            // nearest-integer rounding
            mpz_class a = a_f >= 0 ? mpz_class(a_f + 0.5) : mpz_class(a_f - 0.5);
            mpz_class b = b_f >= 0 ? mpz_class(b_f + 0.5) : mpz_class(b_f - 0.5);

            // Verify exact hit via Gaussian-square test on z = n + i*(2ab)
            mpz_class y_check = 2 * a * b;
            mpz_class aa, bb;
            if (gaussian_integer_sqrt(n, y_check, aa, bb)) {
                a_out = aa; b_out = bb;
                return true;
            }
            // If close but missed (rounding), continue a few more steps
        }

        w = w_next;
    }

    // Try final rounding check anyway
    mpf_class a_f = w.re, b_f = w.im;
    mpz_class a = a_f >= 0 ? mpz_class(a_f + 0.5) : mpz_class(a_f - 0.5);
    mpz_class b = b_f >= 0 ? mpz_class(b_f + 0.5) : mpz_class(b_f - 0.5);
    mpz_class y_check = 2 * a * b;
    mpz_class aa, bb;
    if (gaussian_integer_sqrt(n, y_check, aa, bb)) { a_out = aa; b_out = bb; return true; }

    return false;
}

// ---------- High-level factorization via vertical-line Newton ----------
static bool factor_via_vertical_newton(const mpz_class& n, mpz_class& f1, mpz_class& f2) {
    if (n <= 1) return false;
    if (n % 2 == 0) { f1 = 2; f2 = n/2; return true; }

    set_working_prec_bits(n);

    // Initial guess for w ~ sqrt(n): w0 = (sqrt(n), 1)
    mpf_class n_f(n);
    mpf_class sqrt_n; mpf_sqrt(sqrt_n.get_mpf_t(), n_f.get_mpf_t());

    // tolerances and attempts
    mpf_class tol = 1e-30;
    const int max_iter = 200;

    // Try a few small imaginary seeds around b0 in {1,3,5,7,9}, and y0 around 0
    std::vector<long> seed_b = {1, 3, 5, 7, 9, -1, -3, -5};
    std::vector<long> seed_y = {0, 2, -2, 4, -4};

    for (long b0 : seed_b) {
        for (long y0 : seed_y) {
            BigC w0(sqrt_n, mpf_class(b0));
            mpf_class y_init = y0;
            mpz_class a, b;
            if (vertical_line_newton(n, w0, y_init, tol, max_iter, a, b)) {
                // Found sqrt(n + i*(2ab)) = a + bi  => factors (a-b, a+b)
                mpz_class p = a - b;
                mpz_class q = a + b;
                if (p > 1 && q > 1 && p*q == n) { f1 = p; f2 = q; return true; }
                // If signs swapped, try (b-a, b+a) (shouldn't be needed, but safe)
                mpz_class p2 = b - a;
                mpz_class q2 = b + a;
                if (p2 > 1 && q2 > 1 && p2*q2 == n) { f1 = p2; f2 = q2; return true; }
            }
        }
    }

    // As a last resort, try reading a,b from a direct Gaussian test in a small window around rounded sqrt(n)
    // (not strictly necessary, just a tiny safety-net)
    mpz_class a0 = sqrt_n >= 0 ? mpz_class(sqrt_n + 0.5) : mpz_class(sqrt_n - 0.5);
    for (long db = -5; db <= 5; ++db) {
        for (long da = -5; da <= 5; ++da) {
            mpz_class a = a0 + da;
            mpz_class b = db;
            mpz_class y = 2 * a * b;
            mpz_class aa, bb;
            if (gaussian_integer_sqrt(n, y, aa, bb)) {
                mpz_class p = aa - bb;
                mpz_class q = aa + bb;
                if (p > 1 && q > 1 && p*q == n) { f1 = p; f2 = q; return true; }
            }
        }
    }

    return false;
}

// ---------- CLI ----------
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n>\n";
        return 1;
    }
    mpz_class n(argv[1]);

    mpz_class f1, f2;
    if (factor_via_vertical_newton(n, f1, f2)) {
        if (f1 > f2) std::swap(f1, f2);
        std::cout << "factors => " << f1.get_str() << " * " << f2.get_str() << "\n";
    } else {
        std::cout << "No factors found (likely PRIME or algorithm did not converge).\n";
    }
    return 0;
}
