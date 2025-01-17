#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <chrono>

// Function to calculate Newton's square root approximation
mpz_class Newton_sqrt(const mpz_class &n) {
    mpz_class x = n, y = (n + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

// Function to compute f(zz)
mpf_class f(mpf_class zz, const mpz_class& r, const mpz_class& n) {
    mpf_class r_mpf = r.get_d();
    mpf_class product = (r_mpf + zz) * (r_mpf - 3 - zz);
    return product - n.get_d();
}

// Derivative of f(zz)
mpf_class df(mpf_class zz, const mpz_class& r) {
    mpf_class r_mpf = r.get_d();
    return 2 * (r_mpf - 3 - zz) - 2 * zz;
}

// Newton-Raphson method to find zz
mpz_class newton_raphson_zz(const mpz_class& r, const mpz_class& n, int max_iterations = 100, mpf_class tolerance = 1e-10) {
    mpf_class zz, zz_next;
    mpf_set_default_prec(512); // High precision for calculations

    zz = (r.get_d() - 3) / 2;

    for (int i = 0; i < max_iterations; ++i) {
        mpf_class f_val = f(zz, r, n);
        mpf_class df_val = df(zz, r);

        if (abs(df_val) < tolerance) {
            break;
        }

        zz_next = zz - f_val / df_val;

        if (abs(zz_next - zz) < tolerance) {
            break;
        }

        zz = zz_next;
    }

    mpz_class int_zz;
    mpf_class temp_zz(zz);
    temp_zz = floor(temp_zz + 0.5); // Round to nearest integer
    mpz_set_f(int_zz.get_mpz_t(), temp_zz.get_mpf_t());

    return int_zz;
}

// Synchronize starting point with the sieve
void sync_to_sieve(mpz_class &value, const std::vector<int> &sieve) {
    int last_digit = value.get_si() % 10;
    auto it = std::find(sieve.begin(), sieve.end(), last_digit);
    if (it != sieve.end()) {
        value -= last_digit - *it; // Align with the sieve
    } else {
        value += sieve.front() - last_digit; // Use the first sieve element if misaligned
    }
}

// Main SCF function with integrated complex trial multiplication, Fermat's method, and trial division
void SCF(mpz_class n, mpz_class modulo) {
    std::string class_n = (n % 4 == 1) ? "4k+1" : "4k-1";
    std::cout << "Class of n: " << class_n << std::endl;

    // Step 1: Initialization
    mpz_class r = Newton_sqrt(n);
    mpz_class max_im = (n - 9) / 6;
    mpz_class max_real = max_im + 3;
    mpz_class TD = 3; // Trial Division starts from 3
    mpz_class Fermat_real = r + 1;

    std::cout << "Root of n: " << r << std::endl;

    // Step 2: Fermat Sieves based on class of n
    int last_digit = mpz_class(n % 10).get_si();
    std::vector<int> real_sieve, imaginary_sieve;
    std::string sieve_selected;

    if (class_n == "4k-1") { // N = 4k - 1
        sieve_selected = "4k-1 Sieve";
        switch (last_digit) {
            case 1: real_sieve = {0, 0, 4, 6}; imaginary_sieve = {3, 7, 5, 5}; break;
            case 3: real_sieve = {2, 2, 8, 8}; imaginary_sieve = {1, 9, 9, 1}; break;
            case 7: real_sieve = {4, 4, 6, 6}; imaginary_sieve = {3, 7, 7, 3}; break;
            case 9: real_sieve = {0, 0, 2, 8}; imaginary_sieve = {1, 9, 5, 5}; break;
            default: std::cout << "Unexpected last digit for 4k-1 case." << std::endl; return;
        }
    } else { // N = 4k + 1
        sieve_selected = "4k+1 Sieve";
        switch (last_digit) {
            case 1: real_sieve = {1, 5, 5, 9}; imaginary_sieve = {0, 2, 8, 0}; break;
            case 3: real_sieve = {3, 3, 7, 7}; imaginary_sieve = {4, 6, 6, 4}; break;
            case 7: real_sieve = {1, 1, 9, 9}; imaginary_sieve = {2, 8, 8, 2}; break;
            case 9: real_sieve = {3, 5, 5, 7}; imaginary_sieve = {0, 4, 6, 0}; break;
            default: std::cout << "Unexpected last digit for 4k+1 case." << std::endl; return;
        }
    }

    std::cout << "Selected Real Sieve: ";
    for (int i : real_sieve) std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "Selected Imaginary Sieve: ";
    for (int i : imaginary_sieve) std::cout << i << " ";
    std::cout << std::endl;

    // Complex Trial Multiplication Integration
    mpz_class zz = newton_raphson_zz(r, n);
    mpz_class TM_real = r + zz;
    mpz_class TM_imaginary = r - 3 - zz;

    std::cout << "Initial TM_real: " << TM_real << std::endl;
    std::cout << "Initial TM_imaginary: " << TM_imaginary << std::endl;

    sync_to_sieve(TM_real, real_sieve);
    sync_to_sieve(TM_imaginary, imaginary_sieve);

    std::cout << "Updated TM_real: " << TM_real << std::endl;
    std::cout << "Updated TM_imaginary: " << TM_imaginary << std::endl;

    std::vector<mpz_class> TMRS, TMIS, TMdRS, TMdIS;
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        TMRS.push_back(TM_real + ((real_sieve[i] - TM_real.get_si() % 10 + 10) % 10));
        TMdRS.push_back(TMRS.back());
        TMIS.push_back(TM_imaginary + ((imaginary_sieve[i] - TM_imaginary.get_si() % 10 + 10) % 10));
        TMdIS.push_back(TMIS.back());
    }

    auto start = std::chrono::high_resolution_clock::now();
    int iteration = 0;

    while(true) {
        ++iteration;
        bool factor_found = false;

        for (size_t i = 0; i < TMRS.size(); ++i) {
            mpz_class p = TMRS[i] - TMIS[i];
            mpz_class q = TMRS[i] + TMIS[i];
            mpz_class p_d = TMdRS[i] - TMdIS[i];
            mpz_class q_d = TMdRS[i] + TMdIS[i];

            mpz_class N = p * q;
            mpz_class N_d = p_d * q_d;

            if (N == n) {
                factor_found = true;
                std::cout << "Factor found in ascending: " << p << " and " << q << std::endl;
                goto end;
            }
            if (N_d == n) {
                factor_found = true;
                std::cout << "Factor found in descending: " << p_d << " and " << q_d << std::endl;
                goto end;
            }

            if (N > n) TMIS[i] += 10; else TMRS[i] += 10;
            if (N_d < n) TMdIS[i] -= 10; else TMdRS[i] -= 10;

            // Fermat Method check after each trial multiplication
            mpz_class b2 = Fermat_real * Fermat_real - n;
            if (b2 >= 0) {
                mpz_class b = Newton_sqrt(b2);
                if (b * b == b2) {
                    std::cout << "Fermat factorization: " << Fermat_real - b << " and " << Fermat_real + b << std::endl;
                    factor_found = true;
                    goto end;
                }
            }
            Fermat_real += 10; // Next number in Fermat sequence

            // Trial Division - check if n is divisible by TD
            if (n % TD == 0) {
                std::cout << "Trial Division factorization: " << TD << " and " << n / TD << std::endl;
                factor_found = true;
                goto end;
            }
            TD += 2; // Only check odd numbers, skip even numbers except 2 which we already handled implicitly

            if (factor_found) break; // Break the inner loop if a factor is found
        }

        // Check if all descending paths are exhausted
        if (std::all_of(TMdRS.begin(), TMdRS.end(), [TMdIS](const mpz_class &real) {
            return std::none_of(TMdIS.begin(), TMdIS.end(), [&real](const mpz_class &imag) {
                return real - imag >= 3;
            });
        })) {
            std::cout << "Number might be prime or factorization not found with current methods." << std::endl;
            goto end;
        }
    }

end:
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
}

int main() {
    mpz_class n, modulo;
    std::cout << "Enter number to factorize: ";
    std::cin >> n;
    std::cout << "Enter modulo: ";
    std::cin >> modulo;
    SCF(n, modulo);
    return 0;
}