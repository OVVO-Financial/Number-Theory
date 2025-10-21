#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <limits>

// Function to calculate Newton's square root approximation
mpz_class Newton_sqrt(const mpz_class &n) {
    mpz_class x = n, y = (n + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

// Newton-Raphson method to find zz based on new equation
mpz_class newton_raphson_zz(const mpz_class& r, const mpz_class& n, int max_iterations = 100, mpf_class tolerance = 1e-10) {
    mpf_class zz, zz_next;
    mpf_set_default_prec(1024); // Higher precision due to quadratic term

    // Initial guess - adjust based on new equation
    zz = sqrt((n.get_d() / (2 * r.get_d() - 3) - 3) / 2);  

    for (int i = 0; i < max_iterations; ++i) {
        mpf_class f_val = (2 * r.get_d() - 3) * (2 * zz * zz + 3) - n.get_d();
        mpf_class df_val = (2 * r.get_d() - 3) * 4 * zz;

        if (abs(df_val) < tolerance) {
            break; // Avoid division by near-zero
        }

        zz_next = zz - f_val / df_val;

        if (abs(zz_next - zz) < tolerance) {
            break; // Converged
        }

        zz = zz_next;
    }

    mpz_class int_zz;
    mpf_class temp_zz(zz);
    temp_zz = floor(temp_zz + 0.5); // Round to nearest integer
    mpz_set_f(int_zz.get_mpz_t(), temp_zz.get_mpf_t());

    return int_zz;
}

// Function to perform complex trial multiplication
void complex_trial_multiplication(mpz_class n) {
    // Determine the class of n
    std::string class_n = (n % 4 == 1) ? "4k+1" : "4k-1";
    std::cout << "Class of n: " << class_n << std::endl;

    // Step 1: Initialization
    mpz_class r = Newton_sqrt(n);
    std::cout << "Square root of n: " << r << std::endl;

    // Use Newton-Raphson method to find zz
    mpz_class zz = newton_raphson_zz(r, n);
    mpz_class TM_real = r + zz;
    mpz_class TM_imaginary = r - 3 - zz;

    std::cout << "Initial TM_real: " << TM_real << std::endl;
    std::cout << "Initial TM_imaginary: " << TM_imaginary << std::endl;

    std::vector<int> real_sieve, imaginary_sieve;
    std::string sieve_selected;
    int last_digit = n.get_si() % 10;

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

    // Initialize vectors
    std::vector<mpz_class> TMRS, TMIS, TMdRS, TMdIS;

    // Adjust TM_real and TM_imaginary for each sieve element
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        mpz_class last_digit_real = TM_real % 10;
        int last_digit_real_int = last_digit_real.get_si();
        mpz_class adjusted_TM_real = TM_real + ((real_sieve[i] - last_digit_real_int + 10) % 10);
        TMRS.push_back(adjusted_TM_real);
        TMdRS.push_back(adjusted_TM_real); // TMdRS starts the same as TMRS
    }

    for (size_t i = 0; i < imaginary_sieve.size(); ++i) {
        mpz_class last_digit_imag = TM_imaginary % 10;
        int last_digit_imag_int = last_digit_imag.get_si();
        mpz_class adjusted_TM_imaginary = TM_imaginary + ((imaginary_sieve[i] - last_digit_imag_int + 10) % 10);
        TMIS.push_back(adjusted_TM_imaginary);
        TMdIS.push_back(adjusted_TM_imaginary); // TMdIS starts the same as TMIS
    }

    std::cout << "TMRS after sync: ";
    for (const auto &val : TMRS) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "TMIS after sync: ";
    for (const auto &val : TMIS) std::cout << val << " ";
    std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    int iteration = 0;
    bool factor_found = false;

    while(true) {
        ++iteration;

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
                break;
            }
            if (N_d == n) {
                factor_found = true;
                std::cout << "Factor found in descending: " << p_d << " and " << q_d << std::endl;
                break;
            }

            // Adjust ascending and descending paths
            if (N > n) TMIS[i] += 10; else TMRS[i] += 10;
            if (N_d < n) TMdIS[i] -= 10; else TMdRS[i] -= 10;
        }

        if (factor_found) break;

        // Stop if all descending paths are exhausted
        if (std::all_of(TMdRS.begin(), TMdRS.end(), [&TMdIS](const mpz_class &real) {
            return std::none_of(TMdIS.begin(), TMdIS.end(), [&real](const mpz_class &imag) {
                return real - imag >= 3;
            });
        })) {
            std::cout << "Prime" << std::endl;
            break;
        }
    }

    // Print vectors only on the final iteration
    std::cout << "Final State:" << std::endl;
    std::cout << "  TMRS: ";
    for (const auto &val : TMRS) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "  TMIS: ";
    for (const auto &val : TMIS) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "  TMdRS: ";
    for (const auto &val : TMdRS) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "  TMdIS: ";
    for (const auto &val : TMdIS) std::cout << val << " ";
    std::cout << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
}

int main() {
    mpz_class n;
    std::cout << "Enter number to factorize: ";
    std::cin >> n;
    complex_trial_multiplication(n);
    return 0;
}