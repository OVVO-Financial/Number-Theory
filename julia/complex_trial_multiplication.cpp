#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>

// Function to calculate Newton's square root approximation
mpz_class Newton_sqrt(const mpz_class &n) {
    mpz_class x = n, y = (n + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
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
    std::cout << "Synchronized value to: " << value << std::endl;
}

// Function to perform complex trial multiplication
void complex_trial_multiplication(mpz_class n) {
    // Determine the class of n
    std::string class_n = (n % 4 == 1) ? "4k+1" : "4k-1";
    std::cout << "Class of n: " << class_n << std::endl;

    // Step 1: Initialization
    mpz_class r = Newton_sqrt(n);
    std::cout << "Square root of n: " << r << std::endl;
    
    // Using mpf_class for a more precise initial approximation
    mpf_class mpf_r = r.get_d();
    mpf_class mpf_zz = (sqrt(mpf_r - 3) + (mpf_r - 3)) / 2;
    mpz_class zz = mpf_zz.get_si();
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

    // Print selected sieves
    std::cout << "Sieve Selected: " << sieve_selected << std::endl;
    std::cout << "Real Sieve: ";
    for (int i : real_sieve) std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "Imaginary Sieve: ";
    for (int i : imaginary_sieve) std::cout << i << " ";
    std::cout << std::endl;

    // Sync TM_real and TM_imaginary to the sieve
    sync_to_sieve(TM_real, real_sieve);
    sync_to_sieve(TM_imaginary, imaginary_sieve);

    std::cout << "Updated TM_real: " << TM_real << std::endl;
    std::cout << "Updated TM_imaginary: " << TM_imaginary << std::endl;

    // Step 3: Square Check Fermat Sieve (not used in this function, but included for completeness)
    std::vector<int> square_sieve;
    for (int i : imaginary_sieve) {
        int sq = (i * i) % 10;
        if (std::find(square_sieve.begin(), square_sieve.end(), sq) == square_sieve.end()) {
            square_sieve.push_back(sq);
        }
    }

    // Correct Calculation for Resulting Reals (not used, but included for completeness)
    std::vector<int> resulting_reals1, resulting_reals2;
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        int p = (real_sieve[i] - imaginary_sieve[i] + 10) % 10; // Addition
        int q = (real_sieve[i] + imaginary_sieve[i]) % 10;      // Subtraction
        resulting_reals1.push_back(p);
        resulting_reals2.push_back(q);
    }

    std::vector<mpz_class> TMRS, TMIS, TMdRS, TMdIS;
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        TMRS.push_back(TM_real + ((real_sieve[i] - TM_real.get_si() % 10 + 10) % 10));  // Ascending real part
        TMdRS.push_back(TMRS.back());  // TMdRS set equal to TMRS initially
        TMIS.push_back(TM_imaginary + ((imaginary_sieve[i] - TM_imaginary.get_si() % 10 + 10) % 10));  // Ascending imaginary part
        // TMdIS is set to be equal to TMIS after sync
        TMdIS.push_back(TMIS.back());  
    }

    auto start = std::chrono::high_resolution_clock::now();
    int iteration = 0;

    // Iterative factorization
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
                break;
            }
            if (N_d == n) {
                factor_found = true;
                std::cout << "Factor found in descending: " << p_d << " and " << q_d << std::endl;
                break;
            }

            // Adjust ascending and descending paths - TMdIS and TMIS diverge after 1st iteration
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