#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <chrono>

// Function to calculate Newton's square root approximation
mpz_class Newton_sqrt(const mpz_class& n) {
    mpz_class x = n, y = (n + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

// Helper function to find the smallest positive shift or the largest negative shift if all are negative
int find_shift(const std::vector<int>& diffs) {
    auto non_negative = std::find_if(diffs.begin(), diffs.end(), [](int diff) { return diff >= 0; });
    if (non_negative != diffs.end()) {
        return *std::min_element(non_negative, diffs.end());
    } else {
        return diffs.back();
    }
}

// Main SCF function with enforced sieve
void SCF(mpz_class n, mpz_class modulo) {
    // Determine class of n
    std::string class_n = (n % 4 == 1) ? "4k+1" : "4k-1";
    std::cout << "Class of n: " << class_n << std::endl;

    // Step 1: Initialization
    mpz_class r = Newton_sqrt(n);
    mpz_class max_im = (n - 9) / 6;
    mpz_class max_real = max_im + 3;
    mpz_class TD = 3;
    mpz_class Fermat_real = r + 1;

    // Step 2: Fermat Sieves based on class of n
    int last_digit = mpz_class(n % 10).get_si();
    std::vector<int> Fermat_real_sieve, real_sieve, imaginary_sieve;
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
        Fermat_real_sieve = real_sieve;
    } else { // N = 4k + 1
        sieve_selected = "4k+1 Sieve";
        switch (last_digit) {
            case 1: real_sieve = {1, 5, 5, 9}; imaginary_sieve = {0, 2, 8, 0}; break;
            case 3: real_sieve = {3, 3, 7, 7}; imaginary_sieve = {4, 6, 6, 4}; break;
            case 7: real_sieve = {1, 1, 9, 9}; imaginary_sieve = {2, 8, 8, 2}; break;
            case 9: real_sieve = {3, 5, 5, 7}; imaginary_sieve = {0, 4, 6, 0}; break;
            default: std::cout << "Unexpected last digit for 4k+1 case." << std::endl; return;
        }
        Fermat_real_sieve = real_sieve;
    }

    // Step 3: Square Check Fermat Sieve
    std::vector<int> square_sieve;
    for (int i : imaginary_sieve) {
        int sq = (i * i) % 10;
        if (std::find(square_sieve.begin(), square_sieve.end(), sq) == square_sieve.end()) {
            square_sieve.push_back(sq);
        }
    }

    // Correct Calculation for Resulting Reals 
    std::vector<int> resulting_reals1, resulting_reals2;
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        int p = (real_sieve[i] - imaginary_sieve[i] + 10) % 10; // Addition
        int q = (real_sieve[i] + imaginary_sieve[i]) % 10;      // Subtraction
        resulting_reals1.push_back(p);
        resulting_reals2.push_back(q);
    }

    // SYNC COMPLEX TRIAL MULTIPLICATION IMAGINARY TO SIEVE STARTING POINT
    mpz_class zz = Newton_sqrt(r - 3); // New initial approximation based on your instruction
    mpz_class TM_real = r + zz;
    mpz_class TM_imaginary = r - 3 - zz;

    // SYNC COMPLEX TRIAL MULTIPLICATION IMAGINARY TO SIEVE STARTING POINT
    int last_digit_im = mpz_class(TM_imaginary % 10).get_si();
    std::vector<int> im_init_diff;
    for (int i : imaginary_sieve) {
        im_init_diff.push_back(i - last_digit_im);
    }

    // Shift to align with the sieve while moving down
    if (std::find_if(im_init_diff.begin(), im_init_diff.end(), [](int diff) { return diff >= 0; }) != im_init_diff.end()) {
        int shift = find_shift(im_init_diff);
        if (shift > 0) {
            TM_imaginary -= shift; // Move down
        } else {
            TM_imaginary += shift;
        }
    } else {
        TM_imaginary += im_init_diff.back();
    }

    // FIND & SYNC COMPLEX TRIAL MULTIPLICATION REAL
    int last_digit_real = mpz_class(TM_real % 10).get_si();
    std::vector<int> real_init_diff;
    for (int i : real_sieve) {
        real_init_diff.push_back(i - last_digit_real);
    }

    // Shift to align with the sieve while moving left
    if (std::find_if(real_init_diff.begin(), real_init_diff.end(), [](int diff) { return diff >= 0; }) != real_init_diff.end()) {
        int shift = find_shift(real_init_diff);
        if (shift > 0) {
            TM_real -= shift; // Move left
        } else {
            TM_real += shift;
        }
    } else {
        TM_real += real_init_diff.back();
    }

    // Step 5: Prepare for factorization
    std::vector<mpz_class> TMRS, TMIS, TMdRS, TMdIS;
    for (int i : real_sieve) {
        TMRS.push_back(TM_real + ((i - TM_real.get_si() % 10 + 10) % 10));  // Ascending real part
        TMdRS.push_back(TM_real - ((i - TM_real.get_si() % 10 + 10) % 10));  // Descending real part
    }
    for (int i : imaginary_sieve) {
        TMIS.push_back(TM_imaginary + ((i - TM_imaginary.get_si() % 10 + 10) % 10));  // Ascending imaginary part
        TMdIS.push_back(TM_imaginary - ((i - TM_imaginary.get_si() % 10 + 10) % 10));  // Descending imaginary part
    }

    // Adjust for initial sync mismatch
    for (size_t i = 0; i < TMRS.size(); ++i) {
        if (TMRS[i] - TMIS[i] <= 1) TMRS[i] += 10; // Ensure ascending stays above
        if (TMdRS[i] - TMdIS[i] <= 1) TMdRS[i] -= 10; // Ensure descending stays below
    }

    // Fermat Difference of Squares with Square Sieve applied
    std::vector<mpz_class> F_realS;
    for (int i : Fermat_real_sieve) {
        F_realS.push_back(Fermat_real + ((i - Fermat_real.get_si() % 10 + 10) % 10));
    }

    auto start = std::chrono::high_resolution_clock::now();

    int iterations = 0;
    while (true) {
        for (size_t i = 0; i < TMRS.size(); ++i) {
            mpz_class p = TMRS[i] - TMIS[i], p_d = TMdRS[i] - TMdIS[i];
            mpz_class q = TMRS[i] + TMIS[i], q_d = TMdRS[i] + TMdIS[i];
            mpz_class N = p * q, N_d = p_d * q_d;

            if (N == n) {
                std::cout << "TM ascending: " << p << " and " << q << std::endl;
                goto end;
            }
            if (N_d == n) {
                std::cout << "TM descending: " << p_d << " and " << q_d << std::endl;
                goto end;
            }

            // Here, we adjust both ascending and descending:
            if (N < n) {
                TMIS[i] += 10;  // Increase imaginary if too high for ascending
            } else {
                TMRS[i] += 10;  // Increase real if too low for ascending
            }
            if (N_d > n) {
                TMdIS[i] -= 10;  // Decrease imaginary if too low for descending
            } else {
                TMdRS[i] -= 10;  // Decrease real if too high for descending
            }
        }

        for (size_t i = 0; i < F_realS.size(); ++i) {
            mpz_class b2 = F_realS[i] * F_realS[i] - n;
            int last_digit_b2 = mpz_class(b2 % 10).get_si();
            if (b2 >= 0 && std::find(square_sieve.begin(), square_sieve.end(), last_digit_b2) != square_sieve.end()) {
                mpz_class b_test = Newton_sqrt(b2);
                if (b_test * b_test == b2) {
                    std::cout << "Fermat: " << F_realS[i] - b_test << " and " << F_realS[i] + b_test << std::endl;
                    goto end;
                }
            }
            F_realS[i] += 10; // Increment by 10 to match the sieve's pattern
        }

        // Trial Division 
        while (TD % 10 == 5) TD += 2; // Skip numbers ending in 5 as they cannot be prime except for 5 itself
        if (n % TD == 0) {
            std::cout << "Trial Division: " << TD << " and " << n / TD << std::endl;
            goto end;
        }
        TD += 2;

        ++iterations;
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