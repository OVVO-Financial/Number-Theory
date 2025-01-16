#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>

// Function to calculate Newton's square root approximation
mpz_class Newton_sqrt(const mpz_class& n) {
    mpz_class x = n, y = (n + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

// Helper function to calculate the square root of a complex number approximation
void complex_sqrt(mpz_class& real, mpz_class& imag, const mpz_class& n, const mpz_class& max_im, const mpz_class& max_real) {
    double r = sqrt(n.get_d() + max_im.get_d() * max_real.get_d());
    double theta = atan2(max_im.get_d(), n.get_d()) / 2.0;
    real = mpz_class(round(r * cos(theta)));
    imag = mpz_class(round(r * sin(theta)));
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
            case 1: real_sieve = {4, 6, 0, 0}; imaginary_sieve = {5, 5, 3, 7}; break;
            case 3: real_sieve = {2, 8, 2, 8}; imaginary_sieve = {1, 9, 9, 1}; break;
            case 7: real_sieve = {4, 6, 4, 6}; imaginary_sieve = {3, 7, 7, 3}; break;
            case 9: real_sieve = {0, 2, 8, 0}; imaginary_sieve = {1, 5, 5, 9}; break;
            default: std::cout << "Unexpected last digit for 4k-1 case." << std::endl; return;
        }
        Fermat_real_sieve = real_sieve;
    } else { // N = 4k + 1
        sieve_selected = "4k+1 Sieve";
        switch (last_digit) {
            case 1: real_sieve = {1, 5, 5, 9}; imaginary_sieve = {0, 2, 8, 0}; break;
            case 3: real_sieve = {3, 7, 3, 7}; imaginary_sieve = {4, 6, 6, 4}; break;
            case 7: real_sieve = {1, 9, 1, 9}; imaginary_sieve = {2, 8, 8, 2}; break;
            case 9: real_sieve = {3, 5, 5, 7}; imaginary_sieve = {0, 4, 6, 0}; break;
            default: std::cout << "Unexpected last digit for 4k+1 case." << std::endl; return;
        }
        Fermat_real_sieve = real_sieve;
    }

    std::cout << "Sieve Selected: " << sieve_selected << std::endl;
    std::cout << "Real Sieve: "; for (int i : real_sieve) std::cout << i << " "; std::cout << std::endl;
    std::cout << "Imaginary Sieve: "; for (int i : imaginary_sieve) std::cout << i << " "; std::cout << std::endl;

    // Step 3: Square Check Fermat Sieve
    std::vector<int> square_sieve;
    for (int i : imaginary_sieve) {
        int sq = (i * i) % 10;
        if (std::find(square_sieve.begin(), square_sieve.end(), sq) == square_sieve.end()) {
            square_sieve.push_back(sq);
        }
    }
    std::cout << "Square Sieve: "; for (int i : square_sieve) std::cout << i << " "; std::cout << std::endl;

    // Correct Calculation for Resulting Reals 
    std::vector<int> resulting_reals1, resulting_reals2;
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        int p = (real_sieve[i] - imaginary_sieve[i] + 10) % 10; // Addition
        int q = (real_sieve[i] + imaginary_sieve[i]) % 10;      // Subtraction
        resulting_reals1.push_back(p);
        resulting_reals2.push_back(q);
    }

    std::cout << "Resulting Reals (p,g) Vector 1: "; for (int i : resulting_reals1) std::cout << i << " "; std::cout << std::endl;
    std::cout << "Resulting Reals (p,g) Vector 2: "; for (int i : resulting_reals2) std::cout << i << " "; std::cout << std::endl;

    // Step 4: Sync Complex Trial Multiplication
    mpz_class TM_imaginary, TM_real;
    complex_sqrt(TM_real, TM_imaginary, n, max_im, max_real);
    int last_digit_im = TM_imaginary.get_si() % 10;
    std::vector<int> im_init_diff;

    for (int i : imaginary_sieve) {
        im_init_diff.push_back(i - last_digit_im);
    }

    int im_shift = 0;
    for (int diff : im_init_diff) {
        if (diff >= 0) {
            im_shift = diff;
            break;
        }
        im_shift = diff;
    }
    TM_imaginary += im_shift;

    if (TM_real < TM_imaginary) TM_real = TM_imaginary + 3;
    int last_digit_real = TM_real.get_si() % 10;
    std::vector<int> real_init_diff;

    for (int i : real_sieve) {
        real_init_diff.push_back(i - last_digit_real);
    }

    int real_shift = 0;
    for (int diff : real_init_diff) {
        if (diff >= 0) {
            real_shift = diff;
            break;
        }
        real_shift = diff;
    }
    TM_real += real_shift;

    // Step 5: Prepare for factorization
    std::vector<mpz_class> TMRS, TMIS, TMdRS, TMdIS;
    for (int i : real_sieve) {
        TMRS.push_back(TM_real + ((i - TM_real.get_si() % 10 + 10) % 10));
        TMdRS.push_back(TM_real + ((i - TM_real.get_si() % 10 + 10) % 10));
    }
    for (int i : imaginary_sieve) {
        TMIS.push_back(TM_imaginary + ((i - TM_imaginary.get_si() % 10 + 10) % 10));
        TMdIS.push_back(TM_imaginary + ((i - TM_imaginary.get_si() % 10 + 10) % 10));
    }

    // Adjust for initial sync mismatch
    for (size_t i = 0; i < TMRS.size(); ++i) {
        if (TMRS[i] - TMIS[i] <= 1) TMRS[i] += 10;
        if (TMdRS[i] - TMdIS[i] <= 1) TMdIS[i] -= 10;
    }

    // Step 6: Main loop for factorization
    int iterations = 0;
    while (iterations < 10000) { // Limit iterations to prevent infinite loop
        for (size_t i = 0; i < TMRS.size(); ++i) {
            mpz_class p = TMRS[i] - TMIS[i], p_d = TMdRS[i] - TMdIS[i];
            mpz_class q = TMRS[i] + TMIS[i], q_d = TMdRS[i] + TMdIS[i];
            mpz_class N = p * q, N_d = p_d * q_d;

            if (N == n) {
                std::cout << "TM ascending: " << p << " and " << q << std::endl;
                return;
            }
            if (N_d == n) {
                std::cout << "TM descending: " << p_d << " and " << q_d << std::endl;
                return;
            }

            if (N > n) {
                TMIS[i] += 10;
            } else {
                TMRS[i] += 10;
            }
            if (N_d < n) {
                TMdIS[i] -= 10;
            } else {
                TMdRS[i] -= 10;
            }
        }

        // Fermat Difference of Squares 
        std::vector<mpz_class> F_realS;
        for (int i : Fermat_real_sieve) {
            F_realS.push_back(Fermat_real + ((i - Fermat_real.get_si() % 10 + 10) % 10));
        }

        for (mpz_class fr : F_realS) {
            mpz_class b2 = fr * fr - n;
            if (b2 >= 0) {
                mpz_class b_test = Newton_sqrt(b2);
                if (b_test * b_test == b2) {
                    std::cout << "Fermat: " << fr - b_test << " and " << fr + b_test << std::endl;
                    return;
                }
            }
        }

        // Trial Division 
        while (TD % 10 == 5) TD += 2; // Skip numbers ending in 5 as they cannot be prime except for 5 itself
        if (n % TD == 0) {
            std::cout << "Trial Division: " << TD << " and " << n / TD << std::endl;
            return;
        }
        TD += 2;

        ++iterations;
    }

    std::cout << "No factors found within the given constraints." << std::endl;
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