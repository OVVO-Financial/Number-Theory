#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

// Function to calculate Newton's square root approximation
mpz_class Newton_sqrt(const mpz_class &n) {
    if (n < 2) return n;

    mpz_class x = n;            // Current guess
    mpz_class x_prev = 0;       // Previous guess
    mpz_class two = 2;

    // Newton iteration: x_{k+1} = (x_k + n / x_k) / 2
    while (true) {
        x_prev = x;
        x = (x + (n / x)) / two;
        if (x == x_prev || x == x_prev - 1) {
            // Once it stabilizes or changes very little, stop.
            if (x * x > n) return x - 1;
            return x;
        }
    }
}

// Newton-Raphson method to find zz
mpz_class newton_raphson_zz(const mpz_class& r, const mpz_class& n, int max_iterations = 100, mpf_class tolerance = 1e-10) {
    mpf_class zz, zz_next;
    mpf_set_default_prec(2048); // Increase precision for accuracy

    // Initial guess: zz = sqrt(n / (2r - 3) - 3) / 2
    zz = sqrt((n.get_d() / (2 * r.get_d() - 3) - 3) / 2);

    for (int i = 0; i < max_iterations; ++i) {
        // f(zz) = (2zz + 3)(2r - 3) - n
        mpf_class f_val = (2 * zz + 3) * (2 * r.get_d() - 3) - n.get_d();

        // f'(zz) = 2(2r - 3)
        mpf_class df_val = 2 * (2 * r.get_d() - 3);

        if (abs(df_val) < tolerance) {
            break; // Avoid division by near-zero
        }

        zz_next = zz - f_val / df_val;

        // Check if the last two estimates are 1 unit apart
        if (abs(zz_next - zz) <= 1) {
            // Verify that the expression straddles n
            mpf_class expr1 = (2 * zz + 3) * (2 * r.get_d() - 3);
            mpf_class expr2 = (2 * zz_next + 3) * (2 * r.get_d() - 3);

            if ((expr1 > n && expr2 < n) || (expr1 < n && expr2 > n)) {
                break; // Converged and satisfies the straddling condition
            }
        }

        zz = zz_next;
    }

    // Round zz to the nearest integer
    mpz_class int_zz;
    mpf_class temp_zz(zz);
    temp_zz = floor(temp_zz + 0.5); // Round to nearest integer
    mpz_set_f(int_zz.get_mpz_t(), temp_zz.get_mpf_t());

    return int_zz;
}

// Function to generate quadratic residues for a given modulo
std::vector<int> generate_quadratic_residues(int modulo) {
    std::vector<int> residues = {0};
    for (int i = 1; i <= modulo / 2; ++i) {
        residues.push_back((i * i) % modulo);
    }
    // Remove duplicates
    std::sort(residues.begin(), residues.end());
    residues.erase(std::unique(residues.begin(), residues.end()), residues.end());
    return residues;
}

// Function to check if a number is a perfect square
bool is_perfect_square(const mpz_class& n) {
    mpz_class sqrt_n = Newton_sqrt(n);
    return (sqrt_n * sqrt_n == n);
}

// Fermat's Factorization with quadratic residues optimization
void fermat_factorization(const mpz_class& n, int modulo, const std::vector<int>& is_square_sieve, const std::vector<int>& Fermat_real_sieve) {
    // Generate quadratic residues
    std::vector<int> residues = generate_quadratic_residues(modulo);

    // Initialize F_realS
    mpz_class F_realS = Newton_sqrt(n) + 1;
    mpz_class max_real = (n - 9) / 6 + 3;

    while (F_realS < max_real) {
        // Compute b^2 = F_realS^2 - n
        mpz_class b2 = F_realS * F_realS - n;

        // Check if b2 % 10 is in is_square_sieve
        int b2_mod10 = mpz_class(b2 % 10).get_si();
        if (std::find(is_square_sieve.begin(), is_square_sieve.end(), b2_mod10) != is_square_sieve.end()) {
            // Check if b2 % modulo is in residues
            int b2_mod = mpz_class(b2 % modulo).get_si();
            if (std::find(residues.begin(), residues.end(), b2_mod) != residues.end()) {
                // Check if b2 is a perfect square
                if (is_perfect_square(b2)) {
                    mpz_class b = Newton_sqrt(b2);
                    std::cout << "Fermat factorization: " << F_realS - b << " and " << F_realS + b << std::endl;
                    return;
                }
            }
        }

        // Update F_realS
        F_realS += 10;
        for (int i = 0; i < modulo / 10; ++i) {
            int F_realS_mod = mpz_class(F_realS % modulo).get_si();
            if (std::find(Fermat_real_sieve.begin(), Fermat_real_sieve.end(), F_realS_mod) == Fermat_real_sieve.end()) {
                F_realS += 10;
            } else {
                break;
            }
        }
    }

    std::cout << "Fermat factorization failed." << std::endl;
}

// Trial Division method
bool trial_division(const mpz_class& n) {
    mpz_class TD = 3; // Start from 3
    while (TD * TD <= n) {
        if (n % TD == 0) {
            std::cout << "Trial Division factorization: " << TD << " and " << n / TD << std::endl;
            return true; // Factor found
        }
        TD += 2; // Only check odd numbers
    }
    std::cout << "Trial Division failed." << std::endl;
    return false; // No factor found
}

// Main SCF function with integrated complex trial multiplication, Fermat's method, and trial division
void SCF(mpz_class n, mpz_class modulo) {
    std::string class_n = (n % 4 == 1) ? "4k+1" : "4k-1";
    std::cout << "Class of n: " << class_n << std::endl;

    // Step 1: Initialization
    mpz_class r = Newton_sqrt(n);
    mpz_class max_im = (n - 9) / 6;
    mpz_class max_real = max_im + 3;

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

    // Generate square sieve from imaginary_sieve
    std::vector<int> square_sieve;
    for (int i : imaginary_sieve) {
        square_sieve.push_back((i * i) % 10);
    }
    // Remove duplicates
    std::sort(square_sieve.begin(), square_sieve.end());
    square_sieve.erase(std::unique(square_sieve.begin(), square_sieve.end()), square_sieve.end());

    std::cout << "Square Sieve: ";
    for (int i : square_sieve) std::cout << i << " ";
    std::cout << std::endl;

    // Complex Trial Multiplication Integration
    mpz_class zz = newton_raphson_zz(r, n);
    mpz_class TM_real = r + zz;
    mpz_class TM_imaginary = r - 3 - zz;

    std::cout << "Initial TM_real: " << TM_real << std::endl;
    std::cout << "Initial TM_imaginary: " << TM_imaginary << std::endl;

    // Clear vectors if they might have been used before
    std::vector<mpz_class> TMRS, TMIS, TMdRS, TMdIS;
    TMRS.clear();
    TMIS.clear();
    TMdRS.clear();
    TMdIS.clear();

    // Adjust TM_real for each sieve element
    for (size_t i = 0; i < real_sieve.size(); ++i) {
        mpz_class last_digit_real = TM_real % 10;
        int last_digit_real_int = last_digit_real.get_si();
        mpz_class adjusted_TM_real = TM_real + ((real_sieve[i] - last_digit_real_int + 10) % 10);
        TMRS.push_back(adjusted_TM_real);
        TMdRS.push_back(adjusted_TM_real); // TMdRS starts the same as TMRS
    }

    // Adjust TM_imaginary for each sieve element
    for (size_t i = 0; i < imaginary_sieve.size(); ++i) {
        mpz_class last_digit_imag = TM_imaginary % 10;
        int last_digit_imag_int = last_digit_imag.get_si();
        mpz_class adjusted_TM_imaginary = TM_imaginary + ((imaginary_sieve[i] - last_digit_imag_int + 10) % 10);

        // Ensure TMIS <= TMRS
        if (adjusted_TM_imaginary > TMRS[i]) {
            adjusted_TM_imaginary -= 10; // Decrease by 10 to satisfy TMIS <= TMRS
        }

        TMIS.push_back(adjusted_TM_imaginary);
        TMdIS.push_back(adjusted_TM_imaginary); // TMdIS starts the same as TMIS
    }

    std::cout << "TMRS after sync: ";
    for (const auto &val : TMRS) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "TMIS after sync: ";
    for (const auto &val : TMIS) std::cout << val << " ";
    std::cout << std::endl;

    // Print initial p's, q's, and n's
    std::cout << "Initial p's: ";
    for (size_t i = 0; i < TMRS.size(); ++i) {
        std::cout << TMRS[i] - TMIS[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Initial q's: ";
    for (size_t i = 0; i < TMRS.size(); ++i) {
        std::cout << TMRS[i] + TMIS[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Initial n's: ";
    for (size_t i = 0; i < TMRS.size(); ++i) {
        mpz_class p = TMRS[i] - TMIS[i];
        mpz_class q = TMRS[i] + TMIS[i];
        std::cout << p * q << " ";
    }
    std::cout << std::endl;

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    // Trial Division
    if (trial_division(n)) {
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
        return; // Stop if a factor is found
    }

    // Fermat's Factorization
    std::vector<int> is_square_sieve = square_sieve; // Use the square sieve
    std::vector<int> Fermat_real_sieve = {1, 3, 7, 9}; // Example for modulo 10
    fermat_factorization(n, modulo.get_si(), is_square_sieve, Fermat_real_sieve);

    // Stop timer
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