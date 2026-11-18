#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <gmpxx.h>

using namespace std;

// --- 1. R^2 Tensor Space Mappings ---

struct Coordinate {
    mpz_class x;
    mpz_class y;
};

Coordinate map_phi(mpz_class a, mpz_class b) { return {a + b, a - b}; }
pair<mpz_class, mpz_class> unmap_phi(Coordinate c) { return {(c.x + c.y) / 2, (c.x - c.y) / 2}; }

mpz_class geometric_norm(Coordinate c) {
    return (c.x * c.x + c.y * c.y) / 2;
}

Coordinate tensor_multiply(Coordinate c1, Coordinate c2, mpz_class c_val) {
    mpz_class a1 = (c1.x + c1.y) / 2;
    mpz_class b1 = (c1.x - c1.y) / 2;
    mpz_class a2 = (c2.x + c2.y) / 2;
    mpz_class b2 = (c2.x - c2.y) / 2;
    
    mpz_class term_A = a1 * a2 - c_val * b1 * b2;
    mpz_class term_B = a1 * b2 + a2 * b1;
    return map_phi(term_A, term_B);
}

bool is_algebraic_square(Coordinate gamma_sq, mpz_class c_val, Coordinate& root) {
    pair<mpz_class, mpz_class> unmapped = unmap_phi(gamma_sq);
    mpz_class A = unmapped.first;
    mpz_class B = unmapped.second;
    
    mpz_class norm = A * A + c_val * B * B;
    if (sgn(norm) < 0 || !mpz_perfect_square_p(norm.get_mpz_t())) return false;
    
    mpz_class root_norm = sqrt(norm);
    mpz_class num = A + root_norm;
    if (num % 2 != 0) num = A - root_norm;
    if (num % 2 != 0) return false;
    
    mpz_class a_sq = num / 2;
    if (sgn(a_sq) < 0 || !mpz_perfect_square_p(a_sq.get_mpz_t())) return false;
    
    mpz_class a = sqrt(a_sq);
    if (a == 0) return false;
    if (B % (2 * a) != 0) return false;
    
    mpz_class b = B / (2 * a);
    root = map_phi(a, b);
    return true;
}

// --- 2. 2D Vector Geometry ---

struct Vec2 { long long a, b; };

long long dot_product(Vec2 v1, Vec2 v2) { return v1.a * v2.a + v1.b * v2.b; }
long long vec_norm(Vec2 v) { return dot_product(v, v); }

void gauss_reduce(Vec2& v1, Vec2& v2) {
    if (vec_norm(v1) > vec_norm(v2)) swap(v1, v2);
    while (true) {
        double proj = (double)dot_product(v1, v2) / (double)vec_norm(v1);
        long long q = round(proj);
        if (q == 0) break; 
        
        Vec2 v_new = {v2.a - q * v1.a, v2.b - q * v1.b};
        if (vec_norm(v_new) >= vec_norm(v1)) {
            v2 = v_new; 
            break;      
        }
        v2 = v1;
        v1 = v_new;
    }
}

// --- 3. Factor Base & Matrix Helpers ---

struct PrimeIdeal {
    int p;
    float log_p;
    vector<int> roots;
};

vector<PrimeIdeal> generate_factor_base(int max_prime, mpz_class c_val) {
    vector<PrimeIdeal> fb;
    long c = c_val.get_si();
    for (int p = 2; p <= max_prime; ++p) {
        bool is_prime = true;
        for (int i = 2; i * i <= p; ++i) {
            if (p % i == 0) { is_prime = false; break; }
        }
        if (!is_prime) continue;
        
        vector<int> roots;
        for (int r = 0; r < p; ++r) {
            if ((r * r + c) % p == 0) roots.push_back(r);
        }
        if (!roots.empty() || p == 2) fb.push_back({p, log2((float)p), roots});
    }
    return fb;
}

bool get_smooth_vector(mpz_class val, const vector<PrimeIdeal>& fb, vector<int>& vec) {
    if (val == 0) return false;
    for (const auto& ideal : fb) {
        int count = 0;
        while (val % ideal.p == 0) {
            count++;
            val /= ideal.p;
        }
        vec.push_back(count % 2);
    }
    return abs(val) == 1;
}

// --- 4. End-to-End GNFS Pipeline ---

int main(int argc, char* argv[]) {
    if (argc < 2) { cout << "Usage: ./gnfs_gmp <N>\n"; return 1; }
    
    mpz_class N(argv[1]);
    mpz_class m = sqrt(N);
    mpz_class c_val = N - (m * m); 
    
    cout << "=== Full GNFS 2D Lattice Sieve Pipeline ===\n";
    cout << "Target N: " << N << "\nBase m: " << m << ", Poly: x^2 + " << c_val << "\n\n";

// Dynamic GNFS Parameter Scaling
    double n_d = N.get_d();
    
    // 1. Dynamic Factor Base Bound
    // Scales the prime limit relative to N. Forces a minimum of 100 for tiny numbers.
    // For N = 8051, this dynamically evaluates to max_prime = 493.
    int max_prime = std::max(100.0, std::sqrt(n_d) * 5.5);
    
    // 2. Dynamic Geometric Search Grid
    // Expands the 2D search window proportionally to ensure the sieve harvests
    // enough perfectly smooth nodes to exceed the GF(2) matrix columns.
    // For N = 8051, this dynamically evaluates to M = 1479.
    int M = max_prime * 3;
    
    vector<PrimeIdeal> fb = generate_factor_base(max_prime, c_val);
    int columns = fb.size() * 2 + 1; // Rational primes + Algebraic primes + Sign bit
    
    vector<vector<float>> sieve_grid(2 * M + 1, vector<float>(M + 1, 0.0));
    cout << "[Phase 1] Executing 2D vector jumps (Gauss reduced)...\n";
    
    for (const auto& ideal : fb) {
        for (int r : ideal.roots) {
            Vec2 v1 = {ideal.p, 0}, v2 = {r, 1};
            gauss_reduce(v1, v2);
            int jump_bound = (M / max(1LL, (long long)sqrt(ideal.p))) + 2; 
            for (int i = -jump_bound; i <= jump_bound; ++i) {
                for (int j = -jump_bound; j <= jump_bound; ++j) {
                    long long a = i * v1.a + j * v2.a;
                    long long b = i * v1.b + j * v2.b;
                    if (a >= -M && a <= M && b >= 1 && b <= M) {
                        sieve_grid[a + M][b] += ideal.log_p;
                    }
                }
            }
        }
    }
    
    cout << "[Phase 2] Harvesting exact smooth relations...\n";
    vector<pair<mpz_class, mpz_class>> smooth_relations;
    vector<vector<int>> matrix;
    
    for (int a = -M; a <= M; ++a) {
        for (int b = 1; b <= M; ++b) {
            if (std::gcd(abs(a), b) != 1) continue;
            
            Coordinate mapped = map_phi(mpz_class(a), mpz_class(b));
            mpz_class exact_norm = geometric_norm(mapped);
            double expected_log = log2(exact_norm.get_d());
            
            // If the sieve marked it as a highly probable hit, exact-factor it
            if (sieve_grid[a + M][b] >= expected_log - 2.5) {
                mpz_class rational_val = a - b * m;
                vector<int> row_vector;
                
                row_vector.push_back(sgn(rational_val) < 0 ? 1 : 0);
                if (!get_smooth_vector(abs(rational_val), fb, row_vector)) continue;
                if (!get_smooth_vector(exact_norm, fb, row_vector)) continue;
                
                smooth_relations.push_back({a, b});
                matrix.push_back(row_vector);
                
                if (matrix.size() > columns + 5) break; 
            }
        }
        if (matrix.size() > columns + 5) break;
    }
    
    if (matrix.size() <= columns) {
        cout << "Error: Not enough exact smooth relations found to form a matrix.\n";
        return 1;
    }
    
    cout << "[Phase 3] Matrix built (" << matrix.size() << "x" << columns << "). Running Gaussian Elimination...\n";
    
    int rows = matrix.size();
    vector<vector<int>> history(rows, vector<int>(rows, 0));
    for (int i = 0; i < rows; ++i) history[i][i] = 1;

    int pivot_row = 0;
    for (int col = 0; col < columns && pivot_row < rows; ++col) {
        int swap_row = pivot_row;
        while (swap_row < rows && matrix[swap_row][col] == 0) swap_row++;
        if (swap_row == rows) continue; 
        
        swap(matrix[pivot_row], matrix[swap_row]);
        swap(history[pivot_row], history[swap_row]);
        
        for (int r = pivot_row + 1; r < rows; ++r) {
            if (matrix[r][col] == 1) {
                for (int c = col; c < columns; ++c) matrix[r][c] ^= matrix[pivot_row][c];
                for (int c = 0; c < rows; ++c) history[r][c] ^= history[pivot_row][c];
            }
        }
        pivot_row++;
    }

    cout << "[Phase 4] Extracting Null Space and applying Tensor Product...\n";
    
    vector<pair<mpz_class, mpz_class>> final_relations;
    for (int r = pivot_row; r < rows; ++r) {
        vector<pair<mpz_class, mpz_class>> candidate_relations;
        for (int i = 0; i < rows; ++i) {
            if (history[r][i] == 1) candidate_relations.push_back(smooth_relations[i]);
        }
        
        Coordinate gamma_sq_algebraic = map_phi(1, 0);
        for (const auto& rel : candidate_relations) {
            gamma_sq_algebraic = tensor_multiply(gamma_sq_algebraic, map_phi(rel.first, rel.second), c_val);
        }
        
        Coordinate dummy_root;
        if (is_algebraic_square(gamma_sq_algebraic, c_val, dummy_root)) {
            final_relations = candidate_relations;
            break;
        }
    }
    
    if (final_relations.empty()) {
        cout << "Failed to find a true algebraic square in the null space.\n";
        return 1;
    }
    
    cout << "Dependency verified. Collapsing " << final_relations.size() << " relations...\n";

    mpz_class U_sq_rational = 1;
    Coordinate gamma_sq_algebraic = map_phi(1, 0);
    for (const auto& rel : final_relations) {
        U_sq_rational *= (rel.first - rel.second * m);
        gamma_sq_algebraic = tensor_multiply(gamma_sq_algebraic, map_phi(rel.first, rel.second), c_val);
    }
    
    mpz_class U = sqrt(abs(U_sq_rational)); 
    Coordinate gamma_root;
    is_algebraic_square(gamma_sq_algebraic, c_val, gamma_root);
    
    pair<mpz_class, mpz_class> gamma_ab = unmap_phi(gamma_root);
    mpz_class V = gamma_ab.first - gamma_ab.second * m;
    
    V = V % N; if (V < 0) V += N;
    U = U % N; if (U < 0) U += N;

    mpz_class factor1, factor2;
    mpz_class diff = abs(V - U), sum = V + U;
    mpz_gcd(factor1.get_mpz_t(), diff.get_mpz_t(), N.get_mpz_t());
    mpz_gcd(factor2.get_mpz_t(), sum.get_mpz_t(), N.get_mpz_t());

    cout << "\n[Final Factors]\n";
    cout << "Factor 1: " << factor1 << "\n";
    cout << "Factor 2: " << factor2 << "\n";
    
    return 0;
}