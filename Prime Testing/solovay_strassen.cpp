/**
 * @file solovay_strassen.cpp
 * @brief Implementation of Solovay-Strassen Primality Testing.
 *
 * @author Nicholas Russell
 */

#include <iostream>
#include <random>
#include <vector>

#include "solovay_strassen.hpp"

/**
 * Compute the Jacobi symbol (a|n):
 * If n = p_1^{e_1}*...*p_r^{e_r} (prime decomposition), then the Jacobi symbol (a|n)
 * is equal to the product of the Legendre symbols (a|p_1)^(e_1)*...* (a|p_r)^(e_r)
 *
 * The following algorithm is based off Bach and Shallit (Algorithmic Number Theory) Section 5.9
 *
 * Complexity (same as Euclid's gcd): O((lg a)(lg n)) bit operations
 *
 * @param a 'a' from (a|n)
 * @param n 'n' odd > 0 from (a|n)
 * @return (a|n): 1 if 'a' is a quadratic residue mod n,
 *               -1 if 'a' is not a quadratic residue mod n (and n doesn't divide a),
 *                0 if n | a
 */
int jacobi(long long int a, long long int n) {
    // First, reduce a modulo n
    a %= n;
    if (a < 0) a = (a + n) % n; // convert to least non-negative residue
    
    int t = 1;
    while (a != 0) {
        while (a % 2 == 0) {
            a = a / 2;
            long long int rem = n % 8;
            if (rem == 3 || rem == 5) t = -t;
        }
        std::swap(a, n);
        if (a % 4 == 3 && n % 4 == 3) t = -t;
        a %= n;
    }
    
    if (n == 1) return t;
    return 0;
}

/**
 * Compute a^e (mod n) using binary exponentiation
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
long long int mod_pow(long long int a, long long int e, long long int n) {
    a %= n;
    long long int result = 1;
    while (e > 0) {
        if (e & 1) {
            result = (result * a) % n;
        }
        a = (a * a) % n;
        e >>= 1;
    }
    return result;
}

/**
 * Solovay-Strassen Prime Test (k iterations)
 * Determine whether 'n' is prime with probability >= 1-(1/2)^k
 *
 * Complexity: O(k (lg n)^3) bit operations
 * (dominated by calculating a^((n-1)/2) mod n with repeated squaring: lg n)
 *
 * @param n number to prime test
 * @param k iterations
 * @return true if n is prime (probability >= 1-(1/2)^k), false if n is definitely composite
 */
bool solovay_strassen(long long int n, int k) {
    // Simple checks first:
    if (n == 2 || n == 3) return true;
    // if n = 1 or n > 2 is even, n is composite (ignore negatives):
    if (n <= 1 || n % 2 == 0) return false;
    
    // Now we know that n is odd and > 3, so (n-1)/2 is an integer >= 1:
    
    // Perform solovay-strassen prime testing for k iterations:
    for (int i = 0; i < k; i++) {
        
        // Generate 'a' uniformly at random from the interval [2, n-2]
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<long long int> distribution(2, n-2);
        long long int a = distribution(generator);
        
        // Compute the jacobi symbol (a / n). If (a / n) = 0 or
        // a^((n-1)/2) != (a / n) (mod n), then n is composite
        long long int jacobian = jacobi(a, n);
        if (jacobian == -1) jacobian = n - 1; // convert to least non-negative residue
        
        // If (a|n) = 0, then there exists some prime divisor p of n such that p | a
        if (jacobian == 0) return false; // n is composite
        
        // The following check uses the fact that if n is an odd prime, then
        // (a|n) = a^((n-1)/2) (mod n). By the contrapositive, we get that
        // if (a|n) != a^((n-1)/2) (mod n), then n is composite
        a = mod_pow(a, (n-1)/2, n);
        if (a != jacobian) return false; // n is composite
    }
    
    // If we reach this point, n is 'probably' prime
    double prime_probability = 1 - std::pow(0.5, k); // 1 - (1/2)^k
    
    std::cout << "n = " << n << " is prime with probability >= 1-(1/2)^" << k << " = " << prime_probability << std::endl;
    return true;
}

int main() {
    
//    std::cout << jacobi(1001, 9907) << std::endl;
    
    std::vector<long long int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 61543};
    std::vector<long long int> not_primes = {1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36};
    
    int k = 4; // num iterations
    std::cout << "Running Solovay-Strassen with k = " << k << " iterations" << std::endl;
    
    std::cout << "\n********** Testing Prime Numbers **********" << std::endl;
    for (const auto& p : primes) {
        bool probably_prime = solovay_strassen(p, k); // prints if p is prime
        if (!probably_prime) { // THIS SHOULD NEVER HAPPEN!!!
            std::cout << "Prime " << p << " incorrectly classified as composite" << std::endl;
        }
    }
    
    std::cout << "\n********** Testing Non-Prime Numbers **********" << std::endl;
    for (const auto& num : not_primes) {
        bool probably_prime = solovay_strassen(num, k);
        if (probably_prime) {
            double failure_probability = std::pow(0.5, k);
            std::cout << "Composite " << num << " incorrectly classified as prime (probability: <= "
                      << "(1/2)^" << k << " = " << failure_probability << ")" << std::endl;
        }
        else {
            std::cout << "Composite " << num << " correctly classified as composite" << std::endl;
        }
    }
    
    return 0;
}
