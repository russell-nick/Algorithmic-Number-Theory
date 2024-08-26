/**
 * @file fermat.cpp
 * @brief Implementation of Fermat Primality Testing.
 *
 * @author Nicholas Russell
 */

#include <iostream>
#include <random>
#include <vector>

#include "fermat.hpp"


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
 * Fermat Prime Test (k iterations)
 * Determine whether 'n' is prime or Carmichael with probability >= 1-(1/2)^k
 *
 * Complexity: O(k (lg n)^3) bit operations
 * (dominated by calculating a^(n-1) mod n with repeated squaring: lg n)
 *
 * @param n number to prime test
 * @param k iterations
 * @return true if n is prime or Carmichael (probability >= 1-(1/2)^k), false if n is definitely composite
 */
bool fermat(long long int n, int k) {
    // Simple checks first:
    if (n == 2 || n == 3) return true;
    // if n = 1 or n > 2 is even, n is composite (ignore negatives):
    if (n <= 1 || n % 2 == 0) return false;
    
    // Perform fermat prime testing for k iterations:
    for (int i = 0; i < k; i++) {
        
        // Generate 'a' uniformly at random from the interval [2, n-2]
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<long long int> distribution(2, n-2);
        long long int a = distribution(generator);
        
        // Exploit Fermat's little theorem and compute a^(n-1) mod n.
        
        // If a^(n-1) != 1 (mod n), then n is composite (by contrapositive of Fermat's little theorem)
        a = mod_pow(a, n-1, n);
        
        if (a != 1) return false; // n is composite
    }
    
    // If we reach this point, n is 'probably' prime (or Carmichael!)
    double prime_probability = 1 - std::pow(0.5, k); // 1 - (1/2)^k
    
    std::cout << "n = " << n << " is prime or Carmichael with probability >= 1-(1/2)^" << k << " = " << prime_probability << std::endl;
    return true;
}

int main() {
    
    std::vector<long long int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 61543};
    std::vector<long long int> not_primes = {1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36};
    
    int k = 4; // num iterations
    std::cout << "Running Fermat Prime Test with k = " << k << " iterations" << std::endl;
    
    std::cout << "\n********** Testing Prime Numbers **********" << std::endl;
    for (const auto& p : primes) {
        bool probably_prime = fermat(p, k); // prints if p is prime
        if (!probably_prime) { // THIS SHOULD NEVER HAPPEN!!!
            std::cout << "Prime " << p << " incorrectly classified as composite" << std::endl;
        }
    }
    
    std::cout << "\n********** Testing Non-Prime Numbers **********" << std::endl;
    for (const auto& num : not_primes) {
        bool probably_prime = fermat(num, k);
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
