/**
 * @file miller_rabin.cpp
 * @brief Implementation of Miller-Rabin Primality Testing.
 *
 * @author Nicholas Russell
 */

#include <iostream>
#include <random>
#include <vector>

#include "miller_rabin.hpp"

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
 * Miller-Rabin Prime Test (k iterations)
 * Determine whether 'n' is prime with probability >= 1-(1/4)^k
 *
 * Complexity: O(k (lg n)^3) bit operations
 *
 * @param n number to prime test
 * @param k iterations
 * @return true if n is prime (probability >= 1-(1/4)^k), false if n is definitely composite
 */
bool miller_rabin(long long int n, int k) {
    // Simple checks first:
    if (n == 2 || n == 3) return true;
    // if n = 1 or n > 2 is even, n is composite (ignore negatives):
    if (n <= 1 || n % 2 == 0) return false;
    
    // Now we know that n > 3 is odd
    
    // Let n-1 = 2^r m, where m is odd (factor out all 2's from n-1):
    long long int r = 0, m = n - 1;
    while (m % 2 == 0) {
        r++;
        m /= 2;
    }
    
    // Perform miller-rabin prime testing for k iterations:
    for (int i = 0; i < k; i++) {
        
        // Generate 'a' uniformly at random from the interval [2, n-2]
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<long long int> distribution(2, n-2); // or [1, n-1] also works?
        long long int a = distribution(generator);
        
        // Now we want to calculate the sequence
        // a^m, a^{2m}, ..., a^{2^(r-1)m}, a^{2^r m} = a^{n-1} (mod n)
        // Do this by finding a^m and squaring this r - 1 times
        a = mod_pow(a, m, n);
        
        // If a^m = 1 or a^m = n-1 (mod n), then n is 'probably' prime
        if (a == 1 || a == n - 1) continue;
        
        bool probably_prime = false; // flag to skip to next k iteration if 'probably' prime
        for (int j = 0; j < r-1; j++) {
            a = (a * a) % n;
            
            // If a^{2^{j+1} m} = n-1 (mod n), then n is 'probably' prime
            if (a == n - 1) {
                probably_prime = true;
                break;
            }
        }
        
        // If n is 'probably' prime, continue to next iteration
        if (probably_prime) {
            continue;
        }
        
        // If we went through the sequence a^m, a^{2m}, ..., a^{2^(r-1)m}, a^{2^r m} = a^{n-1} (mod n)
        // without ever marking n as 'probably' prime, then n is composite
        return false;
    }

    // If we finish all iterations as 'probably' prime, then n is 'probably' prime
    double prime_probability = 1 - std::pow(0.25, k); // 1 - (1/4)^k
    
    std::cout << "n = " << n << " is prime with probability >= 1-(1/4)^" << k << " = " << prime_probability << std::endl;
    return true;
}

int main() {
    
    std::vector<long long int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 61543};
    std::vector<long long int> not_primes = {1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36};
    
    int k = 4; // num iterations
    std::cout << "Running Miller-Rabin with k = " << k << " iterations" << std::endl;
    
    std::cout << "\n********** Testing Prime Numbers **********" << std::endl;
    for (const auto& p : primes) {
        bool probably_prime = miller_rabin(p, k); // prints if p is prime
        if (!probably_prime) { // THIS SHOULD NEVER HAPPEN!!!
            std::cout << "Prime " << p << " incorrectly classified as composite" << std::endl;
        }
    }
    
    std::cout << "\n********** Testing Non-Prime Numbers **********" << std::endl;
    for (const auto& num : not_primes) {
        bool probably_prime = miller_rabin(num, k);
        if (probably_prime) {
            double failure_probability = std::pow(0.25, k);
            std::cout << "Composite " << num << " incorrectly classified as prime (probability: <= "
                      << "(1/4)^" << k << " = " << failure_probability << ")" << std::endl;
        }
        else {
            std::cout << "Composite " << num << " correctly classified as composite" << std::endl;
        }
    }
    
    return 0;
}
