/**
 * @file pollard_p-1.cpp
 * @brief Implementation of Pollard's p-1 Factorization Method.
 *
 * @author Nicholas Russell
 */

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include "pollard_p-1.hpp"

namespace mp = boost::multiprecision;

/**
 * Sieve of Eratosthenes to generate all primes <= a given limit
 *
 * Complexity: O(n(loglogn)) operations
 * Bit Complexity: O(n(logn)(loglogn)) bit operations
 *
 * @param n sieve limit
 * @return vector of all primes <= n
 */
std::vector<long long int> sieve_of_eratosthenes(long long int n) {
    std::vector<bool> is_prime(n+1, true);
    is_prime[0] = is_prime[1] = false;
    for (long long int i = 2; i * i <= n; i++) {
        if (is_prime[i]) {
            for (long long int j = i * i; j <= n; j += i)
                is_prime[j] = false;
        }
    }

    std::vector<long long int> primes;
    primes.reserve(is_prime.size());
    for (int i = 0; i < is_prime.size(); i++) {
        if (is_prime[i]) primes.push_back(i);
    }

    return primes;
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
 * Calculate the gcd of two integers: 'a' and 'b'.
 *
 * Bit Complexity: O((lg a)(lg b)) bit operations.
 *
 * @param a,b
 * @return gcd(a, b)
 */
long long int gcd(long long int a, long long int b) {
    long long int temp;
    while (b > 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/**
 * Generate a random number from Z_n^* uniformly at random
 *
 * Bit Complexity: The expected number of draws is
 *            n / phi(n) = O(loglogn), where phi(n) is Euler's totient function.
 *
 * @param n number specifying order of multiplicative group Z_n^*
 * @return a number chosen uniformly at random from Z_n^*
 */
long long int random_element_Zn_mult(long long int n) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<long long int> distribution(1, n-1);
    long long int x = distribution(generator);

    while (gcd(x, n) != 1) {
        x = distribution(generator);
    }

    return x;
}

/**
 * Compute floor(log_b(n)) (also index of most significant coefficient (bit when b = 2))
 *
 * To calculate log_b(n), we want to find x such that b^x = n. However, to find floor(logb n), we can
 * continuously divide n by b as long as n >= b. If we divide a total of m times, then 1 <= n/(b^m) < b
 * so log_b(b^m) <= log_b(n) < log_b(b^(m+1)) and m <= log_b(n) < m + 1. Thus, floor(log_b(n)) = m.
 *
 * Alternatively, consider the base-b representation n = a_m*b^m + ... + a_1*b + a_0, where
 * a_i in {0, 1, ..., m} and a_m > 0. Then b^(m+1) > m*b^m + ... + m*b + m >= n >= b^m so
 * taking log_b of b^(m+1) > n >= b^m gives m+1 > n >= m, so floor(log_b(n)) = m.
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param n number
 * @param b base
 * @return floor(log_b(n))
 */
long long int floor_log(long long int n, long long int b) {
    long long int result = 0;
    while (n >= b) {
        n /= b;
        result++;
    }
    return result;
}

/**
 * Try to compute a multiple of the order of a group related to the factoring algorithm:
 * Here we want to consider |Z/pZ| = p-1 for some prime p dividing n
 * If E = k(p-1), then x^E = x^(k(p-1)) = (x^(p-1))^k = 1^k = 1 (in Z/pZ)
 * and we have a chance to extract a non-trivial divisor with gcd(x^E - 1, n)
 * since x^E - 1 = 0 (mod p)
 *
 * We use E = prod_{prime p <= B} p^(floor(log B / log p))
 *          = prod_{prime p <= B} p^(floor(log_p(B)) = lcm(1...B)
 *
 * One can also use E = prod_{prime p <= B} p^(floor(log n / log p))
 *                    = prod_{prime p <= B} p^(floor(log_p(n))
 * for better odds of factoring, although E will be much larger
 * and factoring will take a bit longer.
 *
 * @param B smoothness bound
 * @return lcm(1...B) = prod_{prime p <= B} p^(floor(log B / log p))
 */
long long int compute_E(long long int B) {
    long long int E = 1;

    std::vector<long long int> primes = sieve_of_eratosthenes(B);
    for (auto p : primes) {
        long long int exponent = floor_log(B, p); // floor(log_p(B))
        E *= pow(p, exponent); // (exponent will be small, but can still replace pow with binary exp / loop)
    }

    return E;
}

/**
 * Pollard's p-1 Factoring (to split n)
 *
 * Complexity: Using repeated squaring to calculate x^E uses
 *             O(log E) = O(sum_{prime p <= B} (log B / log p) * log p)
 *                      = O((log B)(pi(B)) = O(B)
 *             group operations for a total of O(B (lg n)^2) bit operations to split n.
 *
 *  Note: Pollard's p-1 factorization method relies on the smoothness
 *  of the p-1's for each prime divisor p of n. If the p-1's consist of
 *  large prime factors, then this method will not be efficient and can
 *  devolve into trial divison.
 *
 * @param n number to factor
 * @return a non-trivial divisor of n
 */
long long int pollard_p_1(long long int n) {
    long long int B = pow(n, 0.25);
    long long int E = compute_E(B);

    // Generate 'x' uniformly at random from Z_n^*
    long long int x = random_element_Zn_mult(n);
    std::cout << "x: " << x << std::endl;
    std::cout << "B: " << B << std::endl;
    std::cout << "E: " << E << std::endl;
    
    // Compute x^E (mod n)
    x = mod_pow(x, E, n);
    std::cout << "x^E: " << x << std::endl;
    // Hope 1 < gcd(x^E - 1, n) < n to split n:
    long long int d = gcd(x - 1, n);
    
    // While d = n, try to compute x^E with a new x.
    // If we end up with d == 1, we should increase the smoothness bound B
    while (d == n) {
        std::cout << "Failed to find factor: gcd(x^E - 1, n) = n. Choosing new x." << std::endl;
        x = random_element_Zn_mult(n);
        std::cout << "x: " << x << std::endl;
        
        // Compute x^E (mod n)
        x = mod_pow(x, E, n);
        std::cout << "x^E: " << x << std::endl;
        // Hope 1 < gcd(x^E - 1, n) < n to split n:
        d = gcd(x - 1, n);
    }
    if (d == 1)
        std::cout << "Failed to find factor: gcd(x^E - 1, n) = 1. Try increasing B." << std::endl;
    
    return d;
}


int main() {
    long long int n = 1234567; // E overflows after this
    
    std::cout << "\n********** Factoring with Pollard's p-1 Method **********" << std::endl;
    long long int divisor = pollard_p_1(n);
    if (divisor != n && divisor != 1) {
        std::cout << "***** Pollard's p-1 Method successfully found a non-trivial divisor *****" << std::endl
                  << "n = " << divisor << " * " << (n / divisor) << std::endl;
    } else {
        std::cout << "***** Pollard's p-1 Method failed to find a non-trivial divisor *****" << std::endl
                  << "divisor = " << divisor << std::endl;
    }
    
    return 0;
}
