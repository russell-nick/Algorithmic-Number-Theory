/**
 * @file helpers.cpp
 * @brief Implementation of helper functions required for elliptic curve factoring.
 *
 * @author Nicholas Russell
 */

#include <cmath>
#include <iostream>
#include <random>

#include <boost/multiprecision/cpp_int.hpp>

#include "helpers.hpp"

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
long long int floor_log(mp::cpp_int n, mp::cpp_int b) {
    long long int result = 0;
    while (n >= b) {
        n /= b;
        result++;
    }
    return result;
}

/**
 * Calculate the least non-negative residue of 'a' modulo 'n'
 * C++'s "%" need not do this.
 * Ex: -1 % 5 = -1, but we want -1 % 5 = 4
 *
 * Uses the formula a (mod n) = a - n * floor(a / n)
 *
 * Bit Complexity: O((lg q)(lg n)) where q = floor(a / n)
 * (same complexity as multiplication / division)
 *
 * @param a number to reduce modulo 'n'
 * @param n modulus
 * @return a (mod n)
 */
mp::cpp_int mod(mp::cpp_int a, mp::cpp_int n) {
    if (n == 0) return a;
    mp::cpp_int result = a - n * (a/n); // quotient = a / n; return a - n * quotient
    return result < 0 ? result + n : result;
}

/**
 * Compute a^e (mod n) using binary exponentiation
 *
 * Complexity: O(lg e) multiplications
 *             O((lg e)(lg n)^2) bit operations
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
mp::cpp_int mod_pow(mp::cpp_int a, mp::cpp_int e, mp::cpp_int n) {
    a %= n;
    mp::cpp_int result = 1;
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
mp::cpp_int gcd(mp::cpp_int a, mp::cpp_int b) {
    mp::cpp_int temp;
    while (b > 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/**
 * Use the Extended Euclidean Algorithm to find gcd(a, b) and solve
 * ax + by = gcd(a, b) for x, y
 *
 * Complexity (same as gcd): O(log(min(a, b))) steps
 * Bit Complexity (same as gcd): O((lg a)(lg b)) bit operations
 *
 * @param[in] a,b
 * @param[out] x,y integers x,y such that ax + by = gcd(a, b)
 * @return gcd(a, b)
 */
mp::cpp_int extended_gcd(mp::cpp_int a, mp::cpp_int b, mp::cpp_int &x, mp::cpp_int &y) {
    x = 1;
    y = 0;
    
    if (b == 0) return abs(a);

    mp::cpp_int x1 = 0, y1 = 1;
    mp::cpp_int a1 = a, b1 = b;
    mp::cpp_int quotient, tmp;
    while (b1 != 0) {
        quotient = a1 / b1;

        tmp = a1;
        a1 = b1;
        b1 = tmp - quotient * b1;

        tmp = x;
        x = x1;
        x1 = tmp - quotient * x1;

        tmp = y;
        y = y1;
        y1 = tmp - quotient * y1;
    }
    return abs(a1);
}
//mp::cpp_int extended_gcd(mp::cpp_int a, mp::cpp_int b, mp::cpp_int& x, mp::cpp_int& y) {
//    if (b == 0) {
//        x = 1;
//        y = 0;
//        return abs(a);
//    }
//
//    mp::cpp_int x1, y1;
//    mp::cpp_int d = extended_gcd(b, a % b, x1, y1);
//
//    x = y1;
//    y = x1 - y1 * (a / b);
//    return d;
//}

/**
 * Compute lcm(a_0, a_1, ..., a_n), where
 * nums = {a_0, a_1, ..., a_n}
 *
 * @param nums vector of numbers to take lcm of
 * @return lcm of all integers in 'nums'
 */
mp::cpp_int lcm(const std::vector<mp::cpp_int>& nums) {
    if (nums.size() == 0) return 1;
    
    mp::cpp_int lcm = nums[0];
    for (int i = 1; i < nums.size(); i++) {
        lcm = lcm / gcd(lcm, nums[i]) * nums[i];
    }
    return lcm;
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
bool miller_rabin(mp::cpp_int n, int k) {
    // Simple checks first:
    if (n == 2 || n == 3) return true;
    // if n = 1 or n > 2 is even, n is composite (ignore negatives):
    if (n <= 1 || n % 2 == 0) return false;
    
    // Now we know that n > 3 is odd
    
    // Let n-1 = 2^r m, where m is odd (factor out all 2's from n-1):
    mp::cpp_int r = 0, m = n - 1;
    while (m % 2 == 0) {
        r++;
        m /= 2;
    }
    
    // Perform miller-rabin prime testing for k iterations:
    for (int i = 0; i < k; i++) {
        
        // Generate 'a' uniformly at random from the interval [2, n-2]
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<mp::cpp_int> distribution(2, n-2); // or [1, n-1] also works?
        mp::cpp_int a = distribution(generator);
        
        // Now we want to calculate the sequence
        // a^m, a^{2m}, ..., a^{2^(r-1)m}, a^{2^r m} = a^{n-1} (mod n)
        // Do this by finding a^m and squaring this r - 1 times
        a = mod_pow(a, m, n);
        
        // If a^m = 1 or a^m = n-1 (mod n), then n is 'probably' prime
        if (a == 1 || a == n - 1) continue;
        
        bool probably_prime = false; // flag to skip to next k iteration if 'probably' prime
        for (mp::cpp_int j = 0; j < r-1; j++) {
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
