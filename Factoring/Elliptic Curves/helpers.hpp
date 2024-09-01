/**
 * @file helpers.hpp
 * @brief Header file for helper functions required for elliptic curve factoring.
 *
 * @author Nicholas Russell
 */

#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <boost/multiprecision/gmp.hpp>

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
std::vector<long long int> sieve_of_eratosthenes(long long int n);

/**
 * Compute floor(log_b(n)) (also index of most significant coefficient (bit when b = 2))
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param n number
 * @param b base
 * @return floor(log_b(n))
 */
long long int floor_log(mp::mpz_int n, mp::mpz_int b);

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
mp::mpz_int mod(mp::mpz_int a, mp::mpz_int n);

/**
 * Compute a^e (mod n) using binary exponentiation
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
mp::mpz_int mod_pow(mp::mpz_int a, mp::mpz_int e, mp::mpz_int n);


/**
 * Calculate the gcd of two integers: 'a' and 'b'.
 *
 * Bit Complexity: O((lg a)(lg b)) bit operations.
 *
 * @param a,b
 * @return gcd(a, b)
 */
mp::mpz_int gcd(mp::mpz_int a, mp::mpz_int b);

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
mp::mpz_int extended_gcd(mp::mpz_int a, mp::mpz_int b, mp::mpz_int& x, mp::mpz_int& y);

/**
 * Compute lcm(a_0, a_1, ..., a_n), where
 * nums = {a_0, a_1, ..., a_n}
 *
 * @param nums vector of numbers to take lcm of
 * @return lcm of all integers in 'nums'
 */
mp::mpz_int lcm(const std::vector<mp::mpz_int>& nums);

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
bool miller_rabin(mp::mpz_int n, int k);

#endif /* HELPERS_HPP */
