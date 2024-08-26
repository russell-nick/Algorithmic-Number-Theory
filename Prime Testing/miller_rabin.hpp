/**
 * @file miller_rabin.hpp
 * @brief Header file for Miller-Rabin Primality Testing.
 *
 * @author Nicholas Russell
 */

#ifndef MILLER_RABIN_HPP
#define MILLER_RABIN_HPP

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
bool miller_rabin(long long int n, int k);

#endif /* MILLER_RABIN_HPP */
