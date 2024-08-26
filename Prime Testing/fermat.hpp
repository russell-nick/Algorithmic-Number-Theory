/**
 * @file fermat.hpp
 * @brief Header file for Fermat Primality Testing.
 *
 * @author Nicholas Russell
 */

#ifndef FERMAT_HPP
#define FERMAT_HPP

/**
 * Fermat Prime Test (k iterations)
 * Determine whether 'n' is 'probably' prime or Carmichael
 *
 * Complexity: O(k (lg n)^3) bit operations
 *
 * @param n number to prime test
 * @param k iterations
 * @return true if n is prime or Carmichael (probability >= 1-(1/2)^k), false if n is definitely composite
 */
bool fermat(long long int n, int k);

#endif /* FERMAT_HPP */
