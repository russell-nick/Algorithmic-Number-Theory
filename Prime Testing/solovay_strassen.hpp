/**
 * @file solovay_strassen.hpp
 * @brief Header file for Solovay-Strassen Primality Testing.
 *
 * @author Nicholas Russell
 */

#ifndef SOLOVAY_STRASSEN_HPP
#define SOLOVAY_STRASSEN_HPP

/**
 * Solovay-Strassen Prime Test (k iterations)
 * Determine whether 'n' is prime with probability >= 1-(1/2)^k
 *
 * Complexity: O(k (lg n)^3) bit operations
 *
 * @param n number to prime test
 * @param k iterations
 * @return true if n is prime (probability >= 1-(1/2)^k), false if n is definitely composite
 */
bool solovay_strassen(long long int n, int k);

#endif /* SOLOVAY_STRASSEN_HPP */
