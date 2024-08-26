/**
 * @file binomial.cpp
 * @brief Header file for algorithms related to binomial coefficients.
 *
 * @author Nicholas Russell
 */

#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

/**
 * Calculate the binomial coefficient n choose k
 *
 * @param n
 * @param k
 * @return n choose k
 */
mp::cpp_int binomial_coefficient(mp::cpp_int n, mp::cpp_int k);

/**
 * Determine whether Q is a nontrivial binomial coefficient (with proof)
 * (time polynomial in log Q)
 *
 * Bit Complexity: O((lg Q)^3 * (lg sqrt(Q))^3) bit operations
 *
 * @param Q
 * @return true if Q is a nontrivial binomial coefficient (prints Q = n choose k), false otherwise.
 */
bool is_binomial_coefficient(mp::cpp_int Q);

/**
 * Find all (n, k) pairs such that Q = n choose k is a non-trivial binomial coefficient
 * (time polynomial in log Q)
 *
 * Bit Complexity: O((lg Q)^3 * (lg sqrt(Q))^3) bit operations
 *
 * @param Q
 * return vector containing all (n, k) pairs such that Q = n choose k (nontrivially)
 */
std::vector<std::pair<mp::cpp_int, mp::cpp_int>> all_binomial_coefficients(mp::cpp_int Q);

#endif /* BINOMIAL_HPP */
