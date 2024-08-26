/**
 * @file logs.hpp
 * @brief Header file for algorithms related to logarithms.
 *
 * @author Nicholas Russell
 */

#ifndef LOGS_HPP
#define LOGS_HPP

/**
 * Compute floor(log_b(n)) (also index of most significant coefficient (bit when b = 2))
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param n number
 * @param b base
 * @return floor(log_b(n))
 */
long long int floor_log(long long int n, long long int b);

/**
 * Convert n to base b:
 * n = a_k*b^k + ... + a_1*b + a_0
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param n number
 * @param b base
 * @return coefficients a_k, ..., a_1, a_0 in reverse order
 */
std::vector<long long int> base_b_conversion(long long int n, long long int b);

#endif /* LOGS_HPP */
