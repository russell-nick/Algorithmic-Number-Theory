/**
 * @file gcd_lcm.hpp
 * @brief Header file for various GCD and LCM algorithms.
 *
 * @author Nicholas Russell
 */

#ifndef GCD_LCM_HPP
#define GCD_LCM_HPP

/**
 * Calculate the gcd of two integers: 'a' and 'b'.
 *
 * Complexity: Lame (~1844) showed that the Euclidean algorithm
 * uses <= 5log_10(min(a, b)) steps, so the runtime is O(log(min(a, b))).
 *
 * Bit Complexity: O((lg a)(lg b)) bit operations with Collins' bound.
 *
 * @param a,b
 * @return gcd(a, b)
 */
long long int gcd(long long int a, long long int b);

/**
 * Recursively calculate the gcd of two integers: 'a' and 'b'.
 *
 * Complexity: O(log(min(a, b))) steps
 * Bit Complexity: O((lg a)(lg b)) bit operations
 *
 * @param a,b
 * @return gcd(a, b)
 */
long long int gcd_recursive(long long int a, long long int b);

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
long long int extended_gcd(long long int a, long long int b, long long int& x, long long int& y);

/**
 * Compute lcm(a, b) = ab / gcd(a, b)
 *
 * Complexity (same as gcd): O(log(min(a, b))) divisions
 * Bit Complexity (same as gcd): O((lg a)(lg b)) bit operations.
 *
 * @param a,b
 * @return lcm(a, b)
 */
long long int lcm(long long int a, long long int b);

/**
 * Compute gcd(a_0, a_1, ..., a_n), where
 * nums = {a_0, a_1, ..., a_n}
 *
 * @param nums vector of numbers to take gcd of
 * @return gcd of all integers in 'nums'
 */
long long int gcd(std::vector<long long int>& nums);

/**
 * Compute lcm(a_0, a_1, ..., a_n), where
 * nums = {a_0, a_1, ..., a_n}
 *
 * @param nums vector of numbers to take lcm of
 * @return lcm of all integers in 'nums'
 */
long long int lcm(std::vector<long long int>& nums);

#endif /* GCD_LCM_HPP */
