/**
 * @file elliptic_curve_factoring.hpp
 * @brief Header file for Lenstra's Elliptic Curve Factorization Method.
 *
 * @author Nicholas Russell
 */

#ifndef ELLIPTIC_CURVE_FACTORING_HPP
#define ELLIPTIC_CURVE_FACTORING_HPP

#include <utility>  // pair
#include <boost/multiprecision/cpp_int.hpp>

#include "elliptic_curve.hpp"

namespace mp = boost::multiprecision;


/**
 * Try to compute a multiple of the order of a group related to the factoring algorithm:
 * Here we want to consider |E(Z/pZ)| for some prime p dividing n
 * If E = k|E(Z/pZ)|, then E*P = (0:1:0) and we can extract a non-trivial divisor
 *
 * We use E = lcm(1...B) = prod_{prime p <= B} p^(floor(log B / log p))
 *                       = prod_{prime p <= B} p^(floor(log_p(B))
 *
 * One can also use E = prod_{prime p <= B} p^(floor(log n / log p))
 *                    = prod_{prime p <= B} p^(floor(log_p(n))
 * for better odds of factoring, although E will be much larger
 * and factoring will take a bit longer. lcm(1...B) is a standard choice.
 *
 * @param B smoothness bound
 * @return lcm(1...B) = prod_{prime p <= B} p^(floor(log B / log p))
 *
 */
mp::cpp_int compute_E(const mp::cpp_int& B);

/**
 * "Psuedo-addition" from Lenstra's Elliptic Curve Factorization paper.
 * This is essentially applying mutliplication (repeated squaring) until the
 * group operation fails and stopping early if we find a non-trivial divisor.
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param P first point in addition
 * @param Q second point in addition
 * @return (R, d) pair, where R = P+Q (if no nontrivial divisor is found),
 *                      d is the nontrivial divisor if one is found (if R = inf),
 *                      or d = 0 if no nontrivial divisor is found (R != inf)
 *
 */
std::pair<ECPoint, mp::cpp_int> partial_addition(const ECPoint& P, const ECPoint& Q);

/**
 * "Psuedo-multiplication" from Lenstra's Elliptic Curve Factorization paper
 * (terminate early if we get a nontrivial divisor of n)
 *
 * Complexity: O(lg k) group operations
 *             O((lg k)(lg n)^2) bit operations
 *
 * @param k scalar
 * @param P point on elliptic curve
 * @return (R, d) pair, where R = k*P (if no nontrivial divisor is found: d = 1 or d = n)
 *                      or d is a nontrivial divisor (R != k*P since we exit early)
 *
 */
std::pair<ECPoint, mp::cpp_int> partial_multiplication(mp::cpp_int k, ECPoint P);

/**
 * Lenstra's Elliptic Curve Factoring (to split n)
 *
 * We use the smoothness bound B = L[1/2, 1/sqrt(2)] = exp( sqrt(1/2 * ln(n) * ln(ln(n))) )
 * as mentioned in Lenstra's paper. This is similar to the choice of B for factoring
 * algorithms such as Dixon's algorithm, the Quadratic Sieve, Pollard's p-1 method, etc.
 *
 * Complexity:
 * Using repeated squaring to calculate E*P uses
 * O(log E) = O(B) group operations for a total of O(B (lg n)^2) bit operations to split n.
 *
 * The (conjectured!) expected bit operations to find a factor of n is O(B (lg n)^2) * # Curves to check
 * which works out to be O( exp(sqrt( (2 + o(1)) * ln p ln ln p )) * (lg n)^2) bit operations,
 * where p is the smallest prime divisor of n.
 *
 * @param n number to factor
 * @param B smoothness bound
 * @return a non-trivial divisor of n
 *
 */
mp::cpp_int lenstra_elliptic_curve(const mp::cpp_int& n, const mp::cpp_int& B, const mp::cpp_int& E);

/**
 * Merge the prime decompositions of two non-trivial divisors of n
 *
 * @param f1,f2 prime decompositions of two non-trivial divisors of n
 * @return prime decomposition of f1*f2 as vector of pairs (p_1, e_1), ..., (p_r, e_r),
 *         where f1*f2 = p_1^{e_1} * ... * p_r^{e_r} is the prime decomposition of f1*f2
 */
std::vector<std::pair<mp::cpp_int, mp::cpp_int>> merge_factors(
    const std::vector<std::pair<mp::cpp_int, mp::cpp_int>>& f1,
    const std::vector<std::pair<mp::cpp_int, mp::cpp_int>>& f2);

/**
 * Full factorization using Lenstra's Elliptic Curve method.
 *
 * Before using elliptic curves, we use the following quick checks:
 * (1) Test if the number is prime
 * (2) Perform trial division with small primes
 *
 * If an elliptic curve does not produce a non-trivial divisor, try another one.
 * Once a sufficient number of curves have been tested, increase the smoothness bound.
 * Repeat until a non-trivial divisor d is found and recurse on d and n/d
 *
 * @param n number to factor
 * @return returns vector of pairs (p_1, e_1), ..., (p_r, e_r),
 *         where n = p_1^{e_1} * ... * p_r^{e_r} is the prime decomposition of n
 */
std::vector<std::pair<mp::cpp_int, mp::cpp_int>> factor(mp::cpp_int n);

/**
 * Prints the prime decomposition of a factored number
 *
 * @param n number to factor
 * @param factors vector of pairs (p_1, e_1), ..., (p_r, e_r),
 *        where n = p_1^{e_1} * ... * p_r^{e_r} is the prime decomposition of n
 */
void print_factors(const mp::cpp_int& n, const std::vector<std::pair<mp::cpp_int, mp::cpp_int>>& factors);

#endif /* ELLIPTIC_CURVE_FACTORING_HPP */
