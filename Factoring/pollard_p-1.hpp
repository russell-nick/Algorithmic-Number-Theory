/**
 * @file pollard_p-1.hpp
 * @brief Header file for Pollard's p-1 Factorization Method.
 *
 * @author Nicholas Russell
 */

#ifndef POLLARD_P_1_HPP
#define POLLARD_P_1_HPP

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
long long int compute_E(long long int B);

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
long long int pollard_p_1(long long int n);

#endif /* POLLARD_P_1_HPP */
