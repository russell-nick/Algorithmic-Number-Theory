/**
 * @file fibonacci.hpp
 * @brief Header file for various algorithms to compute fibonacci numbers.
 *
 * @author Nicholas Russell
 */

#ifndef FIBONACCI_HPP
#define FIBONACCI_HPP

/**
 * Compute the n-th fibonacci number naively using its recursive formula:
 * f(n) = f(n-1) + f(n-2)
 * f(0) = 1
 * f(1) = 1
 *
 * Complexity: O(2^n) operations (additions)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_rec(long long int n);

/**
 * Compute the n-th fibonacci number with dynamic programming
 * using its recursive formula:
 * f(n) = f(n-1) + f(n-2)
 * f(0) = 1
 * f(1) = 1
 *
 * Complexity: O(n) operations (additions)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci(long long int n);

/**
 * Compute the n-th fibonacci number with dynamic programming
 * using its recursive formula (with O(1) space):
 * f(n) = f(n-1) + f(n-2)
 * f(0) = 1
 * f(1) = 1
 *
 * Complexity: O(n) operations (additions)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_const_space(long long int n);

/**
 * Compute the n-th fibonacci number with matrix exponentation
 * using the following well-known identity:
 * [f_n+1 f_n] = [1  1] ^n
 * [f_n  f_n-1]    [1  0]
 * which is found by induction after observing that
 * [f_n+1 f_n] = [f_n    f_n-1] * [1  1]
 * [f_n  f_n-1]    [f_n-1 f_n-2]   [1  0]
 *
 * To compute A^n for some matrix A with integer coefficients,
 * we can use repeated squaring in the general linear group GLn(Z).
 * @see exponentation.cpp for more details
 *
 * Complexity: O(log n) operations (additions)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_mat_exp(long long int n);

/**
 * Compute the n-th fibonacci number with Binet's formula:
 * f(n) = (phi^n - psi^n) / sqrt(5), where
 * phi = (1 + sqrt(5)) / 2
 * psi = (1 - sqrt(5)) / 2
 *
 * Note: Due to floating point errors, evaluating this formula directly will
 * result in incorrect fibonacci numbers as n grows.
 *
 * Complexity: O(log n) operations (multiplications for exponents including n)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_binet(long long int n);

/**
 * Compute the n-th fibonacci number with Binet's formula (exact arithmetic):
 * f(n) = (phi^n - psi^n) / sqrt(5), where
 * phi = (1 + sqrt(5)) / 2
 * psi = (1 - sqrt(5)) / 2
 *
 * Note: To get exact arithmetic, we do computations in the extension field Q(sqrt(5))
 * to avoid floating point errors from using sqrt(5) in Binet's formula.
 *
 * Let x = a + b*sqrt(5) and y = c + d*sqrt(5) be two elements in Q(sqrt(5)).
 * Then x + y = (a + c) + (b + d)*sqrt(5)
 *     x – y = (a – c) + (b – d)*sqrt(5)
 *     x * y = (a + b*sqrt(5)) * (c + d*sqrt(5))
 *         = ac + ad*sqrt(5) + bc*sqrt(5) + 5bd
 *         = (ac + 5bd) + (ad + bc)*sqrt(5)
 *
 * Notice that f(n) = (phi^n - psi^n) / sqrt(5), so phi^n - psi^n = 0 + f(n)*sqrt(5).
 * Therefore, we can compute phi^n - psi^n in Q(sqrt(5)) to retrieve f(n).
 *
 * Complexity: O(log n) operations (multiplications for exponents including n)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_binet_exact(long long int n);

#endif /* FIBONACCI_HPP */
