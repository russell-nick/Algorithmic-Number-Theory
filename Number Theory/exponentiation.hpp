/**
 * @file exponentiation.hpp
 * @brief Header file for various exponentiation algorithms.
 *
 * @author Nicholas Russell
 */

#ifndef EXPONENTIATION_HPP
#define EXPONENTIATION_HPP

/**
 * Compute a^e using binary exponentiation (repeated squaring)
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
long long int bin_pow(long long int a, long long int e);

/**
 * Compute a^e using binary exponentiation (repeated squaring) recursively
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
long long int bin_pow_recursive(long long int a, long long int e);

/**
 * Compute a^e using m-ary exponentiation
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param m base (radix) of exponent to work with
 * @return a^e
 */
long long int m_pow(long long int a, long long int e, long long int m);

/**
 * Compute a^e using m-ary exponentiation (recursively)
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param m base (radix) of exponent to work with
 * @return a^e
 */
long long int m_pow_recursive(long long int a, long long int e, long long int m);

/**
 * Compute a^e using binary exponentiation (repeated squaring)
 * where the result is a floating point numnber
 * (a is floating point or e is negative)
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent (can be negative)
 * @return a^e
 */
double bin_pow(double a, long long int e);

/**
 * Compute a^e by repeatedly multiplying by 'a' e times
 *
 * Complexity: O(e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
long long int slow_pow(long long int a, long long int e);

#endif /* EXPONENTIATION_HPP */
