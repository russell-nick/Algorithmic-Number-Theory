/**
 * @file modular.hpp
 * @brief Header file for various algorithms involving modular arithmetic.
 * @author Nicholas Russell
 */

#ifndef MODULAR_HPP
#define MODULAR_HPP

/**
 * Calculate the least non-negative residue of 'a' modulo 'n'
 * C++'s "%" need not do this.
 * Ex: -1 % 5 = -1, but we want -1 % 5 = 4
 *
 * Uses the formula a (mod n) = a - n * floor(a / n)
 * and adds n to get a non-negative result if necessary
 *
 * Bit Complexity: O((lg q)(lg n)) where q = floor(a / n)
 * (same complexity as multiplication / division)
 *
 * @param a number to reduce modulo 'n'
 * @param n modulus
 * @return a (mod n)
 */
long long int mod(long long int a, long long int n);

/**
 * Compute a^e (mod n) using binary exponentiation (repeated squaring)
 *
 * Complexity: O(lg e) multiplications
 *             O((lg e)(lg n)^2) bit operations
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
long long int mod_pow(long long int a, long long int e, long long int n);

/**
 * Compute a^e by repeatedly multiplying by 'a' e times
 *
 * Complexity: O(e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
long long int slow_mod_pow(long long int a, long long int e, long long int n);

/**
 * Safely calculate a * b (mod n) using repeated squaring
 * (Note: this method can still overflow if a + b overflows)
 *
 * Complexity: O(log(min(a, b))) = O(log n) additions
 *
 * @param a
 * @param b
 * @param n modulus
 * @return a * b (mod n)
 */
long long int mod_mult(long long int a, long long int b, long long int n);

/**
 * Safely calculate a * b (mod n) by repeated addition in a loop
 * (Note: this method can still overflow if a + b overflows)
 *
 * Complexity: O(min(a, b)) = O(n) additions
 *
 * @param a
 * @param b
 * @param n modulus
 * @return a * b (mod n)
 */
long long int slow_mod_mult(long long int a, long long int b, long long int n);

/**
 * Compute the Jacobi symbol (a|n):
 * If n = p_1^{e_1}*...*p_r^{e_r} (prime decomposition), then the Jacobi symbol (a|n)
 * is equal to the product of the Legendre symbols (a|p_1)^(e_1)*...* (a|p_r)^(e_r)
 *
 * Complexity (same as Euclid's gcd): O((lg a)(lg n)) bit operations
 *
 * @param a 'a' from (a|n)
 * @param n 'n' odd > 0 from (a|n)
 * @return jacobi smybol (a|n): 1 if 'a' is a quadratic residue mod n (and n doesn't divide a),
 *                             -1 if 'a' is not a quadratic residue mod n (and n doesn't divide a),
 *                              0 if n | a
 */
int jacobi(long long int a, long long int n);

/**
 * Compute the Legendre Symbol using Euler's Criterion
 * (Determine the quadratic character of 'a' modulo 'p')
 *
 * Complexity: Repeated squaring uses O(lg p) multiplications
 * Bit Complexity: O((lg p)^3) bit operations
 *
 * @param a from (a|p)
 * @param p odd prime p from (a|p)
 * @return legendre symbol (a|p): 1 if a is a quadratic residue mod p (and p doesn't divide a),
 *                               -1 if a is a quadratic nonresidue mod p,
 *                                0 if p divides a
 */
int slow_legendre_symbol(long long int a, long long int p);

/**
 * Compute the Legendre Symbol using the standard Jacobi symbol algorithm
 * (Determine the quadratic character of 'a' modulo 'p')
 *
 * Note: When n is an odd prime, the Jacobi symbol (a|n) = Legendre symbol!
 *
 * Bit Complexity: O((lg a)(lg p)) = O((lg p)^2) bit operations
 *
 * @param a a from (a|p)
 * @param p p odd prime p from (a|p)
 * @return legendre symbol (a|p): 1 if a is a quadratic residue mod p (and p doesn't divide a),
 *                               -1 if a is a quadratic nonresidue mod p,
 *                                0 if p divides a
 */
int legendre_symbol(long long int a, long long int p);

/**
 * Determine whether 'a' is a quadratic residue modulo odd prime 'p'
 * (Whether there exists some x in Z_p with x^2 = a (mod p) )
 *
 * Bit Complexity: O((lg p)^2) bit operations
 *
 * @param a number to check for quadratic residue
 * @param p modulus (odd prime)
 * @return true if a is a quadratic residue mod p, false otherwise
 */
bool is_quadratic_residue(long long int a, long long int p);

/**
 * Find the inverse of 'a' modulo 'n'.
 *
 * Complexity (same as gcd): O(log(min(a, n))) divisions
 * Bit Complexity (same as gcd): O((lg a)(lg n)) bit operations.
 *
 * @param a
 * @param n modulus
 * @return inverse of a modulo n, 0 if none exists
 */
long long int mod_inverse(long long int a, long long int n);

/**
 * Solve the linear congruence ax = b (mod n) for x
 *
 * (There exists a solution for x if and only if gcd(a, n) divides b.)
 * If there exists some solution, there are exactly d unique solutions modulo n given by
 * x_0, x_0 + n/d, ..., d_0 + ((d-1)n)/d, where x_0 is a particular solution
 *
 * @param a,b
 * @param n modulus
 * @return vector of all x values (unique modulo n) such that ax = b (mod n) (empty if none)
 */
std::vector<long long int> linear_congruence_solver(long long int a, long long int b, long long int n);

/**
 * Solve the linear diophantine equation ax + by = c for integers x,y
 *
 * (There exists a solution if and only if d := gcd(a, b) divides c.)
 * If there exists some solution, there are infinitely many solutions:
 * x_0 + b/d*t, y_0 - a/d*t, where (x_0, y_0) is a particular solution
 * and t is an integer.
 *
 * @param a,b,c
 * @return some (x, y) pair such that ax + by = c, empty pair if no solution
 */
std::pair<long long int, long long int> linear_diophantine_solver(long long int a, long long int b, long long int c);

/**
 * Solve a system of linear congruences of the following form for x, y:
 * ax + by = r (mod n)
 * cx + dy = s (mod n)
 *
 * There exists a unique solution modulo n if and only if gcd(ad-bc, n) = 1.
 *
 * @param a,b,r,c,d,s
 * @param n modulus
 * @return unique (x, y) pair that solves the system, empty pair if no solution
 */
std::pair<long long int, long long int> linear_congruence_system_solver(
    long long int a, long long int b, long long int r,
    long long int c, long long int d, long long int s,
    long long int n);

#endif /* MODULAR_HPP */
