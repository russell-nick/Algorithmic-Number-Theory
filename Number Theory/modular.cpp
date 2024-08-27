/**
 * @file modular.cpp
 * @brief Implementation of various algorithms involving modular arithmetic.
 *
 * @author Nicholas Russell
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "modular.hpp"
#include "gcd_lcm.hpp"

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
long long int mod(long long int a, long long int n) {
    if (n == 0) return a;
    long long int quotient = a / n;
    return a < 0 ? a - n * quotient + n : a - n * quotient;
}

/*
 * Calculate the least non-negative residue of 'a' modulo 'n'
 * (alternate version of above using built in %)
 */
// long long int mod(long long int a, long long int n) { return (a % n + n) % n; }


/**
 * Compute a^e (mod n) using binary exponentiation (repeated squaring)
 *
 * Let b_k*2^k + ... + b_1*2 + b_0 be the binary representation of
 * the exponent e. We use the fact that a^(b + c) = a^b * a^c
 * to calculate a^e using O(log(e)) multiplications:
    
 * a^e = a^(b_k*2^k + ... + b_1*2 + b_0) = a^(b_k * 2^k)*...*a^(b_1 * 2)*a^(b_0)
 * = a^{(b_k)(2^k)}*...*a^{(b_1)(2)}*a^(b_0)
 * so we can repeatedly square 'a' k-1 times and
 * multiply the result by 'a' if the current bit b_i = 1
 *
 * Complexity: O(lg e) multiplications
 *             O((lg e)(lg n)^2) bit operations
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
long long int mod_pow(long long int a, long long int e, long long int n) {
    // Reduce a modulo n to ensure 0 <= a < n:
    a %= n;
    // Note: If a can be negative, use one of the following to get least non-negative residue
    //       if you want to guarantee the output to be non-negative between 0 and n
    // if (a < 0) a = (a + n) % n;
    // a = mod(a, n) // instead of a %= n;
    
    long long int result = 1;
    while (e > 0) {
        // if current bit in binary rep is = 1
        if (e & 1) { // or e % 2 == 1
            result = (result * a) % n;
        }
        a = (a * a) % n;
        // shift to next bit in exponent:
        e >>= 1; // or e /= 2
    }
    return result;
}

/**
 * Compute a^e (mod n) using binary exponentiation (repeated squaring) recursively
 * using the formula:
 * a^e = a * (a^((e-1)/2))^2 if e is odd,
 * a^e = (a^(e/2))^2 if e is even
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @return a^e (mod n)
 */
long long int mod_pow_recursive(long long int a, long long int e, long long int n) {
    if (e == 0) return 1;
    if (e % 2 == 0) { // e even
        long long int t = mod_pow_recursive(a, e / 2, n);
        return (t * t) % n;
    } else { // e odd
//        long long int t = mod_pow_recursive(a, e - 1, n);
//        return (a * t) % n;
        long long int t = mod_pow_recursive(a, (e - 1) / 2, n);
        return (a * t * t) % n;
    }
    return 0;
}

/**
 * Compute a^e (mod n) using m-ary exponentiation (recursively) using the formula:
 * a^e = a^(e mod m) * ( a^((e - (e mod m)) / m) )^m if e != 0 (mod m)
 * a^e = (a^(e/m))^m if e = 0 (mod m)
 *
 * (If dealing with very large numbers (> long long), one can also precompute
 *  and use a^0, a^1, a^2, ..., a^(m-1) to avoid the use of using something
 *  like binary exponentation for the a^(e mod m) calculation)
 *
 * Note: It is best to let m be a power of 2.
 *
 * Complexity: <= m + log(e) + log_m(e) = O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @param m base (radix) of exponent
 * @return a^e
 */
long long int mod_m_pow_recursive(long long int a, long long int e, long long int n, long long int m) {
    if (e == 0) return 1;
    if (e % m == 0) {
        long long int t = mod_m_pow_recursive(a, e / m, n, m);
        return mod_pow(t, m, n);
    } else {
        long long int rem = e % m;
        long long int t = mod_m_pow_recursive(a, (e - rem) / m, n, m);
        return (mod_pow(a, rem, n) * mod_pow(t, m, n)) % n;
    }
    return 0;
}

/**
 * Compute a^e (mod n) using m-ary exponentiation
 *
 * Let b_k*m^k + ... + b_1*m + b_0 be the m-ary representation of
 * the exponent e. This method is very similar to binary exponentiation,
 * but working in base m instead of base 2.
 *
 * a^e = a^(b_k*m^k + ... + b_1*m + b_0) = a^(b_k * m^k)*...*a^(b_1 * m)*a^(b_0)
 * = a^{(b_k)(m^k)}*...*a^{(b_1)(m)}*a^(b_0)
 * so we can repeatedly raise 'a' to the m-th power k-1 times and
 * multiply the result by 'a^m' if the current coefficient b_i != 0
 *
 *  (If dealing with very large numbers, one can also precompute
 *  and use a^0, a^1, a^2, ..., a^(m-1) to avoid the use of using something
 *  like binary exponentation for the a^(e mod m) calculation)
 *
 * It is best to let m be a power of 2.
 *
 * Complexity: <= m + log(e) + log_m(e) = O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param n modulus
 * @param m base (radix) of exponent to work with
 * @return a^e (mod n)
 */
long long int mod_m_pow(long long int a, long long int e, long long int n, long long int m) {
    a %= n;
    
    // Precompute a^2, a^3, ..., a^m and use if doing left-to-right exponentiation
    // std::vector<long long int> a_pows(m);
//    a_pows[0] = 1; a_pows[1] = a;
//    for (int i = 2; i <= m; i++) {
//        a_pows[i] = (a * a_pows[i-1]) % n;
//    }
    
    long long int result = 1;
    while (e > 0) {
        // If the current base-m coefficient is nonzero (analogous to current bit is 1)
        long long int current_coeff = e % m;
        if (current_coeff != 0) {
            result = (result * mod_pow(a, current_coeff, n)) % n;
        }
        // Compute a^m (mod n) (analogous to squaring with binary exponentiation base m = 2)
        a = mod_pow(a, m, n);
        // Shift to next base-m coefficient in exponent (analogous to shifting to next bit):
        e /= m;
    }
    return result;
}

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
long long int slow_mod_pow(long long int a, long long int e, long long int n) {
    a %= n;
    long long int result = 1;
    for (int i = 0; i < e; i++) {
        result = (result * a) % n;
    }
    return result;
}

/**
 * Safely calculate a * b (mod n) using repeated squaring
 *
 * To multiply two integers a and b, we can use repeated squaring
 * to calculate a^b (or b^a) in the additive group (Z_n, +). Since
 * the operation is commutative, we can take the "exponent" to be
 * min(a, b) and only do O(log(min(a, b))) = O(log n) additions.
 *
 * (Note: this method can still overflow if a + b overflows. To absolutely
 * guarantee safe multiplication, one can instead cast the 64 bit integers
 * to 128 bit integers and directly compute (a*b) % n)
 *
 * @param a
 * @param b
 * @param n modulus
 * @return a * b (mod n)
 */
long long int mod_mult(long long int a, long long int b, long long int n) {
    a %= n;
    b %= n;
    if (a < b) {
        std::swap(a, b);
    }
    long long int result = 0;
    while (b > 0) {
        if (b & 1) {
            result = (result + a) % n;
        }
        a = (a + a) % n;
        b >>= 1;
    }
    return result;
}

/**
 * Safely calculate a * b (mod n) by repeated addition in a loop
 * (Note: as above, this method can still overflow if a + b overflows)
 *
 * Complexity: O(min(a, b)) = O(n) additions
 *
 * @param a
 * @param b
 * @param n modulus
 * @return a * b (mod n)
 */
long long int slow_mod_mult(long long int a, long long int b, long long int n) {
    a %= n;
    b %= n;
    if (a < b) {
        std::swap(a, b);
    }
    long long int result = 0;
    for (int i = 0; i < b; i++) {
        result = (result + a) % n;
    }
    return result;
}

/**
 * Compute the Jacobi symbol (a|n):
 * If n = p_1^{e_1}*...*p_r^{e_r} (prime decomposition), then the Jacobi symbol (a|n)
 * is equal to the product of the Legendre symbols (a|p_1)^(e_1)*...* (a|p_r)^(e_r)
 *
 * The following algorithm is based off Bach and Shallit (Algorithmic Number Theory) Section 5.9
 *
 * Complexity (same as Euclid's gcd): O((lg a)(lg n)) bit operations
 *
 * @param a 'a' from (a|n)
 * @param n 'n' odd > 0 from (a|n)
 * @return jacobi smybol (a|n): 1 if 'a' is a quadratic residue mod n (and n doesn't divide a),
 *                             -1 if 'a' is not a quadratic residue mod n (and n doesn't divide a),
 *                              0 if n | a
 */
int jacobi(long long int a, long long int n) {
    // First, reduce a modulo n
    a %= n;
    if (a < 0) a = (a + n) % n; // convert to least non-negative residue
    
    int t = 1;
    while (a != 0) {
        while (a % 2 == 0) {
            a = a / 2;
            long long int rem = n % 8;
            if (rem == 3 || rem == 5) t = -t;
        }
        std::swap(a, n);
        if (a % 4 == 3 && n % 4 == 3) t = -t;
        a %= n;
    }
    
    if (n == 1) return t;
    return 0;
}


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
int slow_legendre_symbol(long long int a, long long int p) {
    return static_cast<int>( mod_pow(a, (p - 1)/2, p) );
}

/**
 * Compute the Legendre Symbol using the standard Jacobi symbol algorithm
 * (Determine the quadratic character of 'a' modulo 'p')
 *
 * Note: When n is an odd prime, the Jacobi symbol (a|n) = Legendre symbol!
 *
 * To calculate the legendre symbol (a|p) for an odd prime p, this algorithm
 *       is faster than Euler's criterion which would calculate a^((p-1)/2)
 *       with O(lg p) multiplications for a total of O((lg p)^3) bit operations
 *       (when n is an odd prime, jacobi = legendre symbol)
 *
 * Bit Complexity: O((lg a)(lg p)) = O((lg p)^2) bit operations
 *
 * @param a a from (a|p)
 * @param p p odd prime p from (a|p)
 * @return legendre symbol (a|p): 1 if a is a quadratic residue mod p (and p doesn't divide a),
 *                               -1 if a is a quadratic nonresidue mod p,
 *                                0 if p divides a
 */
int legendre_symbol(long long int a, long long int p) {
    return jacobi(a, p);
}

/**
 * Determine whether 'a' is a quadratic residue modulo odd prime 'p'
 * (Whether there exists some x in Z_p with x^2 = a (mod p) )
 *
 * Bit Complexity: O((lg p)^3) bit operations
 *
 * @param a number to check for quadratic residue
 * @param p modulus (odd prime)
 * @return true if a is a quadratic residue mod p, false otherwise
 */
bool is_quadratic_residue(long long int a, long long int p) {
    return legendre_symbol(a, p) == 1;
}

/**
 * Find the inverse of 'a' modulo 'n'.
 * There exists a solution if and only if gcd(a, n) = 1.
 * If there is no solution, 0 is returned
 *
 * Complexity (same as gcd): O(log(min(a, n))) divisions
 * Bit Complexity (same as gcd): O((lg a)(lg n)) bit operations.
 *
 * @param a
 * @param n modulus
 * @return inverse of a modulo n, 0 if none exists
 */
long long int mod_inverse(long long int a, long long int n) {
    long long int x, y;
    long long int d = extended_gcd(a, n, x, y); // ax = 1 mod n -> ax - ny = 1
    if (d != 1) return 0;
    return x % n; // OR mod(x, n) for guaranteed positive result instead of x % n
}

/**
 * Solve the linear congruence ax = b (mod n) for x
 *
 * There exists a solution for x if and only if d := gcd(a, n) divides b.
 * If so, there are exactly d unique solutions modulo n given by
 * x_0, x_0 + n/d, ..., d_0 + ((d-1)n)/d, where x_0 is a particular solution
 *
 * If gcd(a, n) = 1, then x_0 = a^-1*b is the unique solution, but this doesn't
 * work in general. Solving ax = b (mod n) is the same as ax - ny = b for some integer y
 *
 * Use the Euclidean algorithm to find x,y such that ax + ny = gcd(a, n) = d
 * Multiply by b/d to get a(bx/d) + n(by/d) = b
 * Therefore, x_0 = bx/d is a particular solution
 *
 * @param a,b
 * @param n modulus
 * @return vector of all x values (unique modulo n) such that ax = b (mod n) (empty if none)
 */
std::vector<long long int> linear_congruence_solver(long long int a, long long int b, long long int n) {
    long long int x, y;
    long long d = extended_gcd(a, n, x, y);
    // There exists a solution if and only if gcd(a, n) divides b:
    if (b % d != 0) return {};
    std::vector<long long int> solutions;
    
    x = mod(x*(b/d), n); // x = (x * (b/d)) % n;
    for (int i = 0; i < d; i++) {
        solutions.push_back(mod(x + i*(n/d), n)); // O(d ((lg x + n)(lg n)) ) -> O(d(lg n)^2) ?
    }
    
    return solutions;
}

/**
 * Solve the linear diophantine equation ax + by = c for integers x,y
 *
 * There exists a solution if and only if d := gcd(a, b) divides c.
 * If there exists some solution, there are infinitely many solutions:
 * x_0 + b/d*t, y_0 - a/d*t, where (x_0, y_0) is a particular solution
 * and t is an integer.
 *
 * Use the Euclidean algorithm to find x,y such that ax + by = gcd(a, b) = d
 * Since d must divide c, we can multiply by c/d to get a(cx/d) + n(cy/d) = c
 * Therefore, (x_0, y_0) = (cx/d, cy/d) is a particular solution.
 *
 * @param a,b,c
 * @return some (x, y) pair such that ax + by = c, empty pair if no solution
 */
std::pair<long long int, long long int> linear_diophantine_solver(long long int a, long long int b, long long int c) {
    long long int x, y;
    long long int d = extended_gcd(a, b, x, y);
    if (c % d != 0) return {};
    
    // Now we have ax + by = d, so (x_0, y_0) = (cx/d, cy/d) is a particular solution.
    long long int x_0 = c / d * x; // c / d -> O( (lg floor(c/d))(lg d) ) + O( (lg(c/d))(lg x) )
    long long int y_0 = c / d * y;
    
    return {x_0, y_0};
}

/**
 * Solve a system of linear congruences of the following form for x, y:
 * ax + by = r (mod n)
 * cx + dy = s (mod n)
 *
 * There exists a unique solution modulo n if and only if gcd(ad-bc, n) = 1.
 * We want to solve (ad-bc)x = dr-bs (mod n) for x
 *        and solve (ad-bc)y = as-cr (mod n) for y
 * (Based off Thm 4.9 in Burton's Elementary Number Theory)
 *
 * @param a,b,r,c,d,s
 * @param n modulus
 * @return unique (x, y) pair that solves the system, empty pair if no solution
 */
std::pair<long long int, long long int> linear_congruence_system_solver(
    long long int a, long long int b, long long int r,
    long long int c, long long int d, long long int s,
    long long int n)
{
    long long int adbc = mod(a*d-b*c, n); // log(ad - bc) = O(log ad) = O(log(n^2)) = O(log n)
    if (gcd(adbc, n) != 1) return {};
    
    // Solve (ad-bc)x = dr-bs (mod n) for x
    auto x_solutions = linear_congruence_solver(adbc, mod(d*r-b*s, n), n);
    long long int x = x_solutions[0];
    
    // Solve (ad-bc)y = as-cr (mod n) for y
    auto y_solutions = linear_congruence_solver(adbc, mod(a*s-c*r, n), n);
    long long int y = y_solutions[0];
    
    return {x, y};
}

// TODO: Chinese Remainder Theorem Solver


/*
 * Examples running some of the code / comparing execution times
 * are given below in main
 */
int main()
{
    /*
     * Modular exponentiation
     */
    // long long int a = 5, b = 123456789, n = 12;
    long long int a = 5, b = 123456789, n = 12345;
    std::cout << "Let a = " << a << ", b = " << b << ", n = " << n << std::endl;
    std::cout << "Calculating a^e (mod n) = " << a << "^" << b << " mod " << n << std::endl;
    
    std::cout << "\n********** Calculating a^e (mod n) with Binary Exponentiation **********" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    long long int result = mod_pow(a, b, n);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by binary exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^b (mod n) with Binary exponentiation (recursive) **********"
              << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = mod_pow_recursive(a, b, n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by m-ary exponentiation: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^b (mod n) with m-ary exponentiation **********" << std::endl;
    long long int m = 128;
    std::cout << "Let m = " << m << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = mod_m_pow(a, b, n, m);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by m-ary exponentiation: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^b (mod n) with m-ary exponentiation (recursive) **********"
              << std::endl;
    m = 128;
    std::cout << "Let m = " << m << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = mod_m_pow_recursive(a, b, n, m);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by m-ary exponentiation (recursive): " << duration.count() << " nanoseconds"
              << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^e (mod n) with looping **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = slow_mod_pow(a, b, n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by loop exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
         
    std::cout << "\n********** Calculating a^e (mod n) with std::pow **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = (long long int) std::pow(a, b) % n;
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by library exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    /*
     * Modular multiplication
     */
    // int a = 50, b = 50000, n = 12;
    a = 50000; b = 1844674407370955101LL; n = 123456789;
    
    std::cout << "\nLet a = " << a << ", b = " << b << ", n = " << n << std::endl;
    
    std::cout << "\n********** Calculating a * b (mod n) with repeated squaring **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = mod_mult(a, b, n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by fast multiplication: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a * b (mod n) with looping **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = slow_mod_mult(a, b, n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by loop mult: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
         
    std::cout << "\n********** Calculating a * b (mod n) with standard '*' **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = (a * b) % n;
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by standard mult: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    /*
     * Solving linear congruences
     */
    std::cout << "\n********** Solving Linear Congruences **********" << std::endl;
    a = 3; b = 9; n = 12;
    std::cout << "Solving " << a << "x = " << b << " (mod " << n << ")" << std::endl;
    // Solve 3x = 9 (mod 12)
    auto solutions = linear_congruence_solver(a, b, n);
    std::cout << "All solutions (unique modulo " << n << "):" << std::endl;
    for (const auto& x : solutions) {
        std::cout << "x: " << x << std::endl;
    }
    
    /*
     * Solving linear diophantine equations
     */
    std::cout << "\n********** Solving Linear Diophantine Equations **********" << std::endl;
    a = 172; b = 20; n = 1000;
    std::cout << "Solving " << a << "x + " << b << "y = " << n << std::endl;
    // Solve 172x + 20y = 1000
    auto [x_0, y_0] = linear_diophantine_solver(172, 20, 1000);
    std::cout << "x_0: " << x_0 << " y_0: " << y_0 << std::endl;
    long long int d = gcd(a, b);
    std::cout << "\nAll solutions are given by\n" <<
        "x = " << x_0 << " + (" << b << "/" << d << ")*t = " << x_0 << " + " << b/d << "t\n" <<
        "y = " << y_0 << " - (" << a << "/" << d << ")*t = " << y_0 << " - " << a/d << "t"   << std::endl;
    
    /*
     * Solving systems of linear congruences of the form
     * ax + by = r (mod n)
     * cx + dy = s (mod n)
     */
    std::cout << "\n********** Solving Systems of Linear Congruences **********" << std::endl;
    a = 7; b = 3; d = 5; n = 16;
    long long int c = 2, r = 10, s = 9;
    std::cout << "Solving: \n" <<
        a << "x + " << b << "y = " << r << " (mod " << n << ")\n" <<
        c << "x + " << d << "y = " << s << " (mod " << n << ")"   << std::endl;
    // Solve 7x + 3y = 10 (mod 16)
    //       2x + 5y = 9 (mod 16)
    auto [x, y] = linear_congruence_system_solver(a, b, r, c, d, s, n);
    std::cout << "x = " << x << ", y = " << y << std::endl;
    
    return 0;
}
