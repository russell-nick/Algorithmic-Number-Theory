/**
 * @file fibonacci.cpp
 * @brief Implementation for various algorithms to compute fibonacci numbers.
 *
 * @author Nicholas Russell
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "fibonacci.hpp"
#include "rational.hpp"

/**
 * Compute a^e using binary exponentiation (repeated squaring)
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
double bin_pow(double a, long long int e) {
    double result = 1.0;
    while (e > 0) {
        if (e & 1) {
            result *= a;
        }
        a *= a;
        e >>= 1;
    }
    return result;
}

/**
 * Compute a * b in the extension field Q(sqrt(5))
 *
 * // TODO: Add description
 */
std::pair<double, double> mult_ext_field(const std::pair<double, double> a, const std::pair<double, double> b) {
    // ac + 5bd, ad + bc
    return {a.first * b.first + 5 * a.second * b.second, a.first * b.second + a.second * b.first};
}
std::pair<Rational, Rational> mult_ext_field_rational(const std::pair<Rational, Rational>& a,
                                                      const std::pair<Rational, Rational>& b) {
    // ac + 5bd, ad + bc
    return {a.first * b.first + a.second * b.second * 5, a.first * b.second + a.second * b.first};
}

/**
 * Compute a^e using binary exponentiation in the extension field Q(sqrt(5))
 *
 * // TODO: Add description
 */
std::pair<double, double> bin_pow_ext_field(std::pair<double, double> a, long long int e) {
    // Initialize result to identity element:
    std::pair<double, double>  result = {1, 0};
    while (e > 0) {
        if (e & 1) {
            result = mult_ext_field(result, a);
        }
        a = mult_ext_field(a, a);
        e >>= 1;
    }
    return result;
}
std::pair<Rational, Rational> bin_pow_ext_field_rational(std::pair<Rational, Rational> a, long long int e) {
    // Initialize result to identity element:
    std::pair<Rational, Rational> result = {Rational(1, 1), Rational(0, 1)};
    while (e > 0) {
        if (e & 1) {
            result = mult_ext_field_rational(result, a);
        }
        a = mult_ext_field_rational(a, a);
        e >>= 1;
    }
    return result;
}

/**
 * Compute mat1 * mat2 for 2x2 matrices mat1 and mat2
 *
 * // TODO: Add description
 */
std::vector<std::vector<long long int>> mat2d_mult(const std::vector<std::vector<long long int>>& mat1,
                                                   const std::vector<std::vector<long long int>>& mat2) {
    long long int a = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[1][0];
    long long int b = mat1[0][0] * mat2[0][1] + mat1[0][1] * mat2[1][1];
    long long int c = mat1[1][0] * mat2[0][0] + mat1[1][1] * mat2[1][0];
    long long int d = mat1[1][0] * mat2[0][1] + mat1[1][1] * mat2[1][1];
    return {{a, b}, {c, d}};
}

/**
 * Compute mat^e for a 2x2 matrix 'mat'
 *
 * // TODO: Add description
 */
std::vector<std::vector<long long int>> mat2d_bin_pow(std::vector<std::vector<long long int>>& mat, long long int e) {
    // Initialize result to identity element:
    std::vector<std::vector<long long int>> result = {{1, 0}, {0, 1}};
    while (e > 0) {
        if (e & 1) {
            result = mat2d_mult(result, mat);
        }
        mat = mat2d_mult(mat, mat);
        e >>= 1;
    }
    return result;
}

/**
 * Compute the n-th fibonacci number naively using its recursive formula:
 * f(n) = f(n-1) + f(n-2)
 * f(0) = 0
 * f(1) = 1
 *
 * Complexity: O(2^n) operations (additions)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_rec(long long int n) {
    if (n == 0) return 0;
    if (n == 1) return 1;
    return fibonacci_rec(n - 1) + fibonacci_rec(n - 2);
}

/**
 * Compute the n-th fibonacci number with dynamic programming
 * using its recursive formula:
 * f(n) = f(n-1) + f(n-2)
 * f(0) = 0
 * f(1) = 1
 *
 * Complexity: O(n) operations (additions)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci(long long int n) {
    std::vector<long long int> fib(n + 1);
    fib[0] = 0;
    fib[1] = 1;
    
    for (long long int i = 2; i <= n; i++) {
        fib[i] = fib[i - 1] + fib[i - 2];
    }
    
    return fib[n];
}

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
long long int fibonacci_const_space(long long int n) {
    if (n == 0) return 0;
    if (n == 1) return 1;
    
    long long int fib_prev1 = 1;
    long long int fib_prev2 = 0;
    
    long long int fib = 0;
    for (long long int i = 2; i <= n; i++) {
        fib = fib_prev1 + fib_prev2;
        fib_prev2 = fib_prev1;
        fib_prev1 = fib;
        
    }
    return fib;
}

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
 * Complexity: O(log n) operations (multiplications)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_mat_exp(long long int n) {
    std::vector<std::vector<long long int>> fib_mat = {{1, 1}, {1, 0}};
    fib_mat = mat2d_bin_pow(fib_mat, n);
    return fib_mat[0][1];
}

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
long long int fibonacci_binet(long long int n) {
    double sqrt5 = sqrt(5);
    double phi = (1 + sqrt5) / 2;
    double psi = (1 - sqrt5) / 2;
    double fib = (bin_pow(phi, n) - bin_pow(psi, n)) / sqrt5;
    return static_cast<long long int>(fib);
    
//    // Alternatively:
//    double sqrt5 = sqrt(5);
//    double phi = (1 + sqrt5) / 2;
//    double fib = bin_pow(phi, n) / sqrt5;
//    return static_cast<long long int>(fib);
}

/**
 * Compute the n-th fibonacci number with Binet's formula (exact arithmetic):
 * f(n) = (phi^n - psi^n) / sqrt(5), where
 * phi = (1 + sqrt(5)) / 2
 * psi = (1 - sqrt(5)) / 2
 *
 * Note: To get exact arithmetic, we do computations in the extension field Q(sqrt(5))
 * to avoid the irrational number sqrt(5) in Binet's formula.
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
long long int fibonacci_binet_exact(long long int n) {
    std::pair<double, double> phi = {0.5, 0.5};
    std::pair<double, double> psi = {0.5, -0.5};
    
    std::pair<double, double> phi_npow = bin_pow_ext_field(phi, n);
    std::pair<double, double> psi_npow = bin_pow_ext_field(psi, n);
    
    double fib = phi_npow.second - psi_npow.second;
    return static_cast<long long int>(fib);
}

/**
 * Compute the n-th fibonacci number with Binet's formula (exact arithmetic):
 * f(n) = (phi^n - psi^n) / sqrt(5), where
 * phi = (1 + sqrt(5)) / 2
 * psi = (1 - sqrt(5)) / 2
 *
 * Note: To get exact arithmetic, we do computations in the extension field Q(sqrt(5))
 * to avoid the irrational number sqrt(5) in Binet's formula. This implementation also
 * uses a custom rational number class to avoid all floating point arithmetic.
 *
 * Complexity: O(log n) operations (multiplications for exponents including n)
 *
 * @param n
 * @return n-th fibonacci number
 */
long long int fibonacci_binet_exact_rational(long long int n) {
    std::pair<Rational, Rational> phi = {Rational(1, 2), Rational(1, 2)};
    std::pair<Rational, Rational> psi = {Rational(1, 2), Rational(-1, 2)};

    std::pair<Rational, Rational> phi_npow = bin_pow_ext_field_rational(phi, n);
    std::pair<Rational, Rational> psi_npow = bin_pow_ext_field_rational(psi, n);

    Rational fib = phi_npow.second - psi_npow.second;
    return fib.to_int();
}

int main() {
    long long int n = 80;
    
    std::cout << "Let n = " << n << std::endl;
    std::cout << "Calculating the n-th fibonacci number Fn" <<  std::endl;
    
//    std::cout << "\n********** Calculating Fn with recursion **********" << std::endl;
//    auto start = std::chrono::high_resolution_clock::now();
//    long long int result = fibonacci_rec(n);
//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
//    std::cout << "Time taken by recursion: " << duration.count() << " nanoseconds" << std::endl;
//    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating Fn with dynamic programming **********" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    long long int result = fibonacci(n);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by dynamic programming: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating Fn with dynamic programming (constant space) **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = fibonacci_const_space(n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by dynamic programming (constant space): " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating Fn with matrix binary exponentation **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = fibonacci_mat_exp(n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by matrix binary exponentation: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating Fn with Binet's formula **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = fibonacci_binet(n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by Binet's formula: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating Fn with Binet's formula (exact arithmetic) **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = fibonacci_binet_exact(n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by Binet's formula (exact arithmetic): " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating Fn with Binet's formula (exact arithmetic + rational class) **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = fibonacci_binet_exact_rational(n);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by Binet's formula (exact arithmetic + rational class): " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
//    for (long long int i = 0; i < 100; i++) {
//        std::cout << "N = " << i << std::endl;
//        // std::cout << "rec: " << fibonacci_rec(i) << std::endl;
//        std::cout << "fib: " << fibonacci(i) << std::endl;
//        std::cout << "fib const: " << fibonacci_const_space(i) << std::endl;
//        std::cout << "fib mat: " << fibonacci_mat_exp(i) << std::endl;
//        std::cout << "fib binet: " << fibonacci_binet(i) << std::endl;
//        std::cout << "fib binet exact: " << fibonacci_binet_exact(i) << std::endl;
//        std::cout << "\n";
//    }
    
    std::cout << "\nTesting Rational Class" << std::endl;
    Rational x = Rational();
    Rational y = Rational(10);
    Rational z = Rational(5, 2);
    
    std::cout << "x: " << x << std::endl
              << "y: " << y << std::endl
              << "z: " << z << std::endl;
    
    std::cout << "x num: " << x.numerator() << " x denom: " << x.denominator() << std::endl
              << "y num: " << y.numerator() << " y denom: " << y.denominator() << std::endl
              << "z num: " << z.numerator() << " z denom: " << z.denominator() << std::endl;
    
    std::cout << "++x: " << ++x << std::endl;
    std::cout << "x++: " << x++ << std::endl;
    std::cout << "x: " << x << std::endl;
    
    y *= z;
    std::cout << "y *= z: " << y << std::endl;
    y /= z;
    std::cout << "y /= z: " << y << std::endl;
    
    std::cout << "y + 5: " << y + 5 << std::endl;
    std::cout << "y * 5: " << y * 5 << std::endl;
    // std::cout << "5 * y: " << 5 * y << std::endl;
    std::cout << "y * 5 / 2: " << (y * 5) / 2 << std::endl;
    
    y += Rational(5, 2);
    std::cout << "y += (5/2) : " << y << std::endl;
    y -= Rational(5, 2);
    std::cout << "y -= (5/2) : " << y << std::endl;
    std::cout << "y + (5/2) : " << y + Rational(5, 2) << std::endl;
    std::cout << "y - (5/2) : " << y - Rational(5, 2) << std::endl;
    
    Rational a = Rational(-1, 2);
    std::cout << "a: " << a << std::endl;
    std::cout << "a + 1/2: " << a + Rational(1, 2) << std::endl;
    std::cout << "a * -1/2: " << a * Rational(-1, 2) << std::endl;
    
    return 0;
}
