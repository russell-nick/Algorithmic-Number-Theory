/**
 * @file exponentiation.cpp
 * @brief Implementation of various exponentiation algorithms.
 *
 * @author Nicholas Russell
 */

#include <bitset>
#include <chrono>
#include <cmath>
#include <iostream>

#include "exponentiation.hpp"


/**
 * Compute a^e using binary exponentiation (repeated squaring)
 *
 * Let b_k*2^k + ... + b_1*2 + b_0 be the binary representation of
 * the exponent e. We use the fact that a^(b + c) = a^b * a^c
 * to calculate a^e using O(log(e)) multiplications:
 *
 * a^e = a^(b_k*2^k + ... + b_1*2 + b_0) = a^(b_k * 2^k)*...*a^(b_1 * 2)*a^(b_0)
 * = a^{(b_k)(2^k)}*...*a^{(b_1)(2)}*a^(b_0)
 * so we can repeatedly square 'a' k-1 times and
 * multiply the result by 'a' if the current bit b_i = 1
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
long long int bin_pow(long long int a, long long int e) {
    long long int result = 1;
    while (e > 0) {
        // if current bit in binary rep is = 1
        if (e & 1) { // or e % 2 == 0
            result *= a;
        }
        a *= a;
        // shift to next bit in exponent:
        e >>= 1; // or e /= 2
    }
    return result;
}

/*
 * Compute a^e using binary exponentiation (repeated squaring)
 * by looping through the bitset of the exponent.
 *
 * If using fixed sized integers, one can get the bitset of the exponent
 * and do the same thing as above (bitset holds coefficients b_0, ..., b_k)
 */
// long long int bin_pow(int a, int e) {
//     std::bitset<32> exponent(e);
//     // cout << exponent << endl;
//     // cout << exponent.size() << endl;
//     long long int result = 1;
//     for (int i = 0; i <= exponent.size(); i++) {
//         if (exponent[i] == 1) {
//             result = result * a;
//         }
//         a *= a;
//     }
//     return result;
// }

/**
 * Compute a^e using binary exponentiation (repeated squaring) recursively
 * using the formula:
 * a^e = a * a^(e-1) if e is odd,
 * a^e = (a^(e/2))^2 if e is even
 *
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
long long int bin_pow_recursive(long long int a, long long int e) {
    if (e == 0) return 1;
    if (e % 2 == 0) { // e even
        long long int t = bin_pow_recursive(a, e / 2);
        return t * t;
    } else { // e odd
        long long int t = bin_pow_recursive(a, e - 1);
        return a * t;
    }
    return 0;
}

/**
 * Compute a^e using m-ary exponentiation
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
 * Complexity: O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param m base (radix) of exponent to work with
 * @return a^e
 */
long long int m_pow(long long int a, long long int e, long long int m) {
    long long int result = 1;
    while (e > 0) {
        long long int rem = e % m;
        if (rem != 0) {
            result = result * bin_pow(a, rem);
        }
        a = bin_pow(a, m);
        e /= m;
    }
    return result;
}

/**
 * Compute a^e using m-ary exponentiation (recursively) using the formula
 * a^e = a^(e mod m) * ( a^((e - (e mod m)) / m) )^m if e != 0 (mod m)
 * a^e = (a^(e/m))^m if e = 0 (mod m)
 *
 * (If dealing with very large numbers (> long long), one can also precompute
 *  and use a^0, a^1, a^2, ..., a^(m-1) to avoid the use of using something
 *  like binary exponentation for the a^(e mod m) calculation)
 *
 * It is best to let m be a power of 2.
 *
 * Complexity: <= m + log(e) + log_m(e) = O(lg e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @param m base (radix) of exponent to work with
 * @return a^e
 */
long long int m_pow_recursive(long long int a, long long int e, long long int m) {
    if (e == 0) return 1;
    if (e % m == 0) {
        long long int t = m_pow_recursive(a, e / m, m);
        return bin_pow(t, m);
    } else {
        long long int rem = e % m;
        long long int t = m_pow_recursive(a, (e - rem) / m, m);
        return bin_pow(a, rem) * bin_pow(t, m);
    }
    return 0;
}

/* Note:
 *
 * One can also do binary exponentiation for floating point bases.
 * As long as the exponent is an integer, we can get O(log e) multiplications
 * instead of O(e) multiplications (or instead of using some floating point formula
 * to estimate a^e, such as std::pow)
 */

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
double bin_pow(double a, long long int e) {
    // If e < 0, we are calculating (1/a)^e
    if (e < 0) a = 1 / a;
    
    double result = 1;
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
 * Compute a^e by repeatedly multiplying by 'a' e times
 *
 * Complexity: O(e) multiplications
 *
 * @param a number to raise to exponent
 * @param e exponent >= 0
 * @return a^e
 */
long long int slow_pow(long long int a, long long int e) {
    long long int result = 1;
    for (int i = 0; i < e; i++) {
        result *= a;
    }
    return result;
}

int main()
{
    long long int num = 2;
    long long int exponent = 62;
    std::cout << "Calculating a^e = " << num << "^" << exponent << std::endl;
    
//    int64_t n = 2;
//    int64_t e = 63;
//    int64_t r = 1;
//    for (int i = 0; i < e; i++) {
//        r *= n;
//    }
//    std::cout << "int 64 2^63: " << r << std::endl;

    std::cout << "\n********** Calculating a^e with Binary Exponentiation **********" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    long long int result = bin_pow(num, exponent);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by binary exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^e with Binary Exponentiation (recursive) **********"
              << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = bin_pow_recursive(num, exponent);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by binary exp (recursive): " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^e with m-ary exponentiation **********" << std::endl;
    long long int m = 8;
    std::cout << "Let m = " << m << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = m_pow(num, exponent, m);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by m-ary exp (recursive): " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^e with m-ary exponentiation (recursive) **********"
              << std::endl;
    std::cout << "Let m = " << m << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = m_pow_recursive(num, exponent, m);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by m-ary exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
    
    std::cout << "\n********** Calculating a^e with looping **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = slow_pow(num, exponent);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by loop exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
         
    std::cout << "\n********** Calculating a^e with std::pow **********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    result = std::pow(num, exponent);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken by library exp: " << duration.count() << " nanoseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;

    return 0;
}
