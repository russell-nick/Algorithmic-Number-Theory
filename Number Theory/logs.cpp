/**
 * @file logs.cpp
 * @brief Implementation of algorithms related to logarithms.
 *
 * @author Nicholas Russell
 */

#include <iostream>
#include <vector>

#include "logs.hpp"

/**
 * Compute floor(log_b(n)) (also index of most significant coefficient (bit when b = 2))
 *
 * To calculate log_b(n), we want to find x such that b^x = n. However, to find floor(logb n), we can
 * continuously divide n by b as long as n >= b. If we divide a total of m times, then 1 <= n/(b^m) < b
 * so log_b(b^m) <= log_b(n) < log_b(b^(m+1)) and m <= log_b(n) < m + 1. Thus, floor(log_b(n)) = m.
 *
 * Alternatively, consider the base-b representation n = a_m*b^m + ... + a_1*b + a_0, where
 * a_i in {0, 1, ..., m} and a_m > 0. Then b^(m+1) > m*b^m + ... + m*b + m >= n >= b^m so
 * taking log_b of b^(m+1) > n >= b^m gives m+1 > n >= m, so floor(log_b(n)) = m.
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param n number
 * @param b base
 * @return floor(log_b(n))
 */
long long int floor_log(long long int n, long long int b) {
    long long int result = 0;
    while (n >= b) {
        n /= b;
        result++;
    }
    return result;
}

/**
 * Convert n to base b: (change b to just int b?)
 * n = a_k*b^k + ... + a_1*b + a_0
 *
 * Complexity: O((lg n)^2) bit operations
 *
 * @param n number
 * @param b base
 * @return coefficients a_k, ..., a_1, a_0 in reverse order
 */
std::vector<long long int> base_b_conversion(long long int n, long long int b) {
    std::vector<long long int> base_b_rep;
    while (n > 0) {
        long long int remainder = n % b;
        base_b_rep.push_back(remainder);
        n /= b;
    }
    return base_b_rep;
}

/*
 * Estimate log_b(n) (time polynomial in log n) with absolute error <= 1
 *
 * Complexity: O() bit operations
 *
 * @param n
 * @param b
 * @return log_b(n)
 */
// double log_b();

int main() {
    long long int n = 123456789;
    long long int b = 2;
    
    auto base_b_rep = base_b_conversion(n, b);
    
    std::cout << "********** Converting n to base-b **********" << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "b = " << b << std::endl;
    for (long long int i = base_b_rep.size() - 1; i >= 0; i--)
        std::cout << base_b_rep[i] << " ";
    
    return 0;
}
