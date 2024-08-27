/**
 * @file gcd_lcm.cpp
 * @brief Implementation of various GCD and LCM algorithms.
 *
 * @author Nicholas Russell
 */

#include <iostream>
#include <vector>

#include "gcd_lcm.hpp"

/**
 * Calculate the gcd of two integers: 'a' and 'b'.
 *
 * Complexity: Lame (~1844) showed that the Euclidean algorithm
 * uses <= 5log_10(min(a, b)) steps, so the runtime is O(log(min(a, b))).
 *
 * Bit Complexity: # Steps * # Bit operations for division
 * = O((log(min(a, b))) * (log(max(a, b)))^2) bit operations
 
 * However, using Collins' bound, one gets the tighter bound
 * of O((lg a)(lg b)) bit operations.
 *
 * @param a,b
 * @return gcd(a, b)
 */
long long int gcd(long long int a, long long int b) {
    long long int temp;
    while (b > 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/**
 * Recursively calculate the gcd of two integers: 'a' and 'b'.
 *
 * Complexity: O(log(min(a, b))) steps
 * Bit Complexity: O((lg a)(lg b)) bit operations
 *
 * @param a,b
 * @return gcd(a, b)
 */
long long int gcd_recursive(long long int a, long long int b) {
    if (b == 0) return a;
    return gcd_recursive(b, a % b);
}

/**
 * Use the Extended Euclidean Algorithm to find gcd(a, b) and solve
 * ax + by = gcd(a, b) for x, y
 *
 * Complexity (same as gcd): O(log(min(a, b))) steps
 * Bit Complexity (same as gcd): O((lg a)(lg b)) bit operations
 *
 * @param[in] a,b
 * @param[out] x, y integers x,y such that ax + by = gcd(a, b)
 * @return gcd(a, b)
 */
long long int extended_gcd(long long int a, long long int b, long long int& x, long long int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    
    long long int x1, y1;
    long long int d = extended_gcd(b, a % b, x1, y1);
    
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

/**
 * Compute lcm(a, b) = ab / gcd(a, b)
 *
 * Complexity (same as gcd): O(log(min(a, b))) divisions
 * Bit Complexity (same as gcd): O((lg a)(lg b)) bit operations.
 *
 * @param a,b
 * @return lcm(a, b)
 */
long long int lcm(long long int a, long long int b) {
    return a / gcd(a, b) * b;
}

/**
 * Compute gcd(a_0, a_1, ..., a_n), where
 * nums = {a_0, a_1, ..., a_n}
 *
 * @param nums vector of numbers to take gcd of
 * @return gcd of all integers in 'nums'
 */
long long int gcd(const std::vector<long long int>& nums) {
    if (nums.size() < 2) return 1;
    
    long long int result = gcd(nums[0], nums[1]);
    
    for (int i = 2; i < nums.size(); i++) {
        if (result == 1) return result;
        result = gcd(result, nums[i]);
    }
    
    return result;
}

/**
 * Compute lcm(a_0, a_1, ..., a_n), where
 * nums = {a_0, a_1, ..., a_n}
 *
 * @param nums vector of numbers to take lcm of
 * @return lcm of all integers in 'nums'
 */
long long int lcm(const std::vector<long long int>& nums) {
    if (nums.size() == 0) return 1;
    
    long long int lcm = nums[0];
    for (int i = 1; i < nums.size(); i++) {
        lcm = lcm / gcd(lcm, nums[i]) * nums[i];
    }
    return lcm;
}

//int main() {
//    std::vector<long long int> v = {6, 9, 3, 12};
//    std::cout << gcd(v) << std::endl;
//
//    long long int x, y;
//    long long int d = extended_gcd(6, 9, x, y);
//    std::cout << "d: " << d << std::endl;
//    std::cout << "x: " << x << " y: " << y << std::endl;
//
//    return 0;
//}
