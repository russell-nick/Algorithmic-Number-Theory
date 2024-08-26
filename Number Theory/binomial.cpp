/**
 * @file binomial.hpp
 * @brief Implementation of algorithms related to binomial coefficients.
 *        This mainly focuses on finding all non-trivial (n, k) pairs such
 *        that Q = (n choose k) for a given Q (done in polynomial time).
 *
 * @author Nicholas Russell
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/number.hpp>

#include "binomial.hpp"

namespace mp = boost::multiprecision;

/**
 * Calculate the binomial coefficient n choose k:
 * n! / (k! (n-k)!) = [ n * (n-1) * ... * (n - k + 1)] / k!
 *
 * @param n binomial coefficient n
 * @param k binomial coefficient k
 * @return n choose k
 */
mp::cpp_int binomial_coefficient(mp::cpp_int n, mp::cpp_int k) {
    mp::cpp_int res = 1;
    for (mp::cpp_int i = n - k + 1; i <= n; ++i)
        res *= i;
        
    for (mp::cpp_int i = 2; i <= k; ++i)
        res /= i;
        
    return res;
}

/**
 * Check whether 'Q' is nontrivially a binomial coefficient in polynomial time.
 * That is, using time polynomial in log Q, is Q = (n choose k)
 * for some k, n with 1 < k < n − 1 ?
 *
 * Note that simply iterating Pascal's triangle (using only addition!)
 * or iterating values of 'n' to some bound while performing binary search on
 * 'k' values does not give a polynomial time algorithm.
 *
 * We need to use the following bounds (derivation/proof currently in private homework):
 * 2 < k < floor(2 * log_2(Q))
 * 0 < n < ceil((1 + sqrt(1 + 8*Q)) / 2) (must binary search this interval!)
 *
 * Bit Complexity: O((lg Q)^3 * (lg sqrt(Q))^3) bit operations
 *
 * @param Q number to check
 * @return true if Q is a nontrivial binomial coefficient, false otherwise
 */
bool is_binomial_coefficient(mp::cpp_int Q) {
    // K: Upper bound on values of k to check
    auto K1 = floor(2 * log2(mp::cpp_dec_float_50(Q)));
    mp::cpp_int K = K1.convert_to<mp::cpp_int>();
    std::cout << "K (upper bound on k) = " << K << std::endl;
    
    // N: Upper bound on values of n to check
    auto N1 = ceil( (1 + sqrt( mp::cpp_dec_float_50(1 + 8*Q))) / 2.0);
    mp::cpp_int N = N1.convert_to<mp::cpp_int>();
    std::cout << "N (upper bound on n) = " << N << std::endl;
    
    // L: Lower bound on values of n to check (log_2(Q))
    // This makes no difference on the complexity, but might
    // very slightly help in practice
    mp::cpp_int L = mp::msb(Q); // most significant bit = floor(log_2(Q))
    std::cout << "L (lower bound on n) = " << L << std::endl;
    
    /*
        To avoid using cpp_dec_float, one can use the slightly modified bound
        below (adds +1 to K and N at most). However, this will only save a
        fraction of a second at the cost of potentially doing one more iteration
        (which can be expensive for very large input)
    */
    // // K: Upper bound on values of k to check
    // mp::cpp_int K = 2 * mp::msb(Q) + 2;
    // std::cout << "K = " << K << std::endl;
    
    // // N: Upper bound on values of n to check
    // mp::cpp_int radicand = 1 + 8*Q;
    // mp::cpp_int N = ((1 + (sqrt(radicand) + 1 )) / 2) + 1;
    // std::cout << "N = " << N << std::endl;

    for (mp::cpp_int k = 2; k <= K; k++) {
        // To get polynomial time, we must also binary search on 0, ..., N:
        // or max(log_2(Q), k + 2), ..., N for tiny improvement in practice
        mp::cpp_int low = max(L, k+2), high = N;

        while (low <= high) {
            mp::cpp_int mid = low + (high - low) / 2;
            mp::cpp_int n_choose_k = binomial_coefficient(mid, k);
            
            if (n_choose_k == Q) {
                std::cout << "Q = " << mid << " choose " << k << std::endl;
                if (mid - k != k)
                    std::cout << "Q = " << mid << " choose " << mid - k << std::endl;
                return true;
            }
    
            // If mid choose k < Q, we want a larger n value
            if (n_choose_k < Q)
                low = mid + 1;
    
            // If mid choose k > Q, we want a smaller n value
            else
                high = mid - 1;
        }
        // If no value of n can satisfy Q = n choose k for the given k,
        // then continue to the next k value.
    }
    
    std::cout << "No non-trivial binomial coefficient found" << std::endl;
    return false;
}

/**
 * Find all (n, k) pairs such that Q = n choose k (non-trivially)
 *
 * Bit Complexity: O((lg Q)^3 * (lg sqrt(Q))^3) bit operations
 *
 * @param Q number to check
 * @return vector of all (n, k) pairs such that Q = n choose k (non-trivially)
 */
std::vector<std::pair<mp::cpp_int, mp::cpp_int>> all_binomial_coefficients(mp::cpp_int Q) {
    // Vector holding all (n, k) pairs such that Q = n choose k (non-trivially)
    std::vector<std::pair<mp::cpp_int, mp::cpp_int>> binom_coeffs;
    
    // K: Upper bound on values of k to check
    auto K1 = floor(2 * log2(mp::cpp_dec_float_50(Q)));
    mp::cpp_int K = K1.convert_to<mp::cpp_int>();
    std::cout << "K (upper bound on k) = " << K << std::endl;
    
    // N: Upper bound on values of n to check
    auto N1 = ceil( (1 + sqrt( mp::cpp_dec_float_50(1 + 8*Q))) / 2.0);
    mp::cpp_int N = N1.convert_to<mp::cpp_int>();
    std::cout << "N (upper bound on n) = " << N << std::endl;
    
    // L: Lower bound on values of n to check (log_2(Q))
    mp::cpp_int L = mp::msb(Q);
    std::cout << "L (lower bound on n) = " << L << std::endl;

    std::unordered_set<mp::cpp_int> used_k;
    for (mp::cpp_int k = 2; k <= K; k++) {
        // If we already found a (n, k) pair that uses this k, skip to next k
        if (used_k.count(k) > 0) continue;
        
        // To get polynomial time, we must also binary search on 0, ..., N:
        // or max(log_2(Q), k + 2), ..., N for tiny improvement in practice
        mp::cpp_int low = max(L, k+2), high = N;
        
        while (low <= high) {
            mp::cpp_int mid = low + (high - low) / 2;
            mp::cpp_int n_choose_k = binomial_coefficient(mid, k);
            
            if (n_choose_k == Q) {
                std::cout << "Q = " << mid << " choose " << k << std::endl;
                binom_coeffs.push_back({mid, k});
                used_k.insert(k);
                if (mid - k != k) {
                    std::cout << "Q = " << mid << " choose " << mid - k << std::endl;
                    binom_coeffs.push_back({mid, mid - k});
                    used_k.insert(mid - k);
                }
            }
    
            // If mid choose k < Q, we want a larger n value
            if (n_choose_k < Q)
                low = mid + 1;
    
            // If mid choose k > Q, we want a smaller n value
            else
                high = mid - 1;
        }
        // If no value of n can satisfy Q = n choose k for the given k,
        // then continue to the next k value.
    }
    
    if (binom_coeffs.size() == 0)
        std::cout << "No non-trivial binomial coefficient found" << std::endl;
        
    return binom_coeffs;
}


int main()
{
    // When testing, somewhere around 13,123,110 is the max number that works with long long int
    // (+ no boost ints to calculate binomial coefficient). When using cpp_int to calculate the
    // the binomial coefficients, this went up to around 137,846,528,820. Using boost cpp_int
    // for everything allows arbitrary sized integers.
    //
    // Some example #'s: 1184040, 13123110, 300540195, 137846528820, 118264581564861424
    mp::cpp_int Q("112186277816662845432"); // 21 digits (> 64 bit integer)
    
    // A (n, k) pair for following 150 digit number can be found in approximately 3 minutes.
    // (192818982242 nanoseconds when I first ran it: finds 500 choose 250)
    // mp::cpp_int Q("116744315788277682920934734762176619659230081180311446124100284957811112673608473715666417775521605376810865902709989580160037468226393900042796872256");
    
    std::cout << "********** Determining whether Q is a binomial coefficient **********" << std::endl;
    std::cout << "Q (input) = " << Q << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    is_binomial_coefficient(Q);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken: " << duration.count() << " nanoseconds" << std::endl;
    
    std::cout << "\n********** Finding all (n, k) pairs, where Q = n choose k **********" << std::endl;
    std::cout << "Q (input) = " << Q << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto binom_coeffs = all_binomial_coefficients(Q);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken to find all (n, k) pairs: " << duration.count() << " nanoseconds" << std::endl;
    
    for (auto& [n, k] : binom_coeffs) {
        std::cout << "(n, k): " << "(" << n << ", " << k << ")" << std::endl;
    }
    
    std::cout << "\n********** Finding all (n, k) pairs, where Q = n choose k (multiple pairs) **********"         << std::endl;
    Q = 3003; // 3003, 7140
    std::cout << "Q (input) = " << Q << std::endl;
    start = std::chrono::high_resolution_clock::now();
    binom_coeffs = all_binomial_coefficients(Q);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken to find all (n, k) pairs: " << duration.count() << " nanoseconds" << std::endl;
    
    for (auto& [n, k] : binom_coeffs) {
        std::cout << "(n, k): " << "(" << n << ", " << k << ")" << std::endl;
    }

    return 0;
}

