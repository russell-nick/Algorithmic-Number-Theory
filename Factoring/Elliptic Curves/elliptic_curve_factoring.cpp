/**
 * @file elliptic_curve_factoring.cpp
 * @brief Implementation of Lenstra's Elliptic Curve Factorization Method.
 *
 * Note: For improved factoring, one should use a combination of different techniques
 *       such as trial division up to a small bound, pollard rho, quadratic sieve / ecm, etc.
 *       This implementation currently focuses on the Elliptic Curve Method combined with basic trial divison.
 *
 * @author Nicholas Russell
 */

#include <chrono>
#include <cmath>    // exp, sqrt, log
#include <iostream>
#include <random>
#include <string>
#include <utility>  // pair
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "elliptic_curve_factoring.hpp"
#include "elliptic_curve.hpp"
#include "helpers.hpp"

namespace mp = boost::multiprecision;


/**
 * Try to compute a multiple of the order of a group related to the factoring algorithm:
 * Here we want to consider |E(Z/pZ)| for some prime p dividing n
 * If E = k|E(Z/pZ)|, then E*P = (0:1:0) and we can extract a non-trivial divisor
 *
 * We use E = prod_{prime p <= B} p^(floor(log B / log p))
 *          = prod_{prime p <= B} p^(floor(log_p(B)) = lcm(1...B)
 *
 * One can also use E = prod_{prime p <= B} p^(floor(log n / log p))
 *                    = prod_{prime p <= B} p^(floor(log_p(n))
 * for better odds of factoring, although E will be much larger
 * so factoring could take longer.
 *
 * Notes / Alternate Methods for computing E when n is very large:
 * (1) Implement segmented sieve and compute
 *      E = prod_{prime p <= B} p^(floor(log B / log p))
 *        = prod_{prime p <= B} p^(floor(log_p(B)) = lcm(1...B) directly as before
 * (2) Instead of a sieve, perform primality tests on 1, ..., B (simple, but could be wrong)
 * (3) Compute lcm(1...B) = E = prod_{prime p <= B} p^(floor(log B / log p))
 *     prime-by-prime by starting with 2, 3, then iterating to the next primes as follows:
 *     Keep incrementing by 2 until prime (use miller-rabin primality test)
 *     (this improves on (2) by not having to store the primes <= B and skips primality
 *      testing on every other number)
 * (4) Use (3) and instead of computing E entirely beforehand and then trying to calculate E*P,
 *     Calculate E*P prime by prime:
 *     If 'primes' contains all primes <= B, then
 *     for prime in primes:
 *         Find largest e such that p^e <= B (e = floor(log_p(B))
 *         Compute P := p^e * P
 *         Return d if partial multiplication finds nontrivial factor d
 *     (at the end, we would have computed E*P)
 *
 * (Methods (1) and (2) will require a ton of memory to store all ~ n/log(n) (Prime # Thm)
 * prime numbers <= n when n is large, so they will not be practical.)
 *
 * @param B smoothness bound
 * @return lcm(1...B) = prod_{prime p <= B} p^(floor(log B / log p))
 *
 */
mp::cpp_int compute_E(const mp::cpp_int& B) {
    mp::cpp_int lcm = 1;
    for (mp::cpp_int i = 2; i <= B; i++) {
        lcm = lcm / gcd(lcm, i) * i;
    }
    return lcm;
}
/*
 * Directly compute lcm(1...B) = prod_{prime p <= B} p^(floor(log B / log p))
 */
//mp::cpp_int compute_E(mp::cpp_int B) {
//    mp::cpp_int E = 1;
//
//    std::vector<mp::cpp_int> primes = sieve_of_eratosthenes(B);
//    for (auto p : primes) {
//        long long int exponent = floor_log(B, p); // floor(log_p(n))
//        E *= pow(p, exponent); // replace pow with binary exp / loop
//    }
//    return E;
//}

/**
 * "Psuedo-addition" from Lenstra's Elliptic Curve Factorization paper.
 * This is essentially applying mutliplication (repeated squaring) until the
 * group operation fails and stopping early if we find a non-trivial divisor.
 *
 * Based on the elliptic curve addition law (here, O = (0:1:0) denotes identity):
 * P + O = O + P = P for any point P
 *
 * Let P = (x1:y1:1) and Q = (x2:y2:1) be points on the elliptic curve y^2 = x^3 + ax + b.
 * If Q = -P (-P = (x1:-y1:1)), then P + Q = O. Otherwise, calculate slope lambda and add
 * as follows with P + Q = R = (x3:y3:1):
 *
 * lambda = (y2 - y1)/(x2 - x1) if P != Q (or (y1 - y2)/(x1 - x2) will also end up working)
 * lambda = (3*x1^2 + a)/(2*y1) if P = Q
 *
 * x3 = lambda^2 - x1 - x2
 * y3 = lambda(x1 - x3) - y1
 *
 * Complexity: O((lg n)^2) bit operations (dominated by gcd, multiplication, mod)
 *
 * @param P first point in addition
 * @param Q second point in addition
 * @return (R, d) pair, where R = P+Q (if no nontrivial divisor is found),
 *                      d is the nontrivial divisor if one is found (if R = inf),
 *                      or d = 0 if no nontrivial divisor is found (R != inf)
 *
 */
std::pair<ECPoint, mp::cpp_int> partial_addition(const ECPoint& P, const ECPoint& Q) {
    if (P.inf) return {Q, 0};
    if (Q.inf) return {P, 0};

    mp::cpp_int n = P.E->n;
    mp::cpp_int x, y;
    // Using mod(P.x - Q.x, n) to get the least non-negative residue is required
    // here since we need the modular inverse of (P.x - Q.x)
    mp::cpp_int d = extended_gcd(mod(P.x - Q.x, n), n, x, y);
    
    // Now we have x, y such that (P.x - Q.x)x + ny = d
    // or (P.x - Q.x)x = d (mod n)

    if (d > 1 && d < n) return {P, d}; // non-invertible denominator (+ fails): return d

    if (d == 1) // P.x != Q.x
    {
        // Since d = 1, we have (P.x - Q.x)x = 1 (mod n)
        // so x is the multiplicative inverse of (P.x - Q.x)
        // (denominator in this addition law case)
        mp::cpp_int denom = x % n;                           // Optional: Replace % n with mod function
        mp::cpp_int lambda = ((P.y - Q.y) * denom) % n;      // Optional: Replace % n with mod function
        
        // Get x, y coordinates for P + Q:
        mp::cpp_int x_3 = (lambda * lambda - P.x - Q.x) % n; // Optional: Replace % n with mod function
        mp::cpp_int y_3 = (lambda * (P.x - x_3) - P.y) % n;  // Optional: Replace % n with mod function
        return {ECPoint(x_3, y_3, false, P.E), 0};
    }
    else // d = n => P.x - Q.x = 0 (mod n) => P.x = Q.x
    {
        // Using mod(P.y + Q.y, n) to get the least non-negative residue is required
        // here since we need the modular inverse of (P.y + Q.y)
        d = extended_gcd(mod(P.y + Q.y, n), n, x, y);
        if (d > 1 && d < n) return {P, d}; // non-invertible denominator (+ fails): return d
        
        // If P.y = -Q.y, the denominator for this addition law (P.y + Q.y = 0 mod n) and we
        // get the point at infinity:
        if (d == n) {
            std::cout << "group operation fails: denominator = 0, gcd(P.y + Q.y, n) = n" << std::endl;
            return {ECPoint(P.x, P.y, true, P.E), 0};
        }

        // Since d = 1, we have (P.y + Q.y)x = 1 (mod n)
        // so x is the multiplicative inverse of (P.y + Q.y)
        // (denominator in this addition law case)
        mp::cpp_int denom = x % n;                           // Optional: Replace % n with mod function
        mp::cpp_int lambda = ((3 * P.x * P.x + P.E->a) * denom) % n;
        
        // Get x, y coordinates for P + Q:
        mp::cpp_int x_3 = (lambda * lambda - P.x - Q.x) % n; // Optional: Replace % n with mod function
        mp::cpp_int y_3 = (lambda * (P.x - x_3) - P.y) % n;  // Optional: Replace % n with mod function
        return {ECPoint(x_3, y_3, false, P.E), 0};
    }
}

/**
 * "Psuedo-multiplication" from Lenstra's Elliptic Curve Factorization paper
 * (terminate early if we get a nontrivial divisor of n)
 *
 * This is just binary exponentiation (repeated squaring) in the "group"
 * E(Z/nZ) until we extract a non-trivial divisor from a failed group operation
 * (see "elliptic_curve.h" for an explanation of this "group" which need not be
 *  well-defined for this algorithm to work)
 *
 * Complexity: O(lg k) group operations
 *             O((lg k)(lg n)^2) bit operations
 *
 * @param k scalar
 * @param P point on elliptic curve
 * @return (R, d) pair, where R = k*P (if no nontrivial divisor is found: d = 1 or d = n)
 *                      or d is a nontrivial divisor (R != k*P since we exit early)
 *
 */
std::pair<ECPoint, mp::cpp_int> partial_multiplication(mp::cpp_int k, ECPoint P) {
    // Start with identity
    ECPoint result = ECPoint(P.x, P.y, true, P.E);
    mp::cpp_int d = 0, n = P.E->n;
    while (k > 0) {
        if (k % 2 == 1) { // k & 1
            std::tie(result, d) = partial_addition(result, P);
            // If we found a non-trivial divisor when trying to add result + P:
            if (d > 1 && d < n) return {result, d};
        }
        std::tie(P, d) = partial_addition(P, P);
        // If we found a non-trivial divisor when trying to add P + P:
        if (d > 1 && d < n) return {result, d};
        k /= 2; // k >>= 1;
    }
    return {result, 0};
}

/**
 * Lenstra's Elliptic Curve Factoring (to split n)
 *
 * We use the smoothness bound B = L[1/2, 1/sqrt(2)] = exp( sqrt(1/2 * ln(n) * ln(ln(n))) )
 * as mentioned in Lenstra's paper. This is similar to the choice of B for factoring
 * algorithms such as Dixon's algorithm, the Quadratic Sieve, Pollard's p-1 method, etc.
 *
 * The complete algorithm calculates B,
 * takes E = lcm(1...B) = prod_{prime p <= B} p^(floor(log B / log p)),
 * chooses a random elliptic curve and point P, tries to compute E*P,
 * and finds a non-trivial divisor if the group operation
 * (Elliptic curve addition over Z/nZ) fails (denominator in addition law will be
 *  non-invertible modulo n).
 *
 * Complexity:
 * Using repeated squaring to calculate E*P uses
 * O(log E) = O(log (prod_{prime p <= B} p^(floor(log B / log p)) )
 *          = O(sum_{prime p <= B} (log(B)/log(p))*log(p) )
 *          = O(sum_{prime p <= B} log(B) )
 *          = O(pi(B) * log(B)) = O( (B / log(B)) * log(B)) = O(B) by the Prime # Thm
 * group operations for a total of O(B (lg n)^2) bit operations to split n.
 *
 * The (conjectured!) expected bit operations to find a factor of n is O(B (lg n)^2) * # Curves to check
 * which works out to be O( exp(sqrt( (2 + o(1)) * ln p ln ln p )) * (lg n)^2) bit operations,
 * where p is the smallest prime divisor of n.
 *
 * @param n number to factor
 * @param B smoothness bound
 * @return a non-trivial divisor of n
 *
 */
mp::cpp_int lenstra_elliptic_curve(const mp::cpp_int& n, const mp::cpp_int& B, const mp::cpp_int& E) {

    /*
     * Generate a random elliptic curve with Weierstrass form y^2 = x^3 + ax + b
     * and a random point on the curve:
     * Uniformly at random choose a, x, y in Zn
     * Set b = y^2 - x^3 + ax
     * If the curve is singular (4a^3 + 27b^2 = 0), try again
     */
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<mp::cpp_int> distribution(0, n-1);
    mp::cpp_int a = distribution(generator), x = distribution(generator), y = distribution(generator);
    mp::cpp_int b = (y * y - (x * x * x) - a*x) % n;

    while ( (4*a*a*a + 27*b*b) % n == 0) {
        a = distribution(generator);
        x = distribution(generator);
        y = distribution(generator);
        b = (y * y - (x * x * x) - a*x) % n;
    }

    std::cout << "\n********** Generating Elliptic Curve y^2 = x^3 + ax + b **********\n"
              << "Factoring n = " << n << "\n"
              << "a: " << a << " b: " << b << "\n"
              << "P = (x:y:1) = (" << x << ":" << y << ":1) " << std::endl;
    
    // Try gcd(4a^3 + 27b^2, n) to see if we are lucky (gcd != 1). Otherwise, continue.
    mp::cpp_int d1 = gcd(mod(4*a*a*a + 27*b*b, n), n);
    if (d1 < 1 && d1 < n) {
        std::cout << "Lucky split before EC factoring!: 1 < gcd(discriminant, n) < n" << std::endl;
        return d1;
    }
    
    /*
     * Here, B = L[1/2, 1/sqrt(2)] = exp( sqrt(1/2 * ln(n) * ln(ln(n))) )
     * and we slowly increase the bound as we do more unsuccessful iterations
     *
     * We also take E = prod_{prime p <= B} p^(floor(log B / log p)) = lcm(1...B)
     *
     * Both are passed as arguments to avoid recomputation if the current
     * elliptic curve fails to find a non-trivial factor.
     */
    std::cout << "B: " << B << std::endl;
    //mp::cpp_int E = compute_E(B);
    //std::cout << "E: " << E << std::endl;
    
    EllipticCurve* C = new EllipticCurve(a, b, n);
    ECPoint P = ECPoint(x, y, false, C);

    /*
     * Alternatively, one could do this prime by prime
     * instead of calculating E = lcm(1...B) and attempting E*P:
     * If 'primes' contains all primes <= B, then
     * for prime in primes:
     *      Find largest e such that p^e <= B (e = floor(log_p(B))
     *      Compute P := p^e * P
     *      Return d if partial multiplication finds nontrivial factor d
     *
     * (at the end, we would have computed E*P)
     * (for large n, see (3) in the comment for compute_E to do this probabilistically
     *  without storing all primes <= B)
     */
    auto [EP, d] = partial_multiplication(E, P);
    std::cout << "Partial E*P: " << EP << " d: " << d << std::endl;

    delete C;
    return d;
}

/**
 * Merge the prime decompositions of two non-trivial divisors of n
 *
 * Since we have to merge the factors of the nontrivial divisors of n anyway, we
 * might as well perform a merge sort at the same cost to sort the prime factors.
 * (realistically, this will be negligible since n will most likely
 *  not have many distinct prime factors)
 *
 * @param f1,f2 prime decompositions of two non-trivial divisors of n
 * @return prime decomposition of f1*f2 as vector of pairs (p_1, e_1), ..., (p_r, e_r),
 *         where f1*f2 = p_1^{e_1} * ... * p_r^{e_r} is the prime decomposition of f1*f2
 */
std::vector<std::pair<mp::cpp_int, mp::cpp_int>> merge_factors(
    const std::vector<std::pair<mp::cpp_int, mp::cpp_int>>& f1,
    const std::vector<std::pair<mp::cpp_int, mp::cpp_int>>& f2
) {
    std::vector<std::pair<mp::cpp_int, mp::cpp_int>> factors;
    long long int i = 0, j = 0, n1 = f1.size(), n2 = f2.size();
    if (f1.size() == 0) return f2;
    if (f2.size() == 0) return f1;
    
    factors.reserve(n1 + n2);
    while (i < n1 && j < n2) {
        if (f1[i].first < f2[j].first) {
            factors.push_back(f1[i]);
            i++;
        } else if (f2[j].first < f1[i].first) {
            factors.push_back(f2[j]);
            j++;
        } else {
            factors.push_back({f1[i].first, f1[i].second + f2[j].second});
            i++; j++;
        }
    }
    
    while (i < n1) {
        factors.push_back(f1[i]);
        i++;
    }
    while (j < n2) {
        factors.push_back(f2[j]);
        j++;
    }
    return factors;
}

/**
 * Full factorization using Lenstra's Elliptic Curve method.
 *
 * Before using elliptic curves, we use the following quick checks:
 * (1) Test if the number is prime
 * (2) Perform trial division with small primes
 *
 * If an elliptic curve does not produce a non-trivial divisor, try another one.
 * Once a sufficient number of curves have been tested, increase the smoothness bound.
 * Repeat until a non-trivial divisor d is found and recurse on d and n/d
 *
 * @param n number to factor
 * @return vector of pairs (p_1, e_1), ..., (p_r, e_r),
 *         where n = p_1^{e_1} * ... * p_r^{e_r} is the prime decomposition of n
 */
std::vector<std::pair<mp::cpp_int, mp::cpp_int>> factor(mp::cpp_int n) {
    // If prime with probability >= 1 - (1/4)^k:
    int k = 8;
    if (miller_rabin(n, k)) {
        return {{n, 1}};
    }
    
    // TODO: Implement perfect power testing here:
    
    // Trial division by primes <= 100,000
    // (one can also just replace the for loop with trial division by 2 <= i <= 100,000)
    std::string small_factor_decomposition = "";
    std::vector<std::pair<mp::cpp_int, mp::cpp_int>> small_factors;
    auto small_primes = sieve_of_eratosthenes(100000);
    for (long long int p : small_primes) {
        long long int exp = 0;
        while (n % p == 0) {
            exp++;
            n /= p;
        }
        if (exp > 0) {
            small_factor_decomposition += std::to_string(p) + "^" + std::to_string(exp) + " * ";
            small_factors.push_back({p, exp});
            if (exp == 1) {
                std::cout << "Non-trivial prime divisor " << p << " was found!" << std::endl;
            } else {
                std::cout << "Non-trivial prime power divisor " << p << "^" << exp << " was found!" << std::endl;
            }
        }
    }
    
    // Check if we fully factored n with trial division:
    if (small_factors.size() > 0 && n == 1) {
        // Remove "* " from end of the prime decomposition output:
        small_factor_decomposition.pop_back(); small_factor_decomposition.pop_back();
        std::cout << "n was quickly factored by trial division!" << std::endl;
        std::cout << "n = " << small_factor_decomposition << std::endl;
        return small_factors;
    }
    
    // If n was partially factored with trial division:
    if (small_factors.size() > 0) {
        std::cout << "Initial decomposition after trial division with small (<= 100,000) primes:\nn = "
                  << small_factor_decomposition << n << std::endl;
        auto f = factor(n);
        return merge_factors(small_factors, f);
    }

    /*
     * Compute the smoothness bound here (instead of inside lenstra's elliptic curve factoring)
     * to avoid recomputations when the current curve fails to find a non-trivial divisor:
     *
     * Choose B = L[1/2, 1/sqrt(2)] = exp( sqrt(1/2 * ln(n) * ln(ln(n))) )
     * and slowly increase the bound as we do more iterations
     */
    auto B1 = exp(sqrt(0.5 * log(mp::cpp_dec_float_50(n)) * log(log(mp::cpp_dec_float_50(n)))));
    mp::cpp_int B = B1.convert_to<mp::cpp_int>();
    
    // If n is not prime, find a non-trivial divisor:
    mp::cpp_int E = compute_E(B);
    mp::cpp_int divisor = lenstra_elliptic_curve(n, B, E);
    //mp::cpp_int divisor = lenstra_elliptic_curve(n, B);
    int num_curves = 1;
    while (divisor == 0) {
        std::cout << "Failed to find non-trivial divisor. Generating new elliptic curve." << std::endl;
        divisor = lenstra_elliptic_curve(n, B, E);
        //divisor = lenstra_elliptic_curve(n, B);
        num_curves++;
        if (num_curves > 20) {
            B *= 2;
            std::cout << "Increasing B: " << B << std::endl;
            num_curves = 0;
        }
    }
    mp::cpp_int divisor2 = n / divisor;
    std::cout << "Split n: " << n << " = " << divisor << " * " << divisor2 << std::endl;
    
    // Recurse to find prime factors of the nontrivial divisors 'divisor' and 'divisor2'
    auto factors1 = factor(divisor);
    auto factors2 = factor(divisor2);
    
    // Merge the prime factors of divisor and divisor2:
    return merge_factors(factors1, factors2);
}

/**
 * Prints the prime decomposition of a factored number
 *
 * @param n number to factor
 * @param factors vector of pairs (p_1, e_1), ..., (p_r, e_r),
 *        where n = p_1^{e_1} * ... * p_r^{e_r} is the prime decomposition of n
 */
void print_factors(const mp::cpp_int& n, const std::vector<std::pair<mp::cpp_int, mp::cpp_int>>& factors) {
    std::cout << "\n********** Factored n **********" << std::endl;
    std::cout << n << " = ";
    for (int i = 0; i < factors.size() - 1; i++) {
        std::cout << factors[i].first << "^" << factors[i].second << " * ";
    }
    std::cout << factors.back().first << "^" << factors.back().second << std::endl;
}

int main() {
    
    //mp::cpp_int n("398883434337287");
    //mp::cpp_int n("46167045131415113");
    //mp::cpp_int n("64211816600515193");
    //mp::cpp_int n("168541512131094651323");
    //mp::cpp_int n("631211032315670776841");
    //mp::cpp_int n("4132846513818654136451");
    //mp::cpp_int n("4516511326451341281684513");
    //mp::cpp_int n("3146531246531241245132451321");
    //mp::cpp_int n("4269021180054189416198169786894227");
    //mp::cpp_int n("7853316850129");
    
    //mp::cpp_int n("18644474572985789975837398427254139850324935126641");
    //mp::cpp_int n("610703718544355446139790085966"); // 2^1 * 542237^1 * 645187^1 * 49688141^1 * 17566007477^1
    //mp::cpp_int n("1847294298523"); // 1007557 * 1833439 * 6788563
    
    /*
     * Factor with Lenstra's Elliptic Curve Method
     */
    std::cout << "\n********** Factoring with Elliptic Curve **********" << std::endl;
    mp::cpp_int n("41498298392926620982104401176");
    std::cout << n << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto factors = factor(n);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken: " << duration.count() << " nanoseconds" << std::endl;
    print_factors(n, factors);

    return 0;
}
