/**
 * @file rational.hpp
 * @brief Header file for Rational number class.
 *
 * @author Nicholas Russell
 */

#ifndef RATIONAL_HPP
#define RATIONAL_HPP

#include <ostream>

/**
 * @brief Class to support rational number arithmetic.
 *
 * Note: The Rational class is used to get exact arithmetic when
 * dealing with rational numbers, as opposed to using floating point
 * numbers and accumulating error. Currently, this is used to support
 * calculating the n-th fibonacci number using Binet's formula in the
 * extension field Q(sqrt(5)).
 *
 * @see fibonacci.cpp for use case example.
 *
 */
class Rational {
private:
    long long int num;
    long long int denom;
public:
    Rational();
    Rational(long long int n);
    Rational(long long int n, long long int d);
    
    long long int numerator() const;
    long long int denominator() const;
    long long int to_int() const;
    
    Rational& operator+= (long long int x);
    Rational& operator-= (long long int x);
    Rational& operator*= (long long int x);
    Rational& operator/= (long long int x);
    
    friend Rational operator+ (const Rational& r, long long int x);
    friend Rational operator- (const Rational& r, long long int x);
    friend Rational operator* (const Rational& r, long long int x);
    friend Rational operator/ (const Rational& r, long long int x);
    friend Rational operator+ (long long int x, const Rational& r);
    friend Rational operator- (long long int x, const Rational& r);
    friend Rational operator* (long long int x, const Rational& r);
    friend Rational operator/ (long long int x, const Rational& r);
//    Rational operator+ (long long int x) const;
//    Rational operator- (long long int x) const;
//    Rational operator* (long long int x) const;
//    Rational operator/ (long long int x) const;
    
    Rational& operator+= (const Rational& other);
    Rational& operator-= (const Rational& other);
    Rational& operator*= (const Rational& other);
    Rational& operator/= (const Rational& other);
    
    Rational operator+ (const Rational& other) const;
    Rational operator- (const Rational& other) const;
    Rational operator* (const Rational& other) const;
    Rational operator/ (const Rational& other) const;
    
    // Prefix
    const Rational& operator++(); // or Rational& operator++(); ?
    const Rational& operator--(); // or Rational& operator--(); ?
    
    // Postfix
    Rational operator++(int);
    Rational operator--(int);
    
    friend std::ostream& operator<<(std::ostream& os, const Rational& obj);
};

#endif /* RATIONAL_HPP */
