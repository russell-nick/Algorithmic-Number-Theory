/**
 * @file rational.cpp
 * @brief Implementation for Rational number class.
 *
 * @author Nicholas Russell
 */

#include "rational.hpp"
#include <cmath>

long long int gcd(long long int a, long long int b) {
    a = abs(a); b = abs(b);
    long long int temp;
    while (b > 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

Rational::Rational() : num(0), denom(1) {}
Rational::Rational(long long int n) : num(n), denom(1) {}
Rational::Rational(long long int n, long long int d) : num(n), denom(d) {}

long long int Rational::numerator() const {return num;}
long long int Rational::denominator() const {return denom;}
long long int Rational::to_int() const {return static_cast<long long int>(num / denom);};

Rational& Rational::operator+= (long long int x) {
    num += (x * denom);
    return *this;
}
Rational& Rational::operator-= (long long int x) {
    num -= (x * denom);
    return *this;
}
Rational& Rational::operator*= (long long int x) {
    long long int d = gcd(denom, x);
    num *= (x / d);
    denom /= d;
    return *this;
}
Rational& Rational::operator/= (long long int x) {
    long long int d = gcd(num, x);
    num /= d;
    denom *= (x / d);
    return *this;
}

//Rational Rational::operator+ (long long int x) const {
//    Rational temp = Rational(*this);
//    temp.num += (x * temp.denom);
//    return temp;
//}
//Rational Rational::operator- (long long int x) const {
//    Rational temp = Rational(*this);
//    temp.num += (x - temp.denom);
//    return temp;
//}
//Rational Rational::operator* (long long int x) const {
//    Rational temp = Rational(*this);
//    long long int d = gcd(temp.denom, x);
//    temp.num *= (x / d);
//    temp.denom /= d;
//    return temp;
//}
//Rational Rational::operator/ (long long int x) const {
//    Rational temp = Rational(*this);
//    long long int d = gcd(temp.num, x);
//    temp.num /= d;
//    temp.denom *= (x / d);
//    return temp;
//}
Rational operator+ (const Rational& r, long long int x) {
    Rational temp = Rational(r);
    temp.num += (x * temp.denom);
    return temp;
}
Rational operator- (const Rational& r, long long int x) {
    Rational temp = Rational(r);
    temp.num += (x - temp.denom);
    return temp;
}
Rational operator* (const Rational& r, long long int x) {
    Rational temp = Rational(r);
    long long int d = gcd(temp.denom, x);
    temp.num *= (x / d);
    temp.denom /= d;
    return temp;
}
Rational operator/ (const Rational& r, long long int x) {
    Rational temp = Rational(r);
    long long int d = gcd(temp.num, x);
    temp.num /= d;
    temp.denom *= (x / d);
    return temp;
}

Rational operator+ (long long int x, const Rational& r) {return r + x;}
Rational operator- (long long int x, const Rational& r) {return r - x;}
Rational operator* (long long int x, const Rational& r) {return r * x;}
Rational operator/ (long long int x, const Rational& r) {return r / x;}

Rational& Rational::operator+= (const Rational& other) {
    // (num / denom) + (other.num / other.denom) = (other.denom*num + denom*other.num) / (denom * other.denom)
    // Naive implementation:
    num = other.denom * num + denom * other.num;
    denom = denom * other.denom;
    
    long long int d = gcd(num, denom);
    num /= d;
    denom /= d;
    
    return *this;
}
Rational& Rational::operator-= (const Rational& other) {
    // (num / denom) + (other.num / other.denom) = (other.denom*num - denom*other.num) / (denom * other.denom)
    // Naive implementation:
    num = other.denom * num - denom * other.num;
    denom = denom * other.denom;
    
    long long int d = gcd(num, denom);
    num /= d;
    denom /= d;
    
    return *this;
}
Rational& Rational::operator*= (const Rational& other) {
//    num *= other.num;
//    denom *= other.denom;
    
    long long int d1 = gcd(num, other.denom);
    long long int d2 = gcd(denom, other.num);

    // (num/denom) * (other.num/other.denom) = (num * other.num) / (denom * other.denom)
    num = (num / d1) * (other.num / d2);
    denom = (denom / d2) * (other.denom / d1);

    return *this;
}
Rational& Rational::operator/= (const Rational& other) {
//    num *= other.denom;
//    denom *= other.num;
    
    long long int d1 = gcd(num, other.num);
    long long int d2 = gcd(denom, other.denom);

    // (num/denom) / (other.num/other.denom) = (num * other.denom) / (denom * other.num)
    num = (num / d1) * (other.denom / d2);
    denom = (denom / d2) * (other.num / d1);

    return *this;
}

Rational Rational::operator+ (const Rational& other) const {
    Rational temp = Rational(*this);
    
    // (num / denom) + (other.num / other.denom) = (other.denom*num + denom*other.num) / (denom * other.denom)
    // Naive implementation:
    temp.num = other.denom * temp.num + temp.denom * other.num;
    temp.denom = temp.denom * other.denom;
    
    long long int d = gcd(temp.num, temp.denom);
    temp.num /= d;
    temp.denom /= d;
    
    return temp;
}
Rational Rational::operator- (const Rational& other) const {
    Rational temp = Rational(*this);
    
    // (num / denom) + (other.num / other.denom) = (other.denom*num - denom*other.num) / (denom * other.denom)
    // Naive implementation:
    temp.num = other.denom * temp.num - temp.denom * other.num;
    temp.denom = temp.denom * other.denom;
    
    long long int d = gcd(temp.num, temp.denom);
    temp.num /= d;
    temp.denom /= d;
    
    return temp;
}
Rational Rational::operator* (const Rational& other) const {
    Rational temp = Rational(*this);
    
    long long int d1 = gcd(temp.num, other.denom);
    long long int d2 = gcd(temp.denom, other.num);

    // (num/denom) * (other.num/other.denom) = (num * other.num) / (denom * other.denom)
    temp.num = (temp.num / d1) * (other.num / d2);
    temp.denom = (temp.denom / d2) * (other.denom / d1);
    
    return temp;
}
Rational Rational::operator/ (const Rational& other) const {
    Rational temp = Rational(*this);
    
    long long int d1 = gcd(temp.num, other.num);
    long long int d2 = gcd(temp.denom, other.denom);

    // (num/denom) / (other.num/other.denom) = (num * other.denom) / (denom * other.num)
    temp.num = (temp.num / d1) * (other.denom / d2);
    temp.denom = (temp.denom / d2) * (other.num / d1);
    
    return temp;
}

const Rational& Rational::operator++() {
    num += denom;
    return *this;
}
const Rational& Rational::operator--() {
    num -= denom;
    return *this;
}
Rational Rational::operator++(int) {
    Rational temp = Rational(*this);
    ++(*this);
    return temp;
}
Rational Rational::operator--(int) {
    Rational temp = Rational(*this);
    --(*this);
    return temp;
}

std::ostream& operator<<(std::ostream& os, const Rational& obj) {
    return os << obj.num << "/" << obj.denom;
}
