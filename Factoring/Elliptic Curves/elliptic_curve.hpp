/**
 * @file elliptic_curve.hpp
 * @brief Header file for Elliptic Curve and Elliptic Curve Point classes.
 *
 * @author Nicholas Russell
 */

#ifndef ELLIPTIC_CURVE_HPP
#define ELLIPTIC_CURVE_HPP

#include <ostream>
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

/**
 * @brief Simple class for elliptic curves with Weierstrass form
 * (field characteristic != 2, 3) y^2 = x^3 + ax + b (mod n)
 *
 * Note: Typically we take elliptic curves over fields, but one can also form
 * elliptic curves over (commutative) rings. Z/nZ is a field if and only if n is prime,
 * so E(Z/nZ) can only be an elliptic curve (over the field Z/nZ) if n is prime.
 *
 * For the sake of Lenstra's elliptic curve factoring method, we don't enforce this,
 * as we can use "psuedo-addition" and apply the group operation until it fails
 * (at which point we get a nontrivial divisor of n and stop early). However, this
 * "pseudo-addition" method will give normal elliptic curve addition if Z/nZ is a field.
 *
 * According to Lenstra's paper, the addition laws are more general for the ring Z/nZ,
 * and we can form an elliptic curve over the commutative ring Z/nZ as long as
 * 6(4a^3 + 27b^2) in (Z/nZ)^*, or gcd(6, n) = 1 and the curve is nonsingular.
 *
 */
class EllipticCurve
{
public:
    mp::cpp_int a, b, n;
    
    EllipticCurve(mp::cpp_int a, mp::cpp_int b, mp::cpp_int n);
    mp::cpp_int discriminant();
    friend std::ostream& operator<<(std::ostream& os, const EllipticCurve& obj);
};

/**
 * @brief Class to hold a point (x:y:1) on an elliptic curve E.
 *
 * Note: As above, E does not actually have to form an elliptic curve. However,
 * the typical elliptic curve addition (over a field) defined here will work
 * if and only if n is prime.
 *
 */
class ECPoint
{
public:
    EllipticCurve* E;   // Elliptic Curve the point is on
    mp::cpp_int x, y;   // Point's x,y coordinates
    bool inf;           // true if point at infinity (identity (0:1:0))
    
    ECPoint(mp::cpp_int x,  mp::cpp_int y, bool inf, EllipticCurve* E);
    //ECPoint operator+(const ECPoint& P); // defined iff n is prime
    friend std::ostream& operator<<(std::ostream& os, const ECPoint& obj);
};

#endif /* ELLIPTIC_CURVE_HPP */
