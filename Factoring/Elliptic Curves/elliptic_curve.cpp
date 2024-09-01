/**
 * @file elliptic_curve.cpp
 * @brief Implementation of Elliptic Curve and Elliptic Curve Point classes.
 *
 * @author Nicholas Russell
 */

#include <ostream>
#include <stdexcept>

#include <boost/multiprecision/gmp.hpp>

#include "elliptic_curve.hpp"

namespace mp = boost::multiprecision;

EllipticCurve::EllipticCurve(mp::mpz_int a, mp::mpz_int b, mp::mpz_int n) : a(a), b(b), n(n) {
    if (discriminant() == 0)
        throw std::invalid_argument("Attempted to initialize an elliptic curve with singular curve");
};

mp::mpz_int EllipticCurve::discriminant() {
    return (4*(this->a * this->a * this->a) + 27*(this->b * this->b)) % this->n;
}

std::ostream& operator<<(std::ostream& os, const EllipticCurve& obj) {
    return os << "y^2 = x^3 + " << obj.a << "x + " << obj.b << " (mod " << obj.n << ")";
}


ECPoint::ECPoint(mp::mpz_int x, mp::mpz_int y, bool inf, EllipticCurve* E) : x(x), y(y), inf(inf), E(E) {};

// TODO: Add mod_inverse function from modular.h and verify that this elliptic curve addition below works.
/*
 * Elliptic curve addition law (here, O = (0:1:0) denotes identity):
 * P + O = O + P = P for any point P
 *
 * Let P = (x1:y1:1) and Q = (x2:y2:1) be points on the elliptic curve y^2 = x^3 + ax + b.
 * If Q = -P (-P = (x1:-y1:1)), then P + Q = O. Otherwise, calculate slope lambda and add
 * as follows with P + Q = R = (x3:y3:1):
 *
 * lambda = (y2 - y1)/(x2 - x1) if P != Q
 * lambda = (3*x1^2 + a)/(2*y1) if P = Q
 *
 * x3 = lambda^2 - x1 - x2
 * y3 = lambda(x1 - x3) - y1
 *
 */
//ECPoint ECPoint::operator+(const ECPoint& P) {
//    // P + O = O + P = P for any point P
//    if (this->inf) return P;
//    if (P.inf) return *this;
//
//    mp::mpz_int n = this->E->n;
//    mp::mpz_int lambda;
//
//    // Not the same point
//    if (this->x != P.x || this->y != P.y) {
//        // (this also convers when x1=x2, y1=-y2: P + -P = O
//        if (this->x == P.x) lambda = 0; // return inf
//        lambda = ((this->y - P.y) * mod_inverse(this->x - P.x, n)) % n;
//
//    } else { // Doubling the same point
//        if (this->y == 0) lambda = 0; // return inf
//        // this will fail if Z/nZ is not a field: ex: Z/4Z: 2 != 0 but 2*2 = 4 = 0 has no inverse
//        lambda = ((3 * this->x * this->x + this->E->a) * mod_inverse(2 * this->y, n)) % n;
//    }
//
//    if (lambda == 0) return ECPoint{this->x, this->y, true, this->E}; // inf
//    // x3 = lambda^2 - x1 - x2
//    // y3 = lambda(x1 - x3) - y1
//    mp::mpz_int x, y;
//    x = (lambda * lambda - this->x - P.x) % n;
//    y = (lambda * (this->x - x) - this->y) % n;
//    return ECPoint{x, y, false, this->E};
//}

std::ostream& operator<<(std::ostream& os, const ECPoint& obj) {
    if (obj.inf) return os << "(0:1:0)";
    return os << "(" << obj.x << ":" << obj.y << ":1)";
}
