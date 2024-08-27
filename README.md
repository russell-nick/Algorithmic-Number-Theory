# Algorithmic Number Theory

### Contents
1. [Project Overview](#project-overview)
2. [Future Work](#future-work)

## Project Overview
This project is a collection of various algorithms related to Algorithmic / Computational Number Theory. The main topics covered are general number theory, factorization, primality testing, discrete logarithms, psuedo-randomness, and other miscellaneous mathematical topics. 

This is mainly for personal / educational purposes, so not all algorithms have fully optimized implementations for practical use. For example, my implementation of Lenstra's Elliptic Curve Factoring method uses Boost's multiprecision library, for the arbitrary precision integer type `cpp_int`, over the GNU multiprecision library (which should be faster) for simplicity. 

*A brief list of currently implemented algorithms are listed below:*

### Current Algorithms
#### Number Theory
- Modular Arithmetic
	- Modulo operation, modular binary exponentiation, modular $m$-ary exponentiation, modular multiplication (binary exponentiation style addition), Jacobi symbol, Legendre symbol, determine if a number is quadratic residue, modular inverse, linear congruence solver, system of linear congruences solver, linear diophantine solver
- GCD / LCM
	- Euclidean algorithm, Extended Euclidean algorithm, LCM, GCD / LCM of a list of numbers
- Exponentiation
	- Binary exponentiation, $m$-ary exponentiation
- Logarithms
	- Base-b integer conversion, calculate $\lfloor\log_b(n)\rfloor$
- Binomial Coefficients
	- Given $Q \in \mathbb{Z}^+$, find all non-trivial $(n, k)$ pairs such that $Q = \binom{n}{k}$ in poly$($lg$(Q))$ time

#### Factorization
- Pollard's p-1 method
- Lenstra's Elliptic Curve method

#### Primality Testing
- Fermat primality testing
- Solovay-Strassen primality testing
- Miller-Rabin primality testing

## Future Work
Two major future plans are implementing more algorithms and potentially refactoring the project into a library. 

The project is currently set up as a collection of algorithms with reused code copied wherever it is needed. This makes it easy to compile and run certain algorithms without having to compile everything in the project together. However, this leads to duplicate code (for example, factoring algorithms use a lot of functions defined throughout the rest of the project) that would be difficult to maintain as the project grows. Therefore, I might later refactor this to be a single library instead of the current "collection of short disjoint programs".

*Below is a tentative list of algorithms and features to implement in the future, separated by topic:*

### To-do List:
#### Number Theory
- Modular Square Roots
	- Square roots modulo p (prime)
		- Tonelli-Shanks algorithm
		- Cipolla's algorithm
	- Square roots modulo p^k (prime power)
		- Hensel lifting
		- Modified Newton's method
	- Square roots modulo 2^k
		- Inductive solution to proving "If $a \equiv 1 \pmod{8}$, then $x^2 \equiv a \pmod{2^k}$ has a solution for x"
		- Taylor series with a quadratically convergent sequence (algorithm from homework)
- Chinese Remainder Theorem solver

#### Factoring
- Integer Factorization
	- Pollard's Rho algorithm
	- Dixon's algorithm (first provably subexponential factoring algorithm)
- Polynomial Factorization (over Finite Fields)
	- Berlekamp's algorithm

#### Discrete Logarithms
- Baby-step Giant-step algorithm
- Pollard's Kangaroo method

#### Pseudorandomness
- Predicting Lehmer random number generators

#### Miscellaneous
- Approximate $\ln(n)$ with absolute error $\le 1$ in poly$($lg$(n))$ time (using Maclaurin series)
- Test for nontrivial symmetry of curves with parametric form $f(t) = \displaystyle \sum_{n \in \mathbb{Z}} a_ne^{int}$, where the $a_n$ are real coefficients and $i = \sqrt{-1}$, and find the maximum value of $m$ such that $f(t)$ has $m$-fold symmetry
- Fast Fibonacci sequence calculations
- Fast multiplication (Karatsuba's algorithm, FFT multiplication, etc.)

#### Refactor project into a library
- Remove duplicate code and "helper" files. 
- Add useful namespaces.
- Add unit tests (using a unit test framework such as CppUnit or GoogleTest).
- Add build / test instructions.
