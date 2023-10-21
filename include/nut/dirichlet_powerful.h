#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Dirichlet Hyperbola extension functions revolving around expanding complicated multiplicative functions in terms of
/// simple ones plus functions restricted to powerful integers.
///
/// For a multiplicative function f, we can express f as g <*> h where g is equal to f at primes, h is zero except at "powerful"
/// numbers (those numbers n where p^2 divides n whenever p divides n), and of course f, g, and h are all multiplicative.
/// When we do this, if g is easy to compute/sum, we can compute/sum f just by iterating over all powerful numbers and evaluating h (and G)
/// at them.  Computing the adjustment h then becomes the hard part, and there are a lot of ways of doing this based on
/// the form of f.  Generally there are a few ways to compute h.  Often, f or at least just h when restricted to prime powers will depend
/// only on the exponent and not the prime.
/// Additionally, we can sometimes get a closed form for h using bell series, which are just a generating function of values of h at prime
/// powers.  In particular, the bell series of h_p(z) = sum(e >= 0, h(p^e)z^e).  In fact, bell series have the nice property that
/// (g <*> h)_p(z) = g_p(z)h_p(z).  However, it also often isn't possible or at least easy to find bell series for h.
///
/// Once again, griff's blog post is a good resource
/// (see here https://gbroxey.github.io/blog/2023/04/30/mult-sum-1.html)
/// (see here https://en.wikipedia.org/wiki/Powerful_number)
/// (see here https://en.wikipedia.org/wiki/Bell_series)
///
/// When we can't find the bell series, we can often compute h at prime powers recursively from f at prime powers, using the formula
/// h(p^e) = f(p^e) - sum(i = 0 ... e-1 : h(p^i)*g(p^(e-i)))
/// This basically lets us find h at all prime powers pretty easily whenever we can find f at all prime powers, which we always can if f is well-defined.
///
/// I'll briefly go over the problem that made me write this library as an example:
/// We want to count how many numbers up to max have exactly 8 divisors.
/// The number of divisors of a number n is d(n), which is a common multiplicative function.
/// Its values at prime powers of course are d(p^e) = e + 1, so we can see that numbers with exactly 8 divisors
/// have three possible forms:
/// n = p^7, n = p^3q, or n = pqr, where p, q, r are (distinct) primes.
/// So we want to find an indicator function for these numbers, that is, a function which is 1 at p^7, p^3q, and pqr, and 0 at other numbers.
/// Using zero divisors, this is easy, but of course, there are no zero divisors in the natural numbers, so we decide to extend the range
/// to be a ring with zero divisors, namely Z[x], the ring of polynomials with integer coefficients, and then we can work mod a power of x.
/// In particular, we will define f(1) = 1, f(p) = x, f(p^3) = x^2, f(p^7) = x^3, and f(p^r) = 0 for other prime powers.
/// Then, notice that if we work in Z[x]/(x^4), this f will be x^3 at n with exactly 8 divisors, and either x, x^2, or 0 otherwise.
/// Thus if we sum this f over all n up to max, we will get a polynomial
/// {#n with exactly 8 divisors}x^3 + {#n with exactly 4 divisors}x^2 + {#n with exactly 2 divisors}x + 1.
/// However, you may have noticed that we've only dicussed techniques for summing multiplicative functions whose range is Z, and while
/// they could be extended to larger "nice" rings like C or Z mod m, it's not easy to extend them to a fancy pants ring like Z[x]/(x^4).
/// So we won't actually find the sum of this f directly, but rather reconstruct it from several concrete values.
/// Consider how f(n) can be considered as a function of x, and even Z[x]/(x^4) can be considered as a function of x too:
/// If we pick an integer, like say 2^16, as a value for x, then polynomial evaluation is a homomorphism from Z[x]/(x^4) to Z mod x^4,
/// in other words, f(n)(2^16) mod 2^64 is a multiplicative function from N to Z mod 2^64.
/// You can understand this as the number b = 2^16 in Z mod 2^64 being a good but not quite perfect representation of x in Z[x]/(x^4):
/// both satisfy the equation x^4 = 0, but unfortunately 2^48 * 2^16 = 0 mod 2^64 while 2^48 * x is not 0 in Z[x]/(x^4).
/// We can use the chinese remainder theorem to combine the information we get from this for multiple b values, but it will still be hard
/// to separate the information about the coefficients on x^3, x^2, and x.
/// If we evaluate sum(n = 1, ..., max : f(n)(b) mod b^4) for enough values of b, then we can reconstruct sum(n = 1, ..., max : f(n)) and find
/// its x^3 coefficient.
/// We could find the coefficients of f using techniques like polynomial interpolation, but to me it makes more sense to rename f to f3
/// and define a second function f2
/// f2 : N -> Z[x]/(x^3)
/// f2(1) = 1, f2(p) = x, f2(p^3) = x^2, f2(p^r) = 0 otherwise
/// Thus the sum of f2(x) is {#n with exactly 4 divisors}x^2 + {#n with exactly 2 divisors}x + 1.
/// So if we call the sums of f3 and f2 F3 and F2 respectively, we can see that F3(max) - F2(max) = {#n with exactly 8 divisors}x^3.
/// Using concretizations f3(n)(b) mod b^4 and f2(n)(b) mod b^3 for multiple different b gives us enough information to find the x^3 coefficient
/// in F3(max) - F2(max), which is of course {#n with exactly 8 divisors}.
/// The only thing we have left to figure out is how to compute the sums of f3(n)(b) and f2(n)(b), but that's where the powerful number adjustment
/// comes in:
/// f3(n)(b) = (db <*> h3,b)(n)
/// Since f3(p)(b) = db(p), db matches f3 at primes and we use db as g in g <*> h, and we could already find h using the recursive formula given above.
/// Let's go over the Bell series anyway, since it is very instructive and actually useful sometimes.
/// (h3,b)_p(z) = (f3(.)(b))_p(z) / db_p(z)
/// From the identity that dk = u <*> d(k-1) and u_p(z) = sum(e >= 0 : z^e), db_p(z) = sum(e >= 0 : z^e)^b.
/// We can represent this as the sum of an infinite geometric series and instead write db_p(z) = 1/(1 - z)^b.
/// Then (f3(.)(b))_p(z) = (1 + bz + b^2z^3 + b^3z^7), so overall
/// (h3,b)_p(z) = (1 + bz + b^2z^3 + b^3z^7)(1 - z)^b


#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include <nut/modular_math.h>
#include <nut/dirichlet.h>

/// Compute the sum of f = g <*> h (from 1 up to g_tbl->x), where f(p) = g(p) and h(p^e) depends only on e
/// @param [out] out: store the result
/// @param [in] m: modulus to reduce result by, or 0 to not reduce
/// @param [in] g_tbl: table of values/sums for g, a simple function that's equal to f at primes
/// @param [in] h_vals: table of values of h at arbitrary prime exponents: h_vals[e] stores the value of h(p^e),
/// since this function requires h(p^e) depend only on e.  h_vals[0] must be 1, and h_vals must have entries up to
/// floor(log_2(g_tbl->x))
/// @return true on success, false on allocation failure
[[gnu::nonnull(1, 3, 4)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 3) NUT_ATTR_ACCESS(read_only, 4)
bool nut_Diri_sum_adjusted_he(int64_t *restrict out, int64_t m, const nut_Diri *restrict g_tbl, const int64_t *restrict h_vals);

/// Compute the sum of f = g <*> h (from 1 up to g_tbl->x), where f(p) = g(p)
/// Equivalently, h(pq) = 0 when p and q are coprime, so we can define h(p^e) (h at prime powers)
/// as a two parameter function.  If h only depends on e, use { @link nut_Diri_sum_adjusted_he}
/// @param [out] out: store the result
/// @param [in] m: modulus to reduce result by, or 0 to not reduce
/// @param [in] g_tbl: table of values/sums for g, a simple function that's equal to f at primes
/// @param [in] h_fn: callback to compute h at prime powers.  p is the prime (base), pp is the eth power of the prime,
/// e is the exponent, and m is forwarded
/// @return true on success, false on allocation failure
[[gnu::nonnull(1, 3, 4)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Diri_sum_adjusted_hpe(int64_t *restrict out, int64_t m, const nut_Diri *restrict g_tbl, int64_t (*h_fn)(uint64_t p, uint64_t pp, uint64_t e, int64_t m));

/// Compute the series quotient h = f / g for two (finite) power series
/// In other words, h will be such that f = g * h (under power series multiplication, not Dirichlet convolution)
/// Similar to polynomial division but cannot have a remainder
/// @param [in] n: length of the power series.  Since this is the same for all three, if one of the series (ie f) is shorter, it must be zero padded
/// @param [in] m: modulus to reduce result by, or zero to not reduce
/// @param [out] h: store the quotient
/// @param [in] f: dividend
/// @param [in] g: divisor (MUST have g[0] = 1)
[[gnu::nonnull(3, 4, 5)]]
NUT_ATTR_ACCESS(read_write, 3) NUT_ATTR_ACCESS(read_only, 4) NUT_ATTR_ACCESS(read_only, 5)
void nut_series_div(uint64_t n, int64_t m, int64_t h[restrict static n], int64_t f[restrict static n], int64_t g[restrict static n]);

