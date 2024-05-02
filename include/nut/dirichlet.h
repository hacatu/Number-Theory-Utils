#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Dirichlet Hyperbola based functions for computing sums of multiplicative functions at a single value quickly
///
/// For a function f defined over the natural numbers, we say f is multiplicative if f(ab) = f(a)f(b) whenever a and b are coprime
/// (have gcd 1).  There are many functions in number theory that satisfy this condition, with one of the most well known being
/// euler's phi function.
///
/// A brief list of common multiplicative functions is:
/// (see here https://en.wikipedia.org/wiki/Multiplicative_function)
/// The dirichlet convolution identity (explained later): I(n) = 1 if n == 1; 0 otherwise
/// Sometimes written epsilon(n), (very confusingly) u(n), or delta(n, 0) (since it is just a kroneker delta function of course)
/// The constant function: u(n) = 1
/// Sometimes written 1(n)
/// The identity function: N(n) = n
/// Sometimes written Id(n) or id(n)
/// The power functions: N_k(n) = n**k
/// Sometimes written Id_k(n)
/// The divisor count function: d(n) = #{k | n} (the number of natural numbers including 1 and n which divide n)
/// The generalized divisor functions: d_k(n) = # k-tuples ks of natural numbers such that ks multiplies out to n
/// The divisor power sum functions: sigma_a(n) = sum(k | n, k**a) (note that a can be any real number)
/// The mobius function: mu(n) = 0 if n is not squarefree; (-1)**Omega(n) otherwise, where Omega(n) is the number of prime factors of n with multiplicity
/// Euler's totient function: phi(n) = #{k = 1 ... n : (k coprime n)}
/// Any constant raised to the power of Omega(n) or omega(n) (the number of prime factors of n counted with and without multiplicity respectively)
/// The number of non-isomorphic abelian groups of order n: a(n)
/// The Ramanujan tau function: tau(n)
/// The unitary divisor power sums: sigma_a^*(n) = sum(d | n where (d coprime n/d), d**a)
/// All dirichlet characters, including the legendre symbol (n/p) for fixed p and gcd(n, m) for fixed m
///
/// Some of these are fully multiplicative, meaning not only do they satisfy f(ab) = f(a)f(b) for a, b coprime, but they satisfy it for all a, b.
/// In particular, I, u, N, and N_k are fully multiplicative.
///
/// We often want to find the sum of a multiplicative function, often denoted by a capital version, for instance phi and Phi, f and F (for a general function, etc).
/// We also often want to be able to analyze multiplicative functions in terms of simpler functions.
///
/// Multiplicative functions are fully determined by their values at prime powers, which often leads to efficient sieve based approaches to calculating them.
/// This header provides some functions for doing this, the nut_euler_sieve family of functions, although these are mainly intended as low level
/// subroutines for the upcoming Dirichlet Hyperbola algorithm.
///
/// However, if we only want to compute the sum of a multiplicative function up to some value n, sieve based methods will basically be O(nlogn),
/// which seems very bad considering how structured multiplicative functions are.
///
/// Indeed, we can often do better, in particular when we want to find the sum H(n) of some multiplicative function h(n) defined as h = f <*> g
/// where <*> denotes Dirichlet convolution.  What is Dirichlet convolution?  It is an operation for combining functions, similar to ordinary
/// multiplication, division, addition, subtraction, composition, and so on.  It is defined as (f <*> g)(n) = sum(k | n, f(k)g(n/k)).
/// This is mostly useful because when f and g are multiplicative f <*> g will be as well, and it lets us build up more complicated functions in terms
/// of simpler ones.
/// (see here https://en.wikipedia.org/wiki/Dirichlet_convolution)
///
/// Obviously, not all operations on multiplicative functions will produce multiplicative functions.  For example, if f is multiplicative, 2f is not,
/// so multiplicative functions are not closed under scalar multiplication, and thus they are not closed under addition.
/// They are closed under multiplication and dirichlet convolution though.
///
/// The Dirichlet Hyperbola algorithm lets us rewrite H(x) = sum(n <= x, h(n)) in terms of F and G as well as f and g, but crucially it lets us evaluate them
/// at fewer values.
///
/// In particular, sum(n <= x, (f <*> g)(n)) = sum(n <= x, ab = n, f(a)g(b)) = sum(ab <= x, f(a)g(b))
/// = sum(a <= A, b <= x/a, f(a)g(b)) + sum(b <= B, a <= x/b, f(a)g(b)) - sum(a <= A, b <= B, f(a)g(b))
/// = sum(a <= A, f(a)G(x/a)) + sum(b <= B, F(x/b)g(b)) - F(A)G(B)
/// where AB = x.
/// (see here https://gbroxey.github.io/blog/2023/04/30/mult-sum-1.html)
/// (see here https://angyansheng.github.io/blog/dirichlet-hyperbola-method)
/// (you may also find this instructive https://en.wikipedia.org/wiki/Marginal_distribution)
/// To accomplish this rearrangement, we first expanded (f <*> g) using the definition of dirichlet convolution.
/// Then, we interpreted the sum both as an a-indexed sum and as a b-indexed sum.
/// We could find the sum either way, but considering both sets us up to cut the work down quadratically.
/// These two ways of interpreting the sum are EXACTLY EQUIVALENT to the trick of transforming an integral of f(x) in terms of x into
/// an integral of f_inverse(y) in terms of y.  For integrals, it is common to just pick whichever interpretation is easier, but for our sum,
/// we benefit from realizing that if we sum f(a)g(b) for (a <= A, b <= x/a) and then separately for (b <= B, a <= x/b), then we double count the
/// rectangular region (a <= A, b <= B).
/// But more importantly, if we pick A = B = sqrt(x), then we've turned one sum over (f <*> g) for n from 1 to x, which we could accomplish in
/// about O(xlogx), into two sums over f*G and F*g for n from 1 to sqrt(x), which we can sometimes accomplish in O(sqrt(x)).
///
/// There is a big catch though, we now need to know values of F and G as well as f and g.
/// Wonderfully though, we only need to know F(x/n) and G(x/n) for integral values of n, and floor(x/n) has at most 2sqrt(x) distinct values.
/// That, arguably, is the real magic, since it lets us store only O(sqrt(x)) values of F and G.
/// When we want to compute only H(x), this is fine, but when we want to repeatedly compute the table of values of H(x/n) for integral values of n,
/// it is often better to make A larger than B, eg A = x**(2/3) and B = x**(1/3).

#include <stdbool.h>
#include <assert.h>

#include <nut/modular_math.h>

/// Wrapper to hold values of some multiplicative function.
/// Stores all values f(n) for n up to (and including) y, then
/// stores values f(x/n) for n less than x/y
typedef struct{
	int64_t x;
	int64_t y, yinv;
	int64_t *buf;
} nut_Diri;

/// Compute the sum of the divisor count function d from 1 to max.
/// @param [in] max: inclusive upper bound of range to compute sum for
/// @param [in] m: modulus to reduce result by, or 0 to skip reducing
/// @return the sum d(1) + ... + d(max)
NUT_ATTR_CONST
uint64_t nut_dirichlet_D(uint64_t max, uint64_t m);

/// Given a table of values of a multiplicative function f, compute (f <*> u)(x) for all x from 1 to n
/// See { @link nut_euler_sieve_conv_u} for more info}
/// @param [in] n: inclusive upper bound of range to compute f <*> u over
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_vals: table of values for f
/// @param [out] f_conv_u_vals: table to store values of f <*> u in
/// @return true on success, false on allocation failure
NUT_ATTR_NONNULL(3, 4)
NUT_ATTR_ACCESS(read_only, 3, 1)
NUT_ATTR_ACCESS(read_write, 4, 1)
bool nut_euler_sieve_conv_u(int64_t n, int64_t m, const int64_t f_vals[static n+1], int64_t f_conv_u_vals[restrict static n+1]);

/// Given a table of values of a multiplicative function f, compute (f <*> N)(x) for all x from 1 to n
/// See { @link nut_euler_sieve_conv_u} for more info}
/// @param [in] n: inclusive upper bound of range to compute f <*> N over
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_vals: table of values for f
/// @param [out] f_conv_N_vals: table to store values of f <*> N in
/// @return true on success, false on allocation failure
NUT_ATTR_NONNULL(3, 4)
NUT_ATTR_ACCESS(read_only, 3, 1)
NUT_ATTR_ACCESS(read_write, 4, 1)
bool nut_euler_sieve_conv_N(int64_t n, int64_t m, const int64_t f_vals[static n+1], int64_t f_conv_N_vals[restrict static n+1]);

/// Given tables of values of multiplicative functions f and g, compute (f <*> g)(x) for all x from 1 to n
/// This uses Euler's sieve, a variant of the sieve of Eratosthenes that only marks off each multiple once.
/// This generally has worse performance than the sieve of Eratosthenes, but some preliminary tests showed
/// that Eratosthenes is only about 14% faster in release mode than Euler, so currently only Euler's sieve
/// is implemented.
/// @param [in] n: inclusive upper bound of range to compute f <*> g over
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_vals: table of values for f
/// @param [in] g_vals: table of values for f
/// @param [out] f_conv_g_vals: table to store values of f <*> g in
/// @return true on success, false on allocation failure
NUT_ATTR_NONNULL(3, 4, 5)
NUT_ATTR_ACCESS(read_only, 3, 1)
NUT_ATTR_ACCESS(read_only, 4, 1)
NUT_ATTR_ACCESS(read_write, 5, 1)
bool nut_euler_sieve_conv(int64_t n, int64_t m, const int64_t f_vals[static n+1], const int64_t g_vals[static n+1], int64_t f_conv_g_vals[restrict static n+1]);

/// Allocate internal buffers for a diri table
/// self->buf will have f(0) through f(y) at indicies 0 through y,
/// and then f(x/1) through f(x/(yinv - 1)) at indicies y + 1 through y + yinv - 1.
/// The first two indicies, f(0) and f(1), are basically dummies though: f(0) should
/// never be relied on, and f(1) should always be 1
/// @param [in] x: inclusive upper bound of the domain of interest
/// @param [in] y: inclusive upper bound of the dense portion of the domain of interest.
/// Will be increased to sqrt(x) if needed, so 0 can be used to explicitly signal you want that behavior
/// @return true on success, false on allocation failure
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(write_only, 1)
bool nut_Diri_init(nut_Diri *self, int64_t x, int64_t y);

/// Copy the values from one diri table to another, which must be initialized
NUT_ATTR_NONNULL(1, 2)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 2)
void nut_Diri_copy(nut_Diri *restrict dest, const nut_Diri *restrict src);

/// Deallocate internal buffers for a diri table
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
void nut_Diri_destroy(nut_Diri *self);

NUT_ATTR_PURE
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_only, 1)
static inline int64_t nut_Diri_get_dense(const nut_Diri *self, int64_t k){
	assert(k >= 0 && k <= self->y);
	return self->buf[k];
}

NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
static inline void nut_Diri_set_dense(nut_Diri *self, int64_t k, int64_t v){
	assert(k >= 0 && k <= self->y);
	self->buf[k] = v;
}

NUT_ATTR_PURE
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_only, 1)
static inline int64_t nut_Diri_get_sparse(const nut_Diri *self, int64_t k){
	assert(k > 0 && k <= self->yinv);
	return self->buf[self->y + k];
}

NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
static inline void nut_Diri_set_sparse(nut_Diri *self, int64_t k, int64_t v){
	assert(k > 0 && k <= self->yinv);
	self->buf[self->y + k] = v;
}

/// Compute the value table for the dirichlet convolution identity I(n) = \{1 if n == 0, 0 otherwise\}
/// Just memset's the dense part, then sets index 1 and y + 1 through y + yinv - 1 to 1 (remember the sparse indicies are sums)
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
void nut_Diri_compute_I(nut_Diri *self);

/// Compute the value table for the unit function u(n) = 1
/// Fills the dense part of the table with 1s, and computes the sparse entries with table[y + k] = U(x/k) = x/k
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
void nut_Diri_compute_u(nut_Diri *self, int64_t m);

/// Compute the value table for the identity function N(n) = n
/// Fills the dense part of the table with increasing numbers, and computes the sparse entries with table[y + k] = sum_N(v = x/k) = v*(v+1)/2
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
void nut_Diri_compute_N(nut_Diri *self, int64_t m);

/// Compute the value table for the mobius function mu(n) (whose sum is called the Mertens function)
/// Requires both an initialized nut_Diri from {@link nut_Diri_init} and a packed table of mobius values from {@link nut_sieve_mobius}.
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] mobius: packed table of mobius values, from {@link nut_sieve_mobius} (with upper bound self->y).
NUT_ATTR_NONNULL(1, 3)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
void nut_Diri_compute_mertens(nut_Diri *restrict self, int64_t m, const uint8_t mobius[restrict static self->y/4 + 1]);

/// Compute the value table for the kth generalized divisor function dk(n)
/// dk(n) = da(n) <*> db(n) whenever a + b = k.
/// This function uses binary exponentiation to compute dk in only log(k) convolutions.
/// However, if all dk's are needed, calling { @link nut_Diri_compute_conv_u}
/// is more efficient
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistent with the inputs
/// @param [in] k: k, ie which generalized divisor function to compute.  dk is u <*> u <*> ... <*> u with k u's
/// @param [in] m: modulus to reduce the result by, or 0 to skip reducing
/// @param [in, out] f_tbl, g_tbl: temporary tables for scratch work, initialized by { @link nut_Diri_init}, fields y, x, and yinv must be set
/// Will still have scratch data stored on return
bool nut_Diri_compute_dk(nut_Diri *restrict self, uint64_t k, int64_t m, nut_Diri *restrict f_tbl, nut_Diri *restrict g_tbl);

/// Compute the value table for the kth power function N^k(n)
/// N^k(n) = n^k, so we can find the sums using https://en.wikipedia.org/wiki/Faulhaber%27s_formula
/// (we actually do this by taking a lower triangular matrix of pascal's triangle and inverting it)
/// This is one of the core components needed to compute the Dirichlet sum table for f where f(p) is a polynomial,
/// the other one being https://en.wikipedia.org/wiki/Jordan%27s_totient_function
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
/// @param [in] k: the power of the power function to compute
/// @param [in] m: modulus to reduce the result by, or 0 to skip reducing
bool nut_Diri_compute_Nk(nut_Diri *restrict self, uint64_t k, int64_t m);

/// Compute the value table for h = f <*> u, the dirichlet convolution of f and u (the unit function u(n) = 1), given the value table for f
/// See { @link nut_Diri_compute_conv } for details, this is just that function but with several specializations due to u being very simple.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistent with the inputs
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_tbl: table for the first operand
NUT_ATTR_NONNULL(1, 3)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
bool nut_Diri_compute_conv_u(nut_Diri *restrict self, int64_t m, const nut_Diri *restrict f_tbl);

/// Compute the value table for h = f <*> N, the dirichlet convolution of f and N (the identity function N(n) = n), given the value table for f
/// See { @link nut_Diri_compute_conv } for details, this is just that function but with several specializations due to N being very simple.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistent with the inputs
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_tbl: table for the first operand
NUT_ATTR_NONNULL(1, 3)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
bool nut_Diri_compute_conv_N(nut_Diri *restrict self, int64_t m, const nut_Diri *restrict f_tbl);

/// Compute the value table for h = f <*> g, the dirichlet convolution of f and g, given value tables for the operands
/// self must have been initialized using { @link nut_Diri_init}, and in particular the lengths and cutoffs for self, f_tbl, and g_tbl
/// must be set and match
/// This uses a sieve to compute h(1) through h(y), then uses dirichlet's hyperbola formula to compute H(x/1) through H(x/(yinv-1))
/// If one of the operands is a simple standard multiplicative function like the unit function u or the mobius function mu, then
/// try to use a specialized function from this library and use that instead.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistent with the inputs
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_tbl: table for the first operand (dirichlet convolution is commutative, so order doesn't matter)
/// @param [in] g_tbl: table for the second operand
NUT_ATTR_NONNULL(1, 3, 4)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
NUT_ATTR_ACCESS(read_only, 4)
bool nut_Diri_compute_conv(nut_Diri *restrict self, int64_t m, const nut_Diri *f_tbl, const nut_Diri *g_tbl);


/// Compute the value table for h such that f = g <*> h, aka h = f </> g where the division is in terms of dirichlet convolution
/// self must have been initialized using { @link nut_Diri_init }, and in particular the lengths and cutoffs for self, f_tbl, and g_tbl
/// must be set and match
/// This is the inverse of <*>, in the sense that if h = f </> g, then h <*> g = f.  In fact, this function is essentially
/// implemented by solving this for H(v) and applying the dirichlet hyperbola method (see the code comments for details).
/// Unlike the `compute_conv_*` family of functions, this needs to allocate space for scratch work, which can cause it to fail if
/// there is not enough memory.  About 16*(self->y + 1) bytes will be allocated.
/// Currently there are not specialized functions for dividing by simple standard functions.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init }.
/// self->y, self->x, and self->yinv must be set ant consistent with the inputs
/// @param [in] m: modulus to reduce results by, or 0 to skip reducing
/// @param [in] f_tbl: table for the first operand (the "numerator" aka "dividend" in the "division")
/// @param [in] g_tbl: the table for the second operand (the "denominator" aka "divisor" in the "division")
bool nut_Diri_convdiv(nut_Diri *restrict self, int64_t m, const nut_Diri *restrict f_tbl, const nut_Diri *restrict g_tbl);

