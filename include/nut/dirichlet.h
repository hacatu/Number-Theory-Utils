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

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

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
/// @return the sum d(1) + ... + d(max)
uint64_t nut_dirichlet_D(uint64_t max);

/// Given a table of values of a multiplicative function f, compute (f <*> u)(x) for all x from 1 to n
/// Currently this uses Euler's sieve, sometimes called a linear sieve, but this is almost certainly slower than
/// a sieve of Eratosthenes in practice
/// @param [in] n: inclusive upper bound of range to compute f <*> u over
/// @param [in] f_vals: table of values for f
/// @param [out] f_conv_u_vals: table to store values of f <*> u in
/// @return true on success, false on allocation failure
bool nut_euler_sieve_conv_u(int64_t n, const int64_t f_vals[static n+1], int64_t f_conv_u_vals[static n+1]);

/// Given a table of values of a multiplicative function f, compute (f <*> N)(x) for all x from 1 to n
/// Currently this uses Euler's sieve, sometimes called a linear sieve, but this is almost certainly slower than
/// a sieve of Eratosthenes in practice
/// @param [in] n: inclusive upper bound of range to compute f <*> N over
/// @param [in] f_vals: table of values for f
/// @param [out] f_conv_N_vals: table to store values of f <*> N in
/// @return true on success, false on allocation failure
bool nut_euler_sieve_conv_N(int64_t n, const int64_t f_vals[static n+1], int64_t f_conv_N_vals[static n+1]);

/// Given tables of values of multiplicative functions f and g, compute (f <*> g)(x) for all x from 1 to n
/// Currently this uses Euler's sieve, sometimes called a linear sieve, but this is almost certainly slower than
/// a sieve of Eratosthenes in practice
/// @param [in] n: inclusive upper bound of range to compute f <*> g over
/// @param [in] f_vals: table of values for f
/// @param [in] g_vals: table of values for f
/// @param [out] f_conv_g_vals: table to store values of f <*> g in
/// @return true on success, false on allocation failure
bool nut_euler_sieve_conv(int64_t n, const int64_t f_vals[static n+1], const int64_t g_vals[static n+1], int64_t f_conv_g_vals[static n+1]);

/// Allocate internal buffers for a diri table
/// self->buf will have f(0) through f(y) at indicies 0 through y,
/// and then f(x/1) through f(x/(yinv - 1)) at indicies y + 1 through y + yinv - 1.
/// The first two indicies, f(0) and f(1), are basically dummies though: f(0) should
/// never be relied on, and f(1) should always be 1
/// @param [in] x: inclusive upper bound of the domain of interest
/// @param [in] y: inclusive upper bound of the dense portion of the domain of interest.
/// Will be increased to sqrt(x) if needed, so 0 can be used to explicitly signal you want that behavior
/// @return true on success, false on allocation failure
bool nut_Diri_init(nut_Diri *self, int64_t x, int64_t y);

/// Deallocate internal buffers for a diri table
void nut_Diri_destroy(nut_Diri *self);

static inline int64_t nut_Diri_get_dense(const nut_Diri *self, int64_t k){
	assert(k >= 0 && k <= self->y);
	return self->buf[k];
}

static inline void nut_Diri_set_dense(nut_Diri *self, int64_t k, int64_t v){
	assert(k >= 0 && k <= self->y);
	self->buf[k] = v;
}

static inline int64_t nut_Diri_get_sparse(const nut_Diri *self, int64_t k){
	assert(k > 0 && k <= self->yinv);
	return self->buf[self->y + k];
}

static inline void nut_Diri_set_sparse(nut_Diri *self, int64_t k, int64_t v){
	assert(k > 0 && k <= self->yinv);
	self->buf[self->y + k] = v;
}

/// Compute the value table for the dirichlet convolution identity I(n) = \{1 if n == 0, 0 otherwise\}
/// Just memset's the dense part, then sets index 1 and y + 1 through y + yinv - 1 to 1 (remember the sparse indicies are sums)
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
void nut_Diri_compute_I(nut_Diri *self);

/// Compute the value table for the unit function u(n) = 1
/// Fills the dense part of the table with 1s, and computes the sparse entries with table[y + k] = U(x/k) = x/k
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
void nut_Diri_compute_u(nut_Diri *self);

/// Compute the value table for the identity function N(n) = n
/// Fills the dense part of the table with increasing numbers, and computes the sparse entries with table[y + k] = sum_N(v = x/k) = v*(v+1)/2
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
void nut_Diri_compute_N(nut_Diri *self);

/// Compute the value table for the mobius function mu(n) (whose sum is called the Mertens function)
/// Requires both an initialized nut_Diri from {@link nut_Diri_init} and a packed table of mobius values from {@link nut_sieve_mobius}.
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
/// @param [in] mobius: packed table of mobius values, from {@link nut_sieve_mobius} (with upper bound self->y).
void nut_Diri_compute_mertens(nut_Diri *self, const uint8_t mobius[self->y/4 + 1]);

/// Compute the value table for h = f <*> u, the dirichlet convolution of f and u (the unit function u(n) = 1), given the value table for f
/// See { @link nut_Diri_compute_conv } for details, this is just that function but with several specializations due to u being very simple.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistent with the inputs
/// @param [in] f_tbl: table for the first operand
bool nut_Diri_compute_conv_u(nut_Diri *self, const nut_Diri *f_tbl);

/// Compute the value table for h = f <*> N, the dirichlet convolution of f and N (the identity function N(n) = n), given the value table for f
/// See { @link nut_Diri_compute_conv } for details, this is just that function but with several specializations due to N being very simple.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistent with the inputs
/// @param [in] f_tbl: table for the first operand
bool nut_Diri_compute_conv_N(nut_Diri *self, const nut_Diri *f_tbl);

/// Compute the value table for h = f <*> g, the dirichlet convolution of f and g, given value tables for the operands
/// self must have been initialized using { @link nut_Diri_init}, and in particular the lengths and cutoffs for self, f_tbl, and g_tbl
/// must be set and match
/// This uses a sieve to compute h(1) through h(y), then uses dirichlet's hyperbola formula to compute H(x/1) through H(x/(yinv-1))
/// If one of the operands is a simple standard multiplicative function like the unit function u or the mobius function mu, then
/// try to use a specialized function from this library and use that instead.
/// @param [in, out] self: the table to store the result in, initialized by { @link nut_Diri_init}.
/// self->y, self->x, and self->yinv must be set and consistend with the inputs
/// @param [in] f_tbl: table for the first operand (dirichlet convolution is commutative, so order doesn't matter)
/// @param [in] g_tbl: table for the second operand
bool nut_Diri_compute_conv(nut_Diri *self, const nut_Diri *f_tbl, const nut_Diri *g_tbl);

