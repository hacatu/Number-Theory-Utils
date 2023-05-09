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

/// Wrapper to hold values of some multiplicative function.
/// Stores all values f(n) for n up to (and including) y, then
/// stores values f(x/n) for n less than x/y
typedef struct{
	int64_t x;
	int64_t y, yinv;
	int64_t *buf;
} diri_table;

/// Compute the sum of the divisor count function d from 1 to max.
/// @param [in] max: inclusive upper bound of range to compute sum for
/// @return the sum d(1) + ... + d(max)
uint64_t dirichlet_D(uint64_t max);

/// Given a table of values of a multiplicative function f, compute (f <*> u)(x) for all x from 1 to n
/// Currently this uses Euler's sieve, sometimes called a linear sieve, but this is almost certainly slower than
/// a sieve of Eratosthenes in practice
/// @param [in] n: inclusive upper bound of range to compute f <*> u over
/// @param [in] f_vals: table of values for f
/// @param [out] f_conv_u_vals: table to store values of f <*> u in
/// @return true on success, false on allocation failure
bool euler_sieve_conv_u(int64_t n, const int64_t f_vals[static n+1], int64_t f_conv_u_vals[static n+1]);

/// Allocate internal buffers for a diri table
/// self->buf will have f(0) through f(y) at indicies 0 through y,
/// and then f(x/1) through f(x/(yinv - 1)) at indicies y + 1 through y + yinv - 1.
/// The first two indicies, f(0) and f(1), are basically dummies though: f(0) should
/// never be relied on, and f(1) should always be 1
/// @param [in] x: inclusive upper bound of the domain of interest
/// @param [in] y: inclusive upper bound of the dense portion of the domain of interest.
/// Will be increased to sqrt(x) if needed, so 0 can be used to explicitly signal you want that behavior
/// @return true on success, false on allocation failure
bool diri_table_init(diri_table *self, int64_t x, int64_t y);

/// Deallocate internal buffers for a diri table
void diri_table_destroy(diri_table *self);

static inline int64_t diri_table_get_dense(const diri_table *self, int64_t k){
	return self->buf[k];
}

static inline void diri_table_set_dense(diri_table *self, int64_t k, int64_t v){
	self->buf[k] = v;
}

static inline int64_t diri_table_get_sparse(const diri_table *self, int64_t k){
	return self->buf[self->y + k];
}

static inline void diri_table_set_sparse(diri_table *self, int64_t k, int64_t v){
	self->buf[self->y + k] = v;
}

/// Compute the value table for the dirichlet convolution identity I(n) = \{1 if n == 0, 0 otherwise\}
/// Just memset's the dense part, then sets index 1 and y + 1 through y + yinv - 1 to 1 (remember the sparse indicies are sums)
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
void compute_I_diri_table(diri_table *self);

/// Compute the value table for the unit function u(n) = 1
/// Fills the dense part of the table with 1s, and computes the sparse entries with table[y + k] = U(x/k) = x/k
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
void compute_u_diri_table(diri_table *self);

/// Compute the value table for the identity function N(n) = n
/// Fills the dense part of the table with increasing numbers, and computes the sparse entries with table[y + k] = sum_N(v = x/k) = v*(v+1)/2
/// @param [in, out] self: the table to store the result in, and take the bounds from.  Must be initialized
void compute_N_diri_table(diri_table *self);

/// Compute the value table for h = f <*> u, the dirichlet convolution of f and u (the unit function u(n) = 1), given the value table for f
/// See { @link compute_conv_diri_table } for details, this is just that function but with several specializations due to u being very simple.
/// @param [in, out] self: the table to store the result in, initialized by { @link diri_table_init}.
/// self->y, self->x, and self->yinv must be set and consistend with the inputs
/// @param [in] f_tbl: table for the first operand
bool compute_conv_u_diri_table(diri_table *self, const diri_table *f_tbl);

/// Compute the value table for h = f <*> g, the dirichlet convolution of f and g, given value tables for the operands
/// self must have been initialized using { @link diri_table_init}, and in particular the lengths and cutoffs for self, f_tbl, and g_tbl
/// must be set and match
/// This uses a sieve to compute h(1) through h(y), then uses dirichlet's hyperbola formula to compute H(x/1) through H(x/(yinv-1))
/// If one of the operands is a simple standard multiplicative function like the unit function u or the mobius function mu, then
/// try to use a specialized function from this library and use that instead.
/// @param [in, out] self: the table to store the result in, initialized by { @link diri_table_init}.
/// self->y, self->x, and self->yinv must be set and consistend with the inputs
/// @param [in] f_tbl: table for the first operand (dirichlet convolution is commutative, so order doesn't matter)
/// @param [in] g_tbl: table for the second operand
bool compute_conv_diri_table(diri_table *self, const diri_table *f_tbl, const diri_table *g_tbl);

