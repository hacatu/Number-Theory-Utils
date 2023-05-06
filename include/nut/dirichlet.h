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

