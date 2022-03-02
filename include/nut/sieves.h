#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Sieve based functions for computing primes in a range or divisor counts/
/// power sums/Euler's totient function/Carmichael's function on all numbers
/// in a range.

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

#include <nut/factorization.h>

/// Compute the maximum number of unique prime divisors a number can have.
/// This has nothing to do with factoring the number and is just a simple binary
/// decision diagram based on comparing the number to 2, 2*3, 2*3*5, etc.
/// @param [in] max: number to find max unique prime divisors of
/// @return the largest number of unique prime divisors any number not exceeding max can possibly have
uint64_t max_prime_divs(uint64_t max) __attribute__((const));

/// Get an entry from a variable pitch array.
/// In a normal array, the element type is a complete type with size known at compile time, or even a variably
/// modified type for which sizeof and ordinary array operations will work as desired.
/// However, it is sometimes useful to have arrays whose elements are structs with flexible length array members,
/// and because flexible length arrays do not have their lengths automatically tracked like static or variable length arrays,
/// we can't use sizeof or ordinary array operations.  In particular, sizeof returns the padded size assuming the flexible
/// length array has zero length, and array operations work as if we had an array of structs where the flexible length array
/// has zero length.  This function instead takes the base pointer of the array, the pitch, and the index, and returns a pointer
/// to a member without bounds checking.
/// @param [in] buf: pointer to the start of the array
/// @param [in] pitch: offset from start of one element to start of next element (currently should be computed as
/// offsetof(type, fla_member) + w*fla_element_size, but a more complex computation should be used if alignment is critically important)
/// @param [in] i: index of element to get
/// @return pointer to the i-th member of an array with given base and pitch
static inline void *pitch_arr_get(void *buf, size_t pitch, uint64_t i) __attribute__((const));
static inline void *pitch_arr_get(void *buf, size_t pitch, uint64_t i){
	return buf + i*pitch;
}

/// Compute the factorization for every number in the range from 0 to max.
/// The factorizations for 0 and 1 are not actually computed.  The factorizations
/// are stored in an array of factors_t structs with capacity w, where w is the maximum
/// number of unique prime divisors of a number not exceeding max.  Thus the pitch of this
/// array is offsetof(factors_t, factors) + *_w*sizeof(dummy->factors[0]),
/// where dummy is some expression with type factors_t*.
/// {@link pitch_arr_get} should be used to handle the returned value.
/// @param [in] max: inclusive upper bound of sieving range in which to factor all numbers
/// @param [out] _w: store w, the maximum number of unique prime divisors of a number not exceeding max
/// @return a pointer to an array of factors_t structs containing the factorization of all numbers not exceeding max
void *sieve_factorizations(uint64_t max, uint64_t *_w);

/// Compute the unique prime factors of every number in the range from 0 to max.
/// The factors for 0 and 1 are not actually computed.  The result is stored in
/// an array of fw_u64arr_t structs with capacity w, where w is the maximum number of unique
/// prime divisors of a number not exceeding max.  Thus the pitch of this array is
/// offsetof(fw_u64arr_t, elems) + *_w*sizeof(uint64_t).
/// {@link pitch_arr_get} should be used to handle the returned value.
/// @param [in] max: inclusive upper bound of sieving range in which to find unique prime factors of all numbers
/// @param [out] _w: maximum numer of unique prime divisors of a number not exceeding max
/// @return a pointer to an array of fw_u64arr_t structs containing lists of unique prime factors for all numbers not exceeding max
void *sieve_factors(uint64_t max, uint64_t *_w);

/// Compute the number of divisors (including 1 and n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_count} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find divisor counts for all numbers
/// @return an array of divisor counts for all numbers in the range
uint64_t *sieve_sigma_0(uint64_t max);

/// Compute the sum of divisors (including 1 and n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_sum} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find divisor sums for all numbers
/// @return an array of divisor sums for all numbers in the range
uint64_t *sieve_sigma_1(uint64_t max);

/// Compute the sum of some power of divisors (including 1 and n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_power_sum} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find divisor power sums for all numbers
/// @param [in] e: power of divisors for summing, eg 0 would produce divisor counts, 1 divisor sums, 2 sums of squares of divisors, etc
/// @return an array of divisor power sums for all numbers in the range
uint64_t *sieve_sigma_e(uint64_t max, uint64_t e);

/// Compute Euler's totient function for every number from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link euler_phi} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find totients for all numbers
/// @return an array of totients for all numbers in the range
uint64_t *sieve_phi(uint64_t max);

/// Compute the Carmichael function for every number from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link carmichael_lambda} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to compute Carmichael for all numbers
/// @return an array of Carmichael function outputs for all numbers in the range
uint64_t *sieve_carmichael(uint64_t max);

/// Compute an array of all primes from 0 to max.
/// @param [in] max: inclusive upper bound of sieving range
/// @param [out] _num_primes: how many primes were found in the range (this pointer cannot be null)
/// @return an array of all primes from 0 to max.
uint64_t *sieve_primes(uint64_t max, uint64_t *_num_primes);

