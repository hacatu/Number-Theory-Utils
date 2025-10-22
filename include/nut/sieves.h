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
/// Arrays returned from sieve functions should be freed when no longer needed.

#include <inttypes.h>
#include <stdlib.h>

#include <nut/modular_math.h>
#include <nut/factorization.h>

/// Compute the maximum number of unique prime divisors a number can have.
/// This has nothing to do with factoring the number and is just a simple binary
/// decision diagram based on comparing the number to 2, 2*3, 2*3*5, etc.
/// @param [in] max: number to find max unique prime divisors of
/// @return the largest number of unique prime divisors any number not exceeding max can possibly have
NUT_ATTR_CONST
uint64_t nut_max_prime_divs(uint64_t max);

/// Compute an upper bound on the number of primes up to max.
/// This uses an inequality involving log derived from the prime number theorem to always get an
/// upper bound.
/// @param [in] max: number to find the number of primes up to
/// @param an upper bound on the number of primes up to max
NUT_ATTR_CONST
uint64_t nut_max_primes_le(uint64_t max);

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
NUT_ATTR_CONST
NUT_ATTR_NONNULL(1)
NUT_ATTR_RETURNS_NONNULL
NUT_ATTR_ACCESS(none, 1)
static inline void *nut_Pitcharr_get(void *buf, size_t pitch, uint64_t i){
	return buf + i*pitch;
}

/// Get a bit from a bitarray.
/// Simply does buf[i/8] & (1ull << (i%8)).
/// @param [in] buf: pointer to bitarray
/// @param [in] i: index of element to get
/// @return false if i-th element is false, true otherwise
NUT_ATTR_PURE
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_only, 1)
static inline bool nut_Bitarray_get(const uint8_t *buf, uint64_t i){
	return buf[i/8] & (UINT8_C(1) << (i%8));
}

/// Set a bit in a bitarray
/// @param [in, out] buf: pointer to bitarray
/// @param [in] i: index of element to set
/// @param [in] v: true to set ith bit, false to clear
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
static inline void nut_Bitarray_set(uint8_t *buf, uint64_t i, bool v){
	if(v){
		buf[i/8] |= 1ull << (i%8);
	}else{
		buf[i/8] &= ~(UINT8_C(1) << (i%8));
	}
}

/// Get an element from an array of bitfields of length 2, aka uint2's.
/// The result will be shifted to the least significant position, that is,
/// only 0, 1, 2, and 3 are possible results.
/// @param [in] buf: pointer to array of bitfields
/// @param [in] i: index of element to get
/// @return i-th element (0, 1, 2, or 3)
NUT_ATTR_PURE
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_only, 1)
static inline uint8_t nut_Bitfield2_arr_get(const uint8_t *buf, uint64_t i){
	return (buf[i/4] >> (i%4*2)) & 3;
}

/// Calculate factorials and inverse factorials for a given upper bound and modulus
/// @param [in] k: factorials[bits + k - 1] is the last factorial that will be computed
/// @param [in] modulus: modulus to reduce result by.  Must be large enough that all inv factorials are actually invertable
/// @param [in] bits: used with k to find the last factorial to compute
/// @param [in] max_denom: inv_factorials[max_denom] is the last one that will be computed
/// @param [out] factorials: output for factorial table
/// @param [out] inv_factorials: output for inverse factorial table
/// @return true on success, false on failure (if the inverse of some factorial can't be found)
bool nut_u64_make_factorial_tbl(uint64_t k, uint64_t modulus, uint64_t bits, uint64_t max_denom, uint64_t factorials[static bits + k], uint64_t inv_factorials[static max_denom + 1]);

/// Compute the factorization for every number in the range from 0 to max.
/// The factorizations for 0 and 1 are not actually computed.  The factorizations
/// are stored in an array of factors_t structs with capacity w, where w is the maximum
/// number of unique prime divisors of a number not exceeding max.  The result is
/// a pitched array.
/// {@link nut_Pitcharr_get} should be used to handle the returned value.
/// {@link get_factorizations_pitch} should be used to get the pitch from the output
/// parameter w.
/// @param [in] max: inclusive upper bound of sieving range in which to factor all numbers
/// @param [out] _w: store w, the maximum number of unique prime divisors of a number not exceeding max
/// @return a pointer to an array of factors_t structs containing the factorization of all numbers not exceeding max,
/// or NULL on allocation failure
[[deprecated("Replace with nut_sieve_smallest_factors_wheel6")]]
NUT_ATTR_NONNULL(2)
NUT_ATTR_MALLOC
NUT_ATTR_ACCESS(write_only, 2)
void *nut_sieve_factorizations(uint64_t max, uint64_t *_w);

/// Get the pitch for a pitched array of factorization structs with w unique prime divisors.
/// This is simply offsetof(factors_t, factors) + w*sizeof(dummy->factors[0]),
/// where dummy is an expression with type factors_t.
/// @param [in] w: The maximal number of unique prime divisors of any potential index of the pitched array.
/// Should be obtained from {@link nut_sieve_factorizations}, {@link max_prime_divisors}, etc.
/// @return The pitch of a pitched array of factorization structs whose flexible length members all have w
/// elements.
NUT_ATTR_CONST
uint64_t nut_get_factorizations_pitch(uint64_t w);

/// Compute the unique prime factors of every number in the range from 0 to max.
/// The factors for 0 and 1 are not actually computed.  The result is stored in
/// an array of nut_u64_Pitcharr structs with capacity w, where w is the maximum number of unique
/// prime divisors of a number not exceeding max.  The pitch of the result may be obtained with
/// {@link get_factors_pitch}.
/// {@link nut_Pitcharr_get} should be used to handle the returned value.
/// @param [in] max: inclusive upper bound of sieving range in which to find unique prime factors of all numbers
/// @param [out] _w: maximum numer of unique prime divisors of a number not exceeding max
/// @return a pointer to an array of nut_u64_Pitcharr structs containing lists of unique prime factors for all numbers not exceeding max,
/// or NULL on allocation failure
[[deprecated("Replace with nut_sieve_smallest_factors_wheel6")]]
NUT_ATTR_NONNULL(2)
NUT_ATTR_MALLOC
NUT_ATTR_ACCESS(write_only, 2)
void *nut_sieve_factors(uint64_t max, uint64_t *_w);

/// Compute the number of distinct prime divisors of every number in the range from 0 to max.
/// The divisors for 0 and 1 are not actually computed.
/// @param [in] max: inclusive upper bound of sieving range in which to find distinct prime divisor counts of all numbers
/// @return a pointer to an array of distinct prime divisor counts for all numbers not exceeding max,
/// or NULL on allocation failure
NUT_ATTR_MALLOC
uint8_t *nut_sieve_omega(uint64_t max);

/// Compute the largest prime factor of every number in the range from 0 to max.
/// Compared to {@link nut_sieve_factors}, this uses up to 30 times less memory, so if
/// the range is very large, most numbers will never actually have their factorizations
/// accessed, or the factorizations will only be accessed about once, this function
/// should be preferred.  The factors of 0 and 1 are not actually computed, their entries in
/// the returned array will be 0.  You can use the resulting table as it is, or convert it to
/// a factorization using {@link nut_fill_factors_from_largest}.
/// @param [in] max: inclusive upper bound of sieving range in which to find largest prime factors of all numbers
/// @return a pointer to an array of largest prime factors for all numbers not exceeding max, or NULL on allocation failure.
[[deprecated("Replace with nut_sieve_smallest_factors_wheel6")]]
NUT_ATTR_MALLOC
uint64_t *nut_sieve_largest_factors(uint64_t max);

/// Use a table of largest prime factors to get the factorization of a number
/// @param [out] out: Factors struct to store result in.  MUST be allocated already, use {@link nut_make_Factors_ub} or
/// {@link nut_max_prime_divs} and {@link nut_make_Factors_w} if needed.
/// @param [in] n: the number to get the factorization of
/// @param [in] largest_factors: table of largest factors, from {@link nut_sieve_largest_factors}
[[deprecated("Replace with nut_fill_factors_from_smallest_wheel6")]]
NUT_ATTR_NONNULL(1, 3)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
void nut_fill_factors_from_largest(nut_Factors *restrict out, uint64_t n, const uint64_t largest_factors[restrict static n + 1]);

/// Compute the smallest prime factor of every number in the range from 0 to max, or 1 for primes
/// Compared to {@link nut_sieve_factors}, this uses up to 60 times less memory, so if
/// the range is very large, most numbers will never actually have their factorizations
/// accessed, or the factorizations will only be accessed about once, this function
/// should be preferred.  The factors of 0 and 1 are not actually computed, their entries in
/// the returned array will be 0.  You can use the resulting table as it is, or convert it to
/// a factorization using {@link nut_fill_factors_from_smallest}.
/// This function is able to store factors as 32 bit integers instead of 64, since the smallest prime factor of a
/// composite number is at most its square root.  This is why we store 1 for primes instead of themselves
/// @param [in] max: inclusive upper bound of sieving range in which to find smallest prime factors of all numbers
/// @return a pointer to an array of smallest prime factors for all numbers not exceeding max, or NULL on allocation failure.
[[deprecated("Replace with nut_sieve_smallest_factors_wheel6")]]
NUT_ATTR_MALLOC
uint32_t *nut_sieve_smallest_factors(uint64_t max);

/// Use a table of smallest prime factors to get the factorization of a number
/// @param [out] out: Factors struct to store result in.  MUST be allocated already, use {@link nut_make_Factors_ub} or
/// {@link nut_max_prime_divs} and {@link nut_make_Factors_w} if needed.
/// @param [in] n: the number to get the factorization of
/// @param [in] smallest_factors: table of smallest factors, from {@link nut_sieve_smallest_factors}
[[deprecated("Replace with nut_fill_factors_from_smallest_wheel6")]]
NUT_ATTR_NONNULL(1, 3)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
void nut_fill_factors_from_smallest(nut_Factors *restrict out, uint64_t n, const uint32_t smallest_factors[restrict static n + 1]);

/// Compute the smallest prime factor of every number coprime to 6 in the range from 0 to max, or 1 for primes
/// Compared to {@link nut_sieve_factors}, this uses up to 180 times less memory, so if
/// the range is very large, most numbers will never actually have their factorizations
/// accessed, or the factorizations will only be accessed about once, this function
/// should be preferred.  The factors of 1 are not actually computed, its entry in
/// the returned array will be 0.  You can use the resulting table as it is, or convert it to
/// a factorization using {@link nut_fill_factors_from_smallest_wheel6}.
/// This function is able to store factors as 32 bit integers instead of 64, since the smallest prime factor of a
/// composite number is at most its square root.  This is why we store 1 for primes instead of themselves
/// @param [in] max: inclusive upper bound of sieving range in which to find smallest prime factors of all numbers
/// @return a pointer to an array of smallest prime factors for all numbers not exceeding max, or NULL on allocation failure.
NUT_ATTR_MALLOC
uint32_t *nut_sieve_smallest_factors_wheel6(uint64_t max);

/// Use a table of smallest prime factors for a wheel of 6 to get the factorization of a number
/// @param [out] out: Factors struct to store result in.  MUST be allocated already, use {@link nut_make_Factors_ub} or
/// {@link nut_max_prime_divs} and {@link nut_make_Factors_w} if needed.
/// @param [in] n: the number to get the factorization of
/// @param [in] smallest_factors: table of smallest factors, from {@link nut_sieve_smallest_factors}
NUT_ATTR_NONNULL(1, 3)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 3)
void nut_fill_factors_from_smallest_wheel6(nut_Factors *restrict out, uint64_t n, const uint32_t smallest_factors[restrict static n/3 + 1]);

/// Get the pitch for a pitched array of {@link nut_u64_Pitcharr} factor lists.
/// This is simply offsetof(nut_u64_Pitcharr, elems) + *_w*sizeof(uint64_t).
/// @param [in] w: The maximal number of unique prime divisors of any potential index of the pitched array.
/// Should be obtained from {@link nut_sieve_factors}, {@link max_prime_divisors}, etc.
/// @return The pitch of a pitched array of factor list structs whose flexible length members all have w
/// elements.
NUT_ATTR_CONST
uint64_t nut_get_factors_pitch(uint64_t w);

/// Compute the number of divisors (including 1 and n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_count} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find divisor counts for all numbers
/// @return an array of divisor counts for all numbers in the range, or NULL on allocation failure
NUT_ATTR_MALLOC
uint64_t *nut_sieve_sigma_0(uint64_t max);

/// Compute the sum of divisors (including 1 and n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_sum} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find divisor sums for all numbers
/// @return an array of divisor sums for all numbers in the range, or NULL on allocation failure
NUT_ATTR_MALLOC
uint64_t *nut_sieve_sigma_1(uint64_t max);

/// Compute the sum of some power of divisors (including 1 and n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_power_sum} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find divisor power sums for all numbers
/// @param [in] e: power of divisors for summing, eg 0 would produce divisor counts, 1 divisor sums, 2 sums of squares of divisors, etc
/// @return an array of divisor power sums for all numbers in the range, or NULL on allocation failure
NUT_ATTR_MALLOC
uint64_t *nut_sieve_sigma_e(uint64_t max, uint64_t e);

/// Compute the generalized divisor function dk(n) (number of k-tuples with product n) for every number n from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link divisor_tuple_count}
/// for all numbers in the range, but without needing to compute or store the factorizations intermediately.
/// Note that dk is multiplicative so dk(mn) = dk(m)dk(n) when m and n are coprime, and dk(p^a) = binom(a + k, k) for prime powers.
/// In other words, dk is exponential in k and this function will overflow if max^k is too large.
/// @param [in] max: inclusive upper bound of sieving range in which to compute generalized divisor function
/// @param [in] k: number of factors per factorization, eg for a prime power p^a we get binom(a + k, k).
/// @param [in] modulus: modulus to reduce results by, or zero to skip reducing
/// @return an array of dk results, or NULL on allocation failure
NUT_ATTR_MALLOC
uint64_t *nut_sieve_dk(uint64_t max, uint64_t k, uint64_t modulus);

/// Compute Euler's totient function for every number from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link euler_phi} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to find totients for all numbers
/// @return an array of totients for all numbers in the range, or NULL on allocation failure
NUT_ATTR_MALLOC
uint64_t *nut_sieve_phi(uint64_t max);

/// Compute the Carmichael function for every number from 0 to max.
/// The results for 0 and 1 are not actually computed.  This effectively computes {@link carmichael_lambda} for
/// all numbers in the range, but without needing to compute or store the factorizations intermediately
/// @param [in] max: inclusive upper bound of sieving range in which to compute Carmichael for all numbers
/// @return an array of Carmichael function outputs for all numbers in the range, or NULL on allocation failure
NUT_ATTR_MALLOC
uint64_t *nut_sieve_carmichael(uint64_t max);

/// Compute the Mobius function for every number from 0 to max.
/// The result is stored in an array of 2 bit integers, which should
/// be accessed by {@link nut_Bitfield2_arr_get}.  That function will return 0 for 0,
/// 1 for 1, and 3 for -1.  2 will not be stored anywhere in bounds in the resulting array.
/// @param [in] max: inclusive upper bound of sieving range in which to compute Mobius for all numbers
/// @return a bitfield array of Mobius function outputs for all numbers in the range, with 3 instead of -1,
/// or NULL on allocation failure
NUT_ATTR_MALLOC
uint8_t *nut_sieve_mobius(uint64_t max);

/// Compute the Mertens function (sum of Mobius function) for every number from 0 to max.
/// Note that this function is signed.
/// @param [in] max: inclusive upper bound of range in which to compute Mertens for all numbers
/// @param [in] mobius: bitfield array of Mobius function outputs (from {@link nut_sieve_mobius}).
/// @return an array of Mertens function outputs for all numbers in the range, or NULL on allocation failure
NUT_ATTR_MALLOC
NUT_ATTR_NONNULL(2)
NUT_ATTR_ACCESS(read_only, 2)
int64_t *nut_compute_mertens_range(uint64_t max, const uint8_t mobius[static max/4 + 1]);

/// Compute a bitarray of whether or not each number from 0 to max is composite.
/// 1 is composite, and 0 is considered composite here.
/// The result should be used with {@link nut_is_composite} since it is packed (only stores bitflags for numbers coprime to 30).
/// @param [in] max: inclusive upper bound of sieving range in which to check compositeness for all numbers
/// @param [out] _num_primes: the number of primes in the range will be stored here.  May be NULL.
/// @return a bitarray of whether or not each number in the range is composite, or NULL on allocation failure
NUT_ATTR_MALLOC
uint8_t *nut_sieve_is_composite(uint64_t max);

/// Check if a number is composite using a packed bitarray from {@link nut_sieve_is_composite}
/// @param [in] n: the number to check if composite
/// @param [in] buf: packed bitarray from {@link nut_sieve_is_composite}
/// @return true if n is composite, false if n is prime
NUT_ATTR_PURE
NUT_ATTR_NONNULL(2)
NUT_ATTR_ACCESS(read_only, 2)
bool nut_is_composite(uint64_t n, const uint8_t buf[static n/30 + 1]);

/// Compute the pi (prime counting) function for every number from 0 to max.
/// The result should be used with {@link nut_compute_pi_from_tables} since pi is only actually calculated at every 240th number
/// since intermediate results can be computed with a single popcount on the packed buf bitarray.
/// @param [in] max: inclusive upper bound of range to compute pi function
/// @param [in] buf: packed bitarray from {@link nut_sieve_is_composite}
/// @return an array of pi values at every 240th number (use {@link nut_compute_pi_from_tables})
NUT_ATTR_NONNULL(2)
NUT_ATTR_MALLOC
NUT_ATTR_ACCESS(read_only, 2)
uint64_t *nut_compute_pi_range(uint64_t max, const uint8_t buf[static max/30 + 1]);

/// Get the value for the pi (prime counting) function for a particular number using precomputed tables.
/// @param [in] n: the number to calculate pi for
/// @param [in] pi_table: array of partial pi values from {@link nut_compute_pi_range}
/// @param [in] buf: packed bitarray from {@link nut_sieve_is_composite}
/// @return the number of primes <= n
NUT_ATTR_PURE
NUT_ATTR_NONNULL(2, 3)
NUT_ATTR_ACCESS(read_only, 2)
NUT_ATTR_ACCESS(read_only, 3)
NUT_ATTR_NO_SAN("vla-bound")
uint64_t nut_compute_pi_from_tables(uint64_t n, const uint64_t pi_table[restrict static n/30], const uint8_t buf[restrict static n/30 + 1]);

/// Compute an array of all primes from 0 to max.
/// @param [in] max: inclusive upper bound of sieving range
/// @param [out] _num_primes: how many primes were found in the range (this pointer cannot be null)
/// @return an array of all primes from 0 to max, or NULL on allocation failure
NUT_ATTR_MALLOC
NUT_ATTR_NONNULL(2)
NUT_ATTR_ACCESS(write_only, 2)
uint64_t *nut_sieve_primes(uint64_t max, uint64_t *_num_primes);

