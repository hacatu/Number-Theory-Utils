#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Functions for dealing with factorization and modular arithmetic on
/// 64 bit integers (not suitable for large integers which may overflow,
/// around 2^31.  a bignum-enabled version may be created to handle this)

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include <nut/modular_math.h>

/// The most distinct prime divisors a uint32_t can have.
#define MAX_PRIMES_32 9
/// The most distinct prime divisors a uint64_t can have.
#define MAX_PRIMES_64 15
/// The most distinct prime divisors a uint128_t can have.
#define MAX_PRIMES_128 25

/// Currently unused.
/// How many trial divisions to do before performing randomized prime tests (Miller Rabin).
/// For 64 bit numbers, we can easily make Miller Rabin deterministic in 7 tests, whereas up to 40
/// would be used to check a large number.  BSW and ECPP are also important for checking large primes,
/// because BSW and Miller Rabin together can provide a higher primality confidence faster, whereas
/// ECPP is much slower but provides an easily verifiable proof.
#define DIVS_BEFORE_PRIME_CHECK 9

/// A prime factorization.
/// Stores a list of prime, power pairs in a flexible length array, along with the number of primes in the factorization.
/// Not intended to be reallocated, hence the absence of a cap field.
typedef struct{
	uint64_t num_primes;
	struct{
		uint64_t prime;
		uint64_t power;
	} factors[];
} factors_t;

/// Configuration/tuning for heuristic factorization.
typedef struct{
	/// Use Pollard-Rho-Brent when factoring numbers at most this large.
	/// 0 to disable Pollard, UINT64_MAX to always use Pollard (note that Pollard is much slower than Lenstra for numbers
	/// larger than 2^20)
	uint64_t pollard_max;
	/// How many iterations of Pollard-Rho-Brent to do in between gcd tests.
	/// 100 or lower is reasonable, this has no effect on whether or not factors will be found, it merely speeds
	/// up iterations at the cost of potentially doing up to 2m extra iterations.
	uint64_t pollard_stride;
	/// Use lenstra ecf when factoring numbers at most this large.
	/// 0 to disable lenstra ecf, UINT64_MAX to always use lenstra above pollard_max.
	/// Currently this is the fastest method implemented for large numbers so any number larger than this might not
	/// be fully factored.  Also note that the factoring algorithms provided may fail with overflow for numbers larger
	/// than 2^30.
	uint64_t lenstra_max;
	/// Factorial smoothness bound for Lenstra ecf.
	/// P <- 2P through P <- kP are computed when searching for a nontrival factor, where k is this bound.
	uint64_t lenstra_bfac;
	/// Currently unused.
	/// Will be the last number for which the quadratic sieve is used, after which the number field sieve is used.
	uint64_t qsieve_max;
} factor_conf_t;

/// Allocate a factors structure that can hold a given number of distinct primes.
/// @param [in] max_primes: the number of distinct primes that should be storable
/// @return pointer to factors structure that can hold max_primes.  pointer needs to be free'd
factors_t *init_factors_t_w(uint64_t max_primes) __attribute__((malloc));

/// Allocate a factors structure that can hold all factors of n, even if it has as many prime factors as possible.
/// 
/// E.g. 210=2*3*5*7 and 2310=2*3*5*7*11 so 210-2309 can have up to 4 distinct prime factors.
/// @param [in] n: number up to which the distinct prime factors of any number should be storable
/// @param [in] num_primes: number of primes in array
/// @param [in] primes: array of primes
/// @return pointer to factors structure that can hold enough distinct primes to factor any number up through n.  pointer needs to be free'd
factors_t *init_factors_t_ub(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes]) __attribute__((malloc));

/// Allocate a copy of a factors struct, but with only enough memory to store its current factors and not its max capacity if higher
/// @param [in] factors: the struct to copy
/// @return a copy of the input or NULL on allocation failure
factors_t *copy_factors_t(const factors_t *factors) __attribute__((malloc));

/// Multiply a factorization into a number.
/// @param [in] factors: pointer to factors struct, as obtained from { @link factor_trial_div}
/// @return product of prime powers described by factors
uint64_t factors_product(const factors_t *factors) __attribute__((pure));

/// Find the number of divisors of a number given its prime factorization, including itself and 1.
/// Works by multiplying power + 1 for all prime factors
/// @param [in] factors: pointer to factors struct, as obtained from { @link factor_trial_div}
/// @return number of divisors
uint64_t divisor_count(const factors_t *factors) __attribute__((pure));

/// Find the sum of divisors of a number given its prime factorization, including itself and 1.
/// Works by multiplying (prime**(power+1)-1)/(prime-1) for all prime factors
/// @param [in] factors: pointer to factors struct, as obtained from { @link factor_trial_div}
/// @return sum of divisors
uint64_t divisor_sum(const factors_t *factors) __attribute__((pure));

/// Find the sum of powers of divisors of a number given its prime factorization, including itself and 1.
/// Note this is NOT the same as the sum of divisors of a number.
/// Works by multiplying (prime**((power+1)*e)-1)/(prime**e-1) for all prime factors, where power is the
/// power of each prime and e is the power of divisors to sum
/// @param [in] factors: pointer to factors struct, as obtained from { @link factor_trial_div}
/// @param [in] power: power of divisors to sum
/// @return sum of divisor powers
uint64_t divisor_power_sum(const factors_t *factors, uint64_t power) __attribute__((pure));

/// Raise a factorization to a power, ie multiply all exponents by a constant.
/// @param [in,out] factors: factorization to raise to a power
/// @param [in] power: power to raise the factorization to
void factors_power(factors_t *factors, uint64_t power);

/// Find Euler's Phi function, the number of coprime numbers less than n.
/// @param [in] factors: the factorization of n for which to compute phi
/// @return phi(n)
uint64_t euler_phi(const factors_t *factors) __attribute__((pure));

/// Find Carmichael's Lambda function, the smallest exponent m so a^m = 1 for all 0 < a < n.
/// Always divides phi(n).
/// @param [in] factors: the factorization of n for which to compute the carmichael function
/// @return lambda(n)
uint64_t carmichael_lambda(const factors_t *factors) __attribute__((pure));

/// Call { @link forall_divisors} with temporarily allocated dfactors and pfactors structs.
int forall_divisors_tmptmp(const factors_t *factors, int (*f)(const factors_t*, uint64_t, void*), void *data);

/// Call a given function on each divisor of a number, given its factorization.
/// @param [in] factors: the factorization of n for which to compute all divisors
/// @param [in] f: callback function.  Arguments are factorization of divisor, divisor, user data (respectively).
/// Should return 0 to continue, or 1 to break.
/// @param [in,out] data: user data (allows arbitrary data to be passed in and out of the callback)
/// @param [in] dfactors: buffer to store internal work.  Should be allocated to hold at least factors->num_primes
/// @param [in] pfactors: buffer to store internal work.  Should be allocated to hlod at least factors->num_primes
/// @return 0 if the callback never returned 1 and all divisors were visited, or 1 if the callback ever returned nonzero.
int forall_divisors(const factors_t *factors, int (*f)(const factors_t*, uint64_t, void*), void *data, factors_t *dfactors, factors_t *pfactors);

/// Call { @link forall_divisors_le} with temporarily allocated dfactors and pfactors structs.
int forall_divisors_le_tmptmp(const factors_t *factors, uint64_t d_max, int (*f)(const factors_t*, uint64_t, void*), void *data);

/// Call a given function on each divisor of a number, given its factorization, skipping divisors over a given bound.
/// @param [in] factors: the factorization of n for which to compute all divisors
/// @param [in] d_max: max divisor to visit.  This should be less than the product of factors, otherwise just use forall_divisors
/// @param [in] f: callback function.  Arguments are factorization of divisor, divisor, user data (respectively).
/// Should return 0 to continue, or 1 to break.
/// @param [in,out] data: user data (allows arbitrary data to be passed in and out of the callback)
/// @param [in] dfactors: buffer to store internal work.  Should be allocated to hold at least factors->num_primes
/// @param [in] pfactors: buffer to store internal work.  Should be allocated to hlod at least factors->num_primes
/// @return 0 if the callback never returned 1 and all divisors were visited, or 1 if the callback ever returned nonzero.
int forall_divisors_le(const factors_t *factors, uint64_t d_max, int (*f)(const factors_t*, uint64_t, void*), void *data, factors_t *dfactors, factors_t *pfactors);

/// Print a factorization of a number.
/// @param [in,out] file: pointer to file to print to
/// @param [in] factors: pointer to factorization struct
/// @return number of characters printed
int factors_fprint(FILE *file, const factors_t *factors);

/// Add a power of some prime to an existing factorization struct.
///
/// If the input factorization is sorted before this function is called,
/// it will still be sorted after, and if the factor to be added is already
/// present its power will be increased rather than creating a separate entry.
/// The factorization struct must have enough space if an additional entry
/// is needed.  Composite numbers or zero powers should not be supplied.
/// @param [in,out] factors: pointer to factorization struct
/// @param [in] m: prime
/// @param [in] k: power
void factors_append(factors_t *factors, uint64_t m, uint64_t k);

/// Add all prime/power pairs in one factorization to another.
///
/// If both inputs are sorted, the result will be sorted and no dupliate entries
/// will be created.  Acts as if { @link factors_append} were called repeatedly.
/// @param [in,out] factors: pointer to factorization struct which should be extended
/// @param [in] factors2: pointer to factorization struct whose entries will be added to the other struct
/// @param [in] k: power to multiply entries of factors2 by, for if some composite number m was discovered where m^k divides n and we have factored m
void factors_combine(factors_t *factors, const factors_t *factors2, uint64_t k);

/// Check if n is prime using a deterministic Miller-Rabin test.
/// 7 partictular bases are used so that no composite number will falsely be reported as prime for the entire 64-bit range
/// @param [in] n: number to check for primality
/// @return true if n is prime, false otherwise
int is_prime_dmr(uint64_t n) __attribute__((const));

/// Factor out all powers of a given array of primes.
///
/// Primesieve is a good general source for primes, but the api for this function is designed
/// to allow giving a specialized list of primes like only primes equal to 1 mod 6 if it is
/// known that all factors have some form, or to provide an empty list when this function is
/// used in a more generic context like {@link factor_heuristic}.
/// @param [in] n: the number to factor
/// @param [in] num_primes: the number of primes in the array
/// @param [in] primes: the array of primes
/// @param [out] factors: pointer to struct where factors will be stored
/// @return n with all factors found and stored in factors divided out.  Thus if n factors completely over the given primes, 1 is returned.
uint64_t factor_trial_div(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes], factors_t *factors) __attribute__((pure));

/// Factor a number using a variety of approaches based on its size.
///
/// If conf->pollard_max is greater than 3, the primes array should include at least 2 and 5 to
/// avoid an infinite loop since {@link factor1_pollard_rho} can't factor 4 or 25.
/// Currently a configuration struct must be passed.  Kraitcheck methods (quadratic sieve and number field sieve) are not implemented
/// because they are useless on 64 bit integers.  Parameters to Pollard-Rho-Brent with gcd aggregation and Lenstra ecf are not tuned
/// by this function.
/// @param [in] n: the number to factor
/// @param [in] num_primes: the number of primes in the array to try trial division on (use 25 if using {@link small_primes})
/// @param [in] primes: array of primes (can use {@link small_primes})
/// @param [in] conf: limits for different algorithms (can use {@link default_factor_conf})
/// @param [out] factors: output
/// @return n with all factors found and stored in factors divided out.  Thus if n factors completely, 1 is returned.
uint64_t factor_heuristic(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes], const factor_conf_t *conf, factors_t *factors);

/// Get the floor of the nth root of a
///
/// Uses bitscan to estimate the base 2 log and in turn nth root of a, then uses Newton's method
/// @param [in] a: the number to take the nth root of
/// @param [in] n: root to take
/// @return floor(a**(1/n)), ie the largest integer x such that x**n <= a
uint64_t u64_nth_root(uint64_t a, uint64_t n) __attribute__((const));

/// Check if a is a perfect power of some integer.
///
/// Used in {@link factor_heuristic}.
/// max is the max exponent a could possibly be, which is limited because it is a uint64_t.
/// If nothing is known about a, then the max exponent it could be is 63, if it were 2**63.
/// However, if we know that the base of a must be larger than some minimum,
/// such as if we know a is not divisible by any primes up to some point,
/// then we can set a lower limit on the max exponent.
/// For example, if a is 2-rough (not divisible by any primes <= 2), then the minimum possible base is 3, so the maximal exponent is 40 instead.
/// For all primes in the small_primes array, the max exponents are as follows:
/// 2-rough -> min base 3 -> max exponent 40
/// 3 -> 5 -> 27
/// 5 -> 7 -> 22
/// 7 -> 11 -> 18
/// 11 -> 13 -> 17
/// 13 -> 17 -> 15
/// 19 -> 23 -> 14
/// 23 -> 29 -> 13
/// 29 -> 31 -> 12
/// 37 -> 41 -> 11
/// 57 -> 59 -> 10
/// 83 -> 89 -> 9
/// and going over 100:
/// 137 -> 139 -> 8
/// 251 -> 257 -> 7
/// 563 -> 569 -> 6
/// 1621 -> 1627 -> 5
/// 7129 -> 7151 -> 4
/// 65521 -> 65537 -> 3
/// 2642239 -> 2642257 -> 2
/// 4294967291 -> 4294967311 -> 1
/// The max exponent can also be brought down if a is known to be small, but depending on context it may or may not be worth computing
/// log(a)/log(p) where p is the largest prime not known to not divide a.
/// @param [in] a: the number to check for being a perfect power
/// @param [in] max: the max exponent a could possibly be.  Use eg 63 if unknown, 11 if a has had all factors of primes up to 47 removed.
/// @param [out] _base: if a is found to be a perfect power, store the base here
/// @param [out] _exp: if a is found to be a perfect power, store the exponent here
/// @return true if a is a perfect power (with exponent max or lower), false otherwise (could mean a is not a perfect power, or could mean a is a perfect power with exponent max + 1 to 61)
bool is_perfect_power(uint64_t a, uint64_t max, uint64_t *_base, uint64_t *_exp);

/// Try to find a factor of a number using Pollard's Rho algorithm with Floyd cycle finding.
/// Note that this will not find factors of 4 or 25 no matter what x is.
/// @param [in] n: number to find a factor of
/// @param [in] x: random value mod n
/// @return a nontrivial factor of n if found, 1 or n otherwise
uint64_t factor1_pollard_rho(uint64_t n, uint64_t x) __attribute__((const));

/// Try to find a factor of a number using Pollard's Rho algorithm with Brent cycle finding and gcd coalescing.
/// Note that his will not find factors of 4 or 25 no matter what x is.
/// m does not affect whether x will lead to a factor of n, it just means each iteration will be
/// cheaper but up to 2m extraneous iterations could be performed.
/// @param [in] n: number to find a factor of
/// @param [in] x: random value mod n
/// @param [in] m: number of iterations per gcd check
/// @return a nontrivial factor of n if found, 1 or n otherwise
uint64_t factor1_pollard_rho_brent(uint64_t n, uint64_t x, uint64_t m) __attribute__((const));

/// Try to find a factor of a number using Lenstra ecf.
///
/// In particular, an affine wierstrass representation is used internally
/// and scalar multiplication is done using a binary "double and add" algorithm.
/// @param [in] n: number to find a factor of
/// @param [in] x, y, a: random values mod n
/// @param [in] B: number of trials before giving up (we compute kP for k from 2 to B)
/// @return a nontrivial factor of n if found, 1 or n otherwise
int64_t factor1_lenstra(int64_t n, int64_t x, int64_t y, int64_t a, int64_t B) __attribute__((const));

/// Same as { @link factor1_lenstra} but using a projective Montgomery curve and Montgomery ladder.
///
/// Instead of an affine Wierstrass curve and binary multiplication.  This is supposed to be
/// faster but seems to be about the same.
/// @param [in] n: number to find a factor of
/// @param [in] x, y, a: random numbers mod n
/// @param [in] B: number of trials before giving up (we compute kP for k from 2 to B)
/// @return a nontrivial factor of n if found, 1 or n otherwise
int64_t factor1_lenstra_montgomery(int64_t n, int64_t x, int64_t y, int64_t a, int64_t B) __attribute__((const));



/// An array containing the 25 primes up to 100.
/// {@link factor1_pollard_rho} and {@link factor1_pollard_rho_brent} can't find factors of 4 or 25 using
/// the default polynomial, so if using {@link factor_heuristic} at least
/// these primes should be used for trial division if the configuration allows pollard rho to be called.
extern const uint64_t small_primes[25];

/// Default value for {@link factor_heuristic}
/// Use if you don't want to tune the parameters or don't care
extern const factor_conf_t default_factor_conf;

