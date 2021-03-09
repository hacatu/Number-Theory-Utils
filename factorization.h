#pragma once

/// @file
/// @author hacatu
/// @version 0.1.0
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

/// Pretty name for signed 128 bit integer.
/// GCC implements 128 bit arithmetic in terms of 64 bit arithmetic.
typedef signed __int128 int128_t;
/// Pretty name for unsigned 128 bit integer.
/// GCC implements 128 bit arithmetic in terms of 64 bit arithmetic.
typedef unsigned __int128 uint128_t;

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
factors_t *init_factors_t_w(uint64_t max_primes);

/// Allocate a factors structure that can hold all factors of n, even if it has as many prime factors as possible.
/// 
/// E.g. 210=2*3*5*7 and 2310=2*3*5*7*11 so 210-2309 can have up to 4 distinct prime factors.
/// @param [in] n: number up to which the distinct prime factors of any number should be storable
/// @param [in] num_primes: number of primes in array
/// @param [in] primes: array of primes
/// @return pointer to factors structure that can hold enough distinct primes to factor any number up through n.  pointer needs to be free'd
factors_t *init_factors_t_ub(uint64_t n, uint64_t num_primes, const uint64_t *primes);

/// Compute nonnegative integral power of integer using binary exponentiation.
/// @param [in] b, e: base and exponent
/// @return b^e, not checked for overflow
static inline uint64_t pow_u64(uint64_t b, uint64_t e) __attribute__((const));
static inline uint64_t pow_u64(uint64_t b, uint64_t e){
	uint64_t r = 1;
	while(e){
		if(e&1){
			r = r*b;
		}
		e >>= 1;
		b *= b;
	}
	return r;
}

/// Multiply a factorization into a number.
/// @param [in] factors: pointer to factors struct, as obtained from { @link factor_trial_div}
/// @return product of prime powers described by factors
uint64_t factors_product(const factors_t *factors) __attribute__((pure));

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

/// Compute nonnegative integral power of a number modulo another using binary exponentiation.
/// @param [in] b, e, n: base, exponent, and modulus
/// @return b^e mod n, computed via binary exponentiation
uint64_t powmod(uint64_t b, uint64_t e, uint64_t n) __attribute__((const));

/// Check if n is prime using a deterministic Miller-Rabin test.
/// 7 partictular bases are used so that no composite number will falsely be reported as prime for the entire 64-bit range
/// @param [in] n: number to check for primality
/// @return true if n is prime, false otherwise
int is_prime_dmr(uint64_t n) __attribute__((const));

/// Generate a (pseudo)random integer uniformly from [a, b).
///
/// Currently calls { @link rand_u64} so this number is generated from /dev/random
/// but this function exists to provide a weaker, pseudorandom number source
/// if this turns out to be a bottleneck.
/// @param [in] a, b: bounds of the interval [a, b)
/// @return (pseudo)random integer uniformly chosen from [a, b)
uint64_t prand_u64(uint64_t a, uint64_t b);

/// Generate a (strong) random integer uniformly from [a, b).
///
/// Currently uses /dev/random via the getrandom function, but this is a Linux api.
/// @param [in] a, b: bounds of the interval [a, b)
/// @return (strong) random integer uniformly chosen from [a, b)
uint64_t rand_u64(uint64_t a, uint64_t b);

/// Factor out all powers of a given array of primes.
///
/// Primesieve is a good general source for primes, but the api for this function is designed
/// to allow giving a specialized list of primes like only primes equal to 1 mod 6 if it is
/// known that all factors have some form, or to provide an empty list when this function is
/// used in a more generic context like { @link factor_heuristic}.
/// @param [in] n: the number to factor
/// @param [in] num_primes: the number of primes in the array
/// @param [in] primes: the array of primes
/// @param [out] factors: pointer to struct where factors will be stored
/// @return n with all factors found and stored in factors divided out.  Thus if n factors completely over the given primes, 1 is returned.
uint64_t factor_trial_div(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes], factors_t *factors) __attribute__((pure));

/// Factor a number using a variety of approaches based on its size.
///
/// Currently a configuration struct must be passed.  Kraitcheck methods (quadratic sieve and number field sieve) are not implemented
/// because they are useless on 64 bit integers.  Parameters to Pollard-Rho-Brent with gcd aggregation and Lenstra ecf are not tuned
/// by this function.
/// @param [in] n: the number to factor
/// @param [in] num_primes: the number of primes in the array to try trial division on
/// @param [in] primes: array of primes
/// @param [in] conf: limits for different algorithms
/// @param [out] factors: output
/// @return n with all factors found and stored in factors divided out.  Thus if n factors completely, 1 is returned.
uint64_t factor_heuristic(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes], const factor_conf_t *conf, factors_t *factors);

/// Try to find a factor of a number using Pollard's Rho algorithm with Floyd cycle finding.
/// @param [in] n: number to find a factor of
/// @param [in] x: random value mod n
/// @return a nontrivial factor of n if found, 1 or n otherwise
uint64_t factor1_pollard_rho(uint64_t n, uint64_t x) __attribute__((const));

/// Try to find a factor of a number using Pollard's Rho algorithm with Brent cycle finding and gcd coalescing.
///
/// m does not affect whether x will lead to a factor of n, it just means each iteration will be
/// cheaper but up to 2m extraneous iterations could be performed.
/// @param [in] n: number to find a factor of
/// @param [in] x: random value mod n
/// @param [in] m: number of iterations per gcd check
/// @return a nontrivial factor of n if found, 1 or n otherwise
uint64_t factor1_pollard_rho_brent(uint64_t n, uint64_t x, uint64_t m) __attribute__((const));

/// Compute d = gcd(a, b) and x, y st. xa + by = d.
/// @param [in] a, b: numbers to find gcd of
/// @param [out] _t, _s: pointers to output x and y to respectively (ignored if NULL)
/// @return d
static inline int64_t egcd(int64_t a, int64_t b, int64_t *_t, int64_t *_s)
#ifndef DOXYGEN
__attribute__((access (write_only, 3), access (write_only, 4)))
#endif
;
static inline int64_t egcd(int64_t a, int64_t b, int64_t *_t, int64_t *_s){
	int64_t r0 = b, r1 = a;
	int64_t s0 = 1, s1 = 0;
	int64_t t0 = 0, t1 = 1;
	while(r1){
		int64_t q = r0/r1, t;
		t = r1;
		r1 = r0 - q*r1;
		r0 = t;
		t = s1;
		s1 = s0 - q*s1;
		s0 = t;
		t = t1;
		t1 = t0 - q*t1;
		t0 = t;
	}
	if(_t){
		*_t = t0;
	}
	if(_s){
		*_s = s0;
	}
	return r0;
}

/// Compute the Euclidean remainder r = a mod n for positive n so that 0 <= r < n.
/// @param [in] a, n: dividend and divisor
/// @return a mod n
static inline int64_t mod(int64_t a, int64_t n) __attribute__((const));
static inline int64_t mod(int64_t a, int64_t n){
	int64_t r = a%n;
	if(r < 0){
		r += n;
	}
	return r;
}

/// Compute n mod pq st n = a mod p and n = b mod q, where p and q are coprime.
/// @param [in] a, p, b, q: Chinese Remainder Theorem parameters.  The residues a and b should not be negative.
/// The moduli p and q should be coprime.
/// @return 0 <= 0 < pq so that n = a mod p and n = b mod q
static inline int64_t lift_crt(int64_t a, int64_t p, int64_t b, int64_t q) __attribute__((const));
static inline int64_t lift_crt(int64_t a, int64_t p, int64_t b, int64_t q){
	int64_t x, y;
	egcd(p, q, &x, &y);
	return mod(b*p%(p*q)*x + a*q%(p*q)*y, p*q);
}

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

/// Compute the Jacobi symbol of n mod k.
///
/// Uses modified euclidean algorithm.
/// @param [in] n, k: Jacobi symbol parameters
/// @return Jacobi symbol (0 if k | n, +1 if n is a quadratic residue mod an odd number of prime divisors of k (with multiplicity), -1 otherwise)
int64_t jacobi(int64_t n, int64_t k) __attribute__((const));

/// Compute a random number mod a prime that is not a quadratic residue.
///
/// This is useful in Shanks's algorithm and others.  This function works by rejection sampling using the Jacobi symbol as a test,
/// which succeeds in 2(p-2)/(p-1) trials on average.  If p is not prime, the Jacobi symbol can be positive for a nonresidue,
/// so not all nonresidues are possible, and the number of trials on average could be larger than the prime p case.
/// @param [in] p: the modulus for which to generate a nonresidue.  Should be prime.
/// @return a quadratic nonresidue mod p
int64_t rand_nr_u64(int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, this function is not guaranteed to terminate.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
int64_t sqrt_shanks(int64_t n, int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, the value returned may not be useful but the
/// function will terminate.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
int64_t sqrt_cipolla(int64_t n, int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, this function is not guaranteed to terminate.
/// If your use case does not guarantee this, call { @link is_prime_dmr} and { @link jacobi} beforehand.
/// uses shortcuts for primes that are 3, 5, or 7 mod 8, otherwise uses Shanks's algorithm
/// unless p-1 is divisible by a sufficiently high power of 2 so Cipolla's algorithm will be
/// faster, in which case it is used.  Only the Shank's branch can fail to terminate,
/// although other branches give useless results if the preconditions are not met.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
int64_t sqrt_mod(int64_t n, int64_t p);

