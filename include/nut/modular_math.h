#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Functions for dealing with modular arithmetic on
/// 64 bit integers (not suitable for large integers which may overflow,
/// around 2^31.  a bignum-enabled version may be created to handle this)

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

/// Pretty name for signed 128 bit integer.
/// GCC implements 128 bit arithmetic in terms of 64 bit arithmetic.
typedef signed __int128 int128_t;
/// Pretty name for unsigned 128 bit integer.
/// GCC implements 128 bit arithmetic in terms of 64 bit arithmetic.
typedef unsigned __int128 uint128_t;

/// A fixed capacity array.
/// The capacity is fixed and the length in use is stored before the array proper.
/// Not intended to be reallocated, hence the absence of a cap field.  Instead intended for an array of these to be returned,
/// for instance from a sieve that finds all distinct prime factors of all numbers in a range.  Such an array is variable
/// pitch and so must be handled as a void pointer with accessor functions for convenience, see {@link sieves.h}.
typedef struct{
	uint64_t len;
	uint64_t elems[];
} nut_u64_Pitcharr;

/// Compute nonnegative integral power of integer using binary exponentiation.
/// @param [in] b, e: base and exponent
/// @return b^e, not checked for overflow
__attribute__((const)) static inline uint64_t nut_u64_pow(uint64_t b, uint64_t e){
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

/// Compute nonnegative integral power of integer using binary exponentiation.
/// @param [in] b, e: base and exponent
/// @return b^e, not checked for overflow
__attribute__((const)) static inline uint128_t nut_u128_pow(uint128_t b, uint64_t e){
	uint128_t r = 1;
	while(e){
		if(e&1){
			r = r*b;
		}
		e >>= 1;
		b *= b;
	}
	return r;
}

/// Compute nonnegative integral power of a number modulo another using binary exponentiation.
/// @param [in] b, e, n: base, exponent, and modulus
/// @return b^e mod n, computed via binary exponentiation
__attribute__((const)) uint64_t nut_u64_powmod(uint64_t b, uint64_t e, uint64_t n);

/// Compute single binomial coefficient semi-naively.
/// Repeatedly does multiplications by n, n-1, n-2, ..., n-k+1 interleaved with divisions by 1, 2, 3, ..., k.
/// There are better ways to do this if we need to know huge binomial coefficients mod some number, or find many
/// binomial coefficients with the same n or k, but this function is sufficient a lot of the time
__attribute__((const)) static inline uint64_t nut_u64_binom(uint64_t n, uint64_t k){
	uint64_t res = 1;
	if(n - k < k){
		k = n - k;
	}
	for(uint64_t i = 0; i < k; ++i){
		res = res*(n - i)/(1 + i);
	}
	return res;
}

/// Generate a (pseudo)random integer uniformly from [a, b).
///
/// Currently calls { @link nut_u64_rand} so this number is generated from /dev/random
/// but this function exists to provide a weaker, pseudorandom number source
/// if this turns out to be a bottleneck.
/// @param [in] a, b: bounds of the interval [a, b)
/// @return (pseudo)random integer uniformly chosen from [a, b)
uint64_t nut_u64_prand(uint64_t a, uint64_t b);

/// Generate a (strong) random integer uniformly from [a, b).
///
/// Currently uses /dev/random via the getrandom function, but this is a Linux api.
/// @param [in] a, b: bounds of the interval [a, b)
/// @return (strong) random integer uniformly chosen from [a, b)
uint64_t nut_u64_rand(uint64_t a, uint64_t b);

/// Compute d = gcd(a, b) and x, y st. xa + by = d.
/// @param [in] a, b: numbers to find gcd of
/// @param [out] _t, _s: pointers to output x and y to respectively (ignored if NULL)
/// @return d
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wunknown-attributes"
#endif
__attribute__((access(write_only, 3), access(write_only, 4))) static inline int64_t nut_i64_egcd(int64_t a, int64_t b, int64_t *_t, int64_t *_s){
#pragma GCC diagnostic pop
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
__attribute__((const)) static inline int64_t nut_i64_mod(int64_t a, int64_t n){
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
__attribute__((const)) static inline int64_t nut_i64_crt(int64_t a, int64_t p, int64_t b, int64_t q){
	int64_t x, y;
	nut_i64_egcd(p, q, &x, &y);
	return nut_i64_mod(b*p%(p*q)*x + a*q%(p*q)*y, p*q);
}

/// Compute the least common multiple of a and b
/// Divides the product by the gcd so can overflow for large arguments
/// @param [in] a, b: numbers to find nut_i64_lcm of
/// @return nut_i64_lcm(a, b)
__attribute__((const)) static inline int64_t nut_i64_lcm(int64_t a, int64_t b){
	return a*b/nut_i64_egcd(a, b, NULL, NULL);
}

/// Compute the Jacobi symbol of n mod k.
///
/// Uses modified euclidean algorithm.
/// @param [in] n, k: Jacobi symbol parameters
/// @return Jacobi symbol (0 if k | n, +1 if n is a quadratic residue mod an odd number of prime divisors of k (with multiplicity), -1 otherwise)
__attribute__((const)) int64_t nut_i64_jacobi(int64_t n, int64_t k);

/// Compute a random number mod a prime that is not a quadratic residue.
///
/// This is useful in Shanks's algorithm and others.  This function works by rejection sampling using the Jacobi symbol as a test,
/// which succeeds in 2(p-2)/(p-1) trials on average.  If p is not prime, the Jacobi symbol can be positive for a nonresidue,
/// so not all nonresidues are possible, and the number of trials on average could be larger than the prime p case.
/// @param [in] p: the modulus for which to generate a nonresidue.  Should be prime.
/// @return a quadratic nonresidue mod p
int64_t nut_i64_rand_nr_mod(int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, this function is not guaranteed to terminate.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
int64_t nut_i64_sqrt_shanks(int64_t n, int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, the value returned may not be useful but the
/// function will terminate.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
__attribute__((const)) int64_t nut_i64_sqrt_cipolla(int64_t n, int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, this function is not guaranteed to terminate.
/// If your use case does not guarantee this, call { @link nut_u64_is_prime_dmr} and { @link nut_i64_jacobi} beforehand.
/// uses shortcuts for primes that are 3, 5, or 7 mod 8, otherwise uses Shanks's algorithm
/// unless p-1 is divisible by a sufficiently high power of 2 so Cipolla's algorithm will be
/// faster, in which case it is used.  Only the Shank's branch can fail to terminate,
/// although other branches give useless results if the preconditions are not met.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
int64_t nut_i64_sqrt_mod(int64_t n, int64_t p);

