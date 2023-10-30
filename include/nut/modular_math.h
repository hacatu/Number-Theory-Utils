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
#include <assert.h>
#include <stdbool.h>

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

/// Wrapper for gcc's `access` function annotation - no clang equivalent, and currently disabled because it is bugged
/// See https://gcc.gnu.org/onlinedocs/gcc/Common-Function-Attributes.html
#if __has_c_attribute(gnu::access) && 0
#define NUT_ATTR_ACCESS(...) [[gnu::access(__VA_ARGS__)]]
#else
#define NUT_ATTR_ACCESS(...)
#endif

/// Prevent clang from reporting spurious errors when a length zero array is passed to a length annotated function parameter
#if __has_c_attribute(clang::no_sanitize)
#define NUT_ATTR_NO_SAN(name) [[clang::no_sanitize(name)]]
#else
#define NUT_ATTR_NO_SAN(name)
#endif

/// Compute nonnegative integral power of integer using binary exponentiation.
/// @param [in] b, e: base and exponent
/// @return b^e, not checked for overflow
[[gnu::const]]
uint64_t nut_u64_pow(uint64_t b, uint64_t e);

/// Compute nonnegative integral power of integer using binary exponentiation.
/// @param [in] b, e: base and exponent
/// @return b^e, not checked for overflow
[[gnu::const]]
uint128_t nut_u128_pow(uint128_t b, uint64_t e);

/// Compute nonnegative integral power of a number modulo another using binary exponentiation.
/// @param [in] b, e, n: base, exponent, and modulus
/// @return b^e mod n, computed via binary exponentiation
[[gnu::const]]
uint64_t nut_u64_powmod(uint64_t b, uint64_t e, uint64_t n);

/// Compute single binomial coefficient semi-naively.
/// Repeatedly does multiplications by n, n-1, n-2, ..., n-k+1 interleaved with divisions by 1, 2, 3, ..., k.
/// Overflows quickly, look into mod m versions if using large inputs.
/// See { @link nut_u64_binom_next} to iterate over values (n choose k), (n choose k+1).
/// See { @link nut_u64_binom_next_mod_2t} to iterate over values mod powers of 2
/// @param [in] n, k: Binomial coefficient arguments
[[gnu::const]]
uint64_t nut_u64_binom(uint64_t n, uint64_t k);

/// Compute the binomial coefficient for (n choose k) given the binomial coefficient for (n choose k-1)
/// Simply uses the recurrence (n choose k) = (n - k + 1)/k * (n choose k-1)
/// The starting point should be (n choose 0) = 1.  To start at k != 0, use
/// { @link nut_u64_binom}.
/// See { @link nut_u64_binom_next_mod_2t} to iterate over values mod powers of 2.
/// @param [in] n, k: Binomial coefficient arguments
/// @param [in] prev: Binomial coefficient value for (n choose k-1)
[[gnu::const]]
uint64_t nut_u64_binom_next(uint64_t n, uint64_t k, uint64_t prev);

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
NUT_ATTR_ACCESS(write_only, 3) NUT_ATTR_ACCESS(write_only, 4)
int64_t nut_i64_egcd(int64_t a, int64_t b, int64_t *restrict _t, int64_t *restrict _s);

/// Find the multiplicative inverse of a mod b
/// If b is a power of 2, use { @link nut_i64_modinv_2t}
/// Just a wrapper around { @link nut_i64_egcd}.
/// If a and b are not coprime, then the value
/// returned will just be the bezout coefficient for a, and will yield gcd(a, b)
/// when multiplied with a mod b instead of 1.
[[gnu::const]]
int64_t nut_i64_modinv(int64_t a, int64_t b);

/// Find the multiplicative inverse of a mod 2**t
/// This uses a hensel/newton like iterative algorithm, described here
/// https://crypto.stackexchange.com/a/47496
/// Basically, we use a lookup table to get the inverse mod 2**8, and then
/// use the fact that ax = 1 mod 2**k --> ax(2-ax) = 1 mod 2**(2k) to lift the inverse to
/// mod 2**16, mode 2**32, etc as needed.  This does cap out at 2**64.
/// @param [in] a: number to invert, MUST be odd.
/// @param [in] t: Exponent of modulus, ie we want to work mod 2**t, MUST be <= 64.  0 and 1 are allowed, but won't give very interesting results.
/// @return b such that a*b = 1 mod 2**t and 0 <= b < 2**t
[[gnu::const]]
uint64_t nut_u64_modinv_2t(uint64_t a, uint64_t t);

/// Compute the Euclidean remainder r = a mod n for positive n so that 0 <= r < n.
/// @param [in] a, n: dividend and divisor
/// @return a mod n
[[gnu::const]]
int64_t nut_i64_mod(int64_t a, int64_t n);

/// Compute n mod pq st n = a mod p and n = b mod q, where p and q are coprime.
/// @param [in] a, p, b, q: Chinese Remainder Theorem parameters.  The residues a and b should not be negative.
/// The moduli p and q should be coprime.
/// @return 0 <= 0 < pq so that n = a mod p and n = b mod q
[[gnu::const]]
int64_t nut_i64_crt(int64_t a, int64_t p, int64_t b, int64_t q);

/// Compute n mod pq st n = a mod p and n = b mod q, where p and q are coprime.
/// @param [in] a, p, b, q: Chinese Remainder Theorem parameters.  The residues a and b should not be negative.
/// The moduli p and q should be coprime.
/// @return 0 <= 0 < pq so that n = a mod p and n = b mod q
[[gnu::const]]
int128_t nut_i128_crt(int64_t a, int64_t p, int64_t b, int64_t q);

/// Compute the least common multiple of a and b
/// Divides the product by the gcd so can overflow for large arguments
/// @param [in] a, b: numbers to find nut_i64_lcm of
/// @return nut_i64_lcm(a, b)
[[gnu::const]]
int64_t nut_i64_lcm(int64_t a, int64_t b);

/// Compute the binomial coefficient for (n choose k) given the binomial coefficient for (n choose k-1)
/// Simply uses the recurrence (n choose k) = (n - k + 1)/k * (n choose k-1)
/// The starting point should be (n choose 0) = 1.  To start at k != 0, use
/// { @link nut_u64_binom}.
/// See { @link nut_u64_binom_next_mod_2t} to iterate over values mod powers of 2.
/// @param [in] n, k: Binomial coefficient arguments
/// @param [in] t: Exponent of modulus, ie we want to work mod 2^t
/// @param [in, out] v2: 2-adic valuation of (n choose k-1), that is, the highest power of 2 dividing (n choose k-1).
/// This will be updated in-place, so to start from (n choose 0) it can just be initialized to 0.
/// @param [in, out] p2: 2-coprime part of (n choose k-1), that is, (n choose k-1) divided by 2^v2,
/// This will be updated in-place, so to start from (n choose 0) it can just be initialized to 1.
/// @return (n choose k) mod 2^t, which is equal to 2^v2 * p2 mod 2^t
[[nodiscard, gnu::nonnull(4, 5)]]
NUT_ATTR_ACCESS(read_write, 4) NUT_ATTR_ACCESS(read_write, 5)
uint64_t nut_u64_binom_next_mod_2t(uint64_t n, uint64_t k, uint64_t t, uint64_t *restrict v2, uint64_t *restrict p2);

/// Compute the Jacobi symbol of n mod k.
///
/// Uses modified euclidean algorithm.
/// @param [in] n, k: Jacobi symbol parameters
/// @return Jacobi symbol (0 if k | n, +1 if n is a quadratic residue mod an odd number of prime divisors of k (with multiplicity), -1 otherwise)
[[gnu::const]]
int64_t nut_i64_jacobi(int64_t n, int64_t k);

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
[[gnu::const]]
int64_t nut_i64_sqrt_shanks(int64_t n, int64_t p);

/// Compute the square root of a quadratic residue mod a prime.
///
/// If n is not a residue or p is not a prime, the value returned may not be useful but the
/// function will terminate.
/// @param [in] n: a quadratic residue mod p
/// @param [in] p: a prime
/// @return r so that r^2 = n mod p
[[gnu::const]]
int64_t nut_i64_sqrt_cipolla(int64_t n, int64_t p);

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
[[gnu::const]]
int64_t nut_i64_sqrt_mod(int64_t n, int64_t p);


/// Compute a constant to use to do modular division faster
/// @param [in] pd: absolute value of divisor
/// @return constant c for use with { @link nut_i32_fastmod}
[[gnu::const]]
uint64_t nut_i32_fastmod_init(uint32_t pd);

/// Compute a constant to use to do modular division faster
/// @param [in] d: divisor
/// @return constant c for use with { @link nut_u32_fastmod}
[[gnu::const]]
uint64_t nut_u32_fastmod_init(uint32_t d);

/// Compute a constant to use to do modular division faster
/// @param [in] pd: absolute value of divisor
/// @return constant c for use with { @link nut_i64_fastmod}
[[gnu::const]]
uint128_t nut_i64_fastmod_init(uint64_t pd);

/// Compute a constant to use to do modular division faster
/// @param [in] d: divisor
/// @return constant c for use with { @link nut_u64_fastmod}
[[gnu::const]]
uint128_t nut_u64_fastmod_init(uint64_t d);

/// Compute n mod d, choosing the smallest signed remainder
/// @param [in] n: dividend
/// @param [in] pd: absolute value of divisor
/// @param [in] c: constant from { @link nut_i32_fastmod_init}
/// @return n mod pd
[[gnu::const]]
int32_t nut_i32_fastmod_trunc(int32_t n, uint32_t pd, uint64_t c);

/// Compute n mod d, choosing the euclidean remainder
/// @param [in] n: dividend
/// @param [in] pd: absolute value of divisor
/// @param [in] c: constant from { @link nut_i32_fastmod_init}
/// @return n mod pd
[[gnu::const]]
int32_t nut_i32_fastmod_floor(int32_t n, uint32_t pd, uint64_t c);

/// Compute n mod d faster using a precomputed constant
/// @param [in] n: dividend
/// @param [in] d: divisor
/// @param [in] c: constant from { @link nut_u32_fastmod_init}
/// @return n mod d
[[gnu::const]]
uint32_t nut_u32_fastmod(uint32_t n, uint32_t d, uint64_t c);

/// Compute n mod d, choosing the smallest signed remainder
/// @param [in] n: dividend
/// @param [in] pd: absolute value of divisor
/// @param [in] c: constant from { @link nut_i64_fastmod_init}
/// @return n mod pd
[[gnu::const]]
int64_t nut_i64_fastmod_trunc(int64_t n, uint64_t pd, uint128_t c);

/// Compute n mod d, choosing the euclidean remainder
/// @param [in] n: dividend
/// @param [in] pd: absolute value of divisor
/// @param [in] c: constant from { @link nut_i64_fastmod_init}
/// @return n mod pd
[[gnu::const]]
int64_t nut_i64_fastmod_floor(int64_t n, uint64_t pd, uint128_t c);

/// Compute n mod d faster using a precomputed constant
/// @param [in] n: dividend
/// @param [in] d: divisor
/// @param [in] c: constant from { @link nut_u64_fastmod_init}
/// @return n mod d
[[gnu::const]]
uint64_t nut_u64_fastmod(uint64_t n, uint64_t d, uint128_t c);

