#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Functions for dealing with polynomials with integer coefficients, especially
/// over finite fields.

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include <nut/modular_math.h>

/// A polynomial with signed integer coefficients.
/// Uses a dense representation, ideal for low degree polynomials and polynomials with almost no
/// zero coefficients.  The coefficient array goes from lowest power to highest power (constant
/// term up).  For the zero polynomial, len should be 1 and the constant coefficient should be 0.
/// For all other polynomials, len should be 1 greater than the degree of the polynomial.
/// In particular, the coefficient on the highest term should not be 0 ({ @link nut_Poly_normalize}).
typedef struct{
	/// Length of the polynomial.
	/// 1 greater than the degree for nonzero polynomials, 1 for the zero polynomial.
	uint64_t len;
	/// Maximum length polynomial the coeffs buffer can currently store.
	uint64_t cap;
	/// Array of coefficients.
	/// Should be managed by { @link nut_Poly_init}, { @link nut_Poly_destroy}, and { @link nut_Poly_ensure_cap}.
	int64_t *coeffs;
} nut_Poly;

/// The roots of a polynomial over a finite field, ignoring multiplicity.
typedef struct{
	/// Number of roots, can be 0 if there are none.
	uint64_t len;
	/// Maximum number of roots the roots buffer can currently store.
	uint64_t cap;
	/// Array of roots.
	/// Should be managed by { @link nut_Roots_init}, { @link nut_Roots_destroy}, and { @link nut_Roots_ensure_cap}.
	int64_t *roots;
} nut_Roots;



/// Initialize a polynomial with at least a given capacity.
///
/// If the given capacity is 0, it defaults to 4.  The polynomial is set to the zero polynomial.
/// @param [out] f: polynomial to initialize
/// @param [in] reserve: initial capacity, defaults to 4 if 0 is given
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Poly_init(nut_Poly *f, uint64_t reserve);

/// Free resources held by a polynomial struct.
/// @param [in,out] f: polynomial to destroy
[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
void nut_Poly_destroy(nut_Poly *f);

/// Asymptotically compare polynomials.
///
/// Does not assume the inputs are normalized.
/// Finds the highest degree term without matching coefficients.
/// @param [in] a, b: polynomials to compare
/// @return -1 if a < b, 1 if a > b, 0 if a == b
[[gnu::pure, gnu::nonnull(1, 2)]]
NUT_ATTR_ACCESS(read_only, 1) NUT_ATTR_ACCESS(read_only, 2)
int nut_Poly_cmp(const nut_Poly *a, const nut_Poly *b);

/// Initialize a root buffer with at least a given capacity.
///
/// If the given capacity is 0, it defaults to 4.  The list will be empty (len=0).
/// @param [out] roots: root buffer to initialize.
/// @param [in] reserve: initial capacity, defaults to 4 if 0 is given
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Roots_init(nut_Roots *roots, uint64_t reserve);

/// Free resources held by a root buffer.
/// @param [in,out] roots: root buffer to destroy
[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
void nut_Roots_destroy(nut_Roots *roots);



/// Extend the capacity of a polynomial struct if needed.
///
/// Generally called at the beginning of most functions.
/// @param [in,out] f: polynomial to extend
/// @param [in] cap: minimum capacity the polynomial should have aferwards
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Poly_ensure_cap(nut_Poly *f, uint64_t cap);

/// Extend the capacity of a polynomial struct if needed and set any new terms to 0.
///
/// @param [in,out] f: polynomial to extend
/// @param [in] len: length to extend poly to
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Poly_zero_extend(nut_Poly *f, uint64_t len);

/// Extend the capacity of a root buffer if needed.
///
/// Called by { @link nut_Poly_roots_modn} and { @link nut_Poly_roots_modn_tmptmp}.
/// @param [in,out] roots: root buffer to extend
/// @param [in] cap: minimum capacity the root buffer should have aferwards
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Roots_ensure_cap(nut_Roots *roots, uint64_t cap);

/// Strip off terms so the leading term is nonzero if needed.
///
/// Generally called at the end of most functions.
/// @param [in,out] f: polynomial to normalize
[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
void nut_Poly_normalize(nut_Poly *f);

/// Reduce coefficients mod n, then normalize.
/// @param [in,out] f: polynomial to normalize
/// @param [in] n: modulus to reduce coefficients by
/// @param [in] use_negatives: if false, coefficients will be in the range [0,n), otherwise [(1-n)/2,n/2)
[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
void nut_Poly_normalize_modn(nut_Poly *f, int64_t n, bool use_negatives);

/// Reduce any exponents over cn
///
/// Generally called by functions that take a cn (carmichael lambda function of n) argument
/// Remove any term x**(k*cn + a) where 0 < a <= cn and add its coefficient to term x**a
[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
void nut_Poly_normalize_exps_modn(nut_Poly *f, uint64_t cn);



/// Copy a polynomial.
/// @param [out] g: destination
/// @param [in] f: source
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1, 2)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2)
bool nut_Poly_copy(nut_Poly *restrict g, const nut_Poly *restrict f);

/// Set a polynomial to a constant (linear terms and up zero).
/// @param [out] f: polynomial to set to constant
/// @param [in] c: constant to set polynomial to
/// @return 1 on success, 0 on failure
[[nodiscard, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Poly_setconst(nut_Poly *f, int64_t c);



/// Evaluate a polynomial mod a number.
///
/// Uses a Horner scheme with modulus taken after every multiply/add step.
/// @param [in] f: polynomial to evaluate
/// @param [in] x: point at which to evaluate
/// @param [in] n: modulus for the result
/// @return f(x) mod n
[[gnu::pure, gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_only, 1)
int64_t nut_Poly_eval_modn(const nut_Poly *f, int64_t x, int64_t n);

/// Print a polynomial to a file.
///
/// Skips zero terms, coefficients of 1 or -1, powers of 0 or 1, and any leading plus sign.
/// Thus the the format is "x^2 + 2x + 1".
/// The zero polynomial is displayed as "0", and a leading negative sign is always displayed as "-"
/// with no trailing space.
/// @param [in,out] file: file to write to
/// @param [in] f: polynomial to print
/// @param [in] vname: string to use for variable, eg "x"
/// @param [in] add: string to use for adding two monomials together, eg "+" or " + "
/// @param [in] sub: string to use for subtracting two monomials, eg "-" or " - "
/// @param [in] pow: string to use for exponents, eg "**" for "x**2" or "^" for "x^2"
/// @param [in] descending: if true, print higher terms down to lower terms, otherwise print lower terms up to higher terms (Taylor Series order)
/// @return the total number of characters printed
[[gnu::nonnull(1, 2, 3, 4, 5, 6)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3) NUT_ATTR_ACCESS(read_only, 4) NUT_ATTR_ACCESS(read_only, 5) NUT_ATTR_ACCESS(read_only, 6)
int nut_Poly_fprint(FILE *restrict file, const nut_Poly *restrict f, const char *vname, const char *add, const char *sub, const char *pow, bool descending);

/// Convert a string to a polynomial (and possibly modulus)
///
/// Ignores whitespace, ignores variable names, allows multiplication signs, duplicate terms, and duplicate minus signs.
/// The format is
///     poly = monomial (+|- monomial)* ("mod" \d+)?
///     monomial = -* ((coeff \*? vpow) | coeff | vpow)
///     coeff = \d+
///     vpow = \w (("**"|"^") \d+)?
/// @param [out] f: polynomial to store output in
/// @param [out] n: pointer to store modulus n in if found
/// @param [in] str: string to parse
/// @param [out] end: pointer to store end of parsed content (first unparsed character, ie typically '\0')
/// @return 0 on failure, 1 if polynomial was parsed without modulus, 2 if both polynomial and modulus were parsed
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(write_only, 2) NUT_ATTR_ACCESS(read_only, 3) NUT_ATTR_ACCESS(write_only, 4)
int nut_Poly_parse(nut_Poly *restrict f, int64_t *restrict n, const char *restrict str, const char *restrict *restrict end);

/// Generate a random polynomial.
///
/// Coefficients are integers chosen uniformly on [0,n).  max_len are chosen and then the polynomial is normalized ({ @link nut_Poly_normalize}).
/// So max_len minus the actual length has bounded geometric distribution with success probability 1/n.
/// @param [out] f: polynomial to randomize
/// @param [in] max_len: number of coefficients to generate before normalizing.  If 0, the output polynomial is always set to zero
/// @param [in] n: modulus for coefficients
/// @return 1 on succes, 0 on failure
[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_write, 1)
bool nut_Poly_rand_modn(nut_Poly *f, uint64_t max_len, int64_t n);



/// Add polynomials f + g mod n and store the result in h.
/// @param [out] h: polynomial in which to store the sum
/// @param [in] f, g: polynomials to add
/// @param [in] n: modulus to reduce result by
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Poly_add_modn(nut_Poly *h, const nut_Poly *f, const nut_Poly *g, int64_t n);

/// Subtract polynomials f - g mod n and store the result in h.
/// @param [out] h: polynomial in which to store the difference f - g mod n
/// @param [in] f, g: polynomials to subtract
/// @param [in] n: modulus to reduce result by
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Poly_sub_modn(nut_Poly *h, const nut_Poly *f, const nut_Poly *g, int64_t n);

/// Take the termwise product of polynomials f and g mod n and store the result in h.
/// @param [out] h: polynomial in which to store the result
/// @param [in] f, g: polynomials to multiply termwise
/// @param [in] n: modulus to reduce result by
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Poly_dot_modn(nut_Poly *h, const nut_Poly *f, const nut_Poly *g, int64_t n);

/// Multiply a polynomial f by a scalar a and store the result in h.
/// @param [out] g: polynomial in which to store a*f mod n
/// @param [in] f: polynomial to multiply by a scalar
/// @param [in] a: scalar to multiply by
/// @param [in] n: modulus to reduce result by
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2)
bool nut_Poly_scale_modn(nut_Poly *g, const nut_Poly *f, int64_t a, int64_t n);



/// Multiply polynomials f*g mod n and store the result in h.
///
/// Uses the basic algorithm
/// @param [out] h: polynomial in which to store f*g mod n
/// @param [in] f, g: polynomials to multiply
/// @param [in] n: modulus to reduce result by
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Poly_mul_modn(nut_Poly *restrict h, const nut_Poly *f, const nut_Poly *g, int64_t n);

/// Raise f**x mod n and store the result in g
///
/// @param [out] g: polynomial in which to store f**x
/// @param [in] f: polynomial to find a power of
/// @param [in] e: exponent to raise f to
/// @param [in] n: modulus to reduce result by
/// @param [in] cn: carmichael lambda function of n, basically modulus for exponents, see {@link carmichael_lambda} and {@link nut_sieve_carmichael}
/// @param [in] tmps: polynomial for scratch work, must have 1 initialized (or at least zeroed out) polynomial ({@link nut_Poly_pow_modn_tmptmp})
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 6)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_write, 6)
bool nut_Poly_pow_modn(nut_Poly *restrict g, const nut_Poly *restrict f, uint64_t e, int64_t n, uint64_t cn, nut_Poly tmps[restrict static 2]);

/// Raise f**x mod n and store the result in g
///
/// Like {@link nut_Poly_pow_modn} except temporary polynomials are allocated and freed internally instead of using
/// supplied temporaries.
/// @param [out] g: polynomial in which to store f**x
/// @param [in] f: polynomial to find a power of
/// @param [in] e: exponent to raise f to
/// @param [in] n: modulus to reduce result by
/// @param [in] cn: carmichael lambda function of n, basically modulus for exponents, see {@link carmichael_lambda} and {@link nut_sieve_carmichael}
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2)
bool nut_Poly_pow_modn_tmptmp(nut_Poly *restrict g, const nut_Poly *restrict f, uint64_t e, int64_t n, uint64_t cn);



/// Compose polynomial f(g(x)) mod n and store the result in h.
///
/// @param [out] h: polynomial in which to store f(g(x)) mod n
/// @param [in] f, g: polynomials to compose
/// @param [in] n: modulus to reduce result by
/// @param [in] cn: carmichael lambda function of n, basically modulus for exponents, see {@link carmichael_lambda} and {@link nut_sieve_carmichael}
/// @param [in] tmps: polynomial for scratch work, must have 1 initialized (or at least zeroed out) polynomial ({@link nut_Poly_compose_modn_tmptmp})
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 3, 6)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3) NUT_ATTR_ACCESS(read_write, 6)
bool nut_Poly_compose_modn(nut_Poly *restrict h, const nut_Poly *f, const nut_Poly *g, int64_t n, uint64_t cn, nut_Poly tmps[restrict static 2]);

/// Compose polynomial f(g(x)) mod n and store the result in h.
///
/// Like {@link nut_Poly_compose_modn} except temporary polynomials are allocated and freed internally instead of using
/// supplied temporaries.
/// @param [out] h: polynomial in which to store f(g(x)) mod n
/// @param [in] f, g: polynomials to compose
/// @param [in] n: modulus to reduce result by
/// @param [in] cn: carmichael lambda function of n, basically modulus for exponents, see {@link carmichael_lambda} and {@link nut_sieve_carmichael}
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Poly_compose_modn_tmptmp(nut_Poly *restrict h, const nut_Poly *f, const nut_Poly *g, int64_t n, uint64_t cn);



/// Divide two polynomials to get a quotient and remainder.
///
/// f = q*g + r mod n.  Uses extended synthetic division.
/// Note that while weaker forms of polynomial division can be defined for
/// weaker forms of rings than fields, the synthetic division algorithm can
/// fail if the coefficient ring is not a field because division by the leading
/// coefficient might fail.  Even if n is composite, if the leading coefficient
/// of the divisor is 1 then this algorithm will not fail, otherwise, this function
/// fails.
/// @param [out] q, r: polynomials in which to store the quotient and remainder such that f = q*g + r mod n
/// @param [in] f, g: polynomials to divide
/// @param [in] n: modulus to reduce result by
/// @return 1 on success, 0 on failure.  Can fail if the leading cofficient of g is not invertable mod n, as well as on allocation failure
[[gnu::nonnull(1, 2, 3, 4)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_write, 2) NUT_ATTR_ACCESS(read_only, 3) NUT_ATTR_ACCESS(read_only, 4)
bool nut_Poly_quotrem_modn(nut_Poly *restrict q, nut_Poly *restrict r, const nut_Poly *restrict f, const nut_Poly *restrict g, int64_t n);



//TODO: use arrays for temporaries instead of single pointers, maybe provide a nut_Poly[1] typedef as is tradition
/// Find the gcd of two polynomials.
///
/// The gcd is reduced to be monic, so if n is not prime this can fail.
/// If d = gcd(f, g), then there is some pair of polynomials a and b so
/// that f = ad and g = bd (d is a common divisor), and for any other
/// common divisor h, d = ch for some polynomial c (d is greatest because
/// it has at least as high degree as any other common divisor.
/// Note that in the integers, we have two choices for gcd, which are associates
/// (meaning one is the negative of the other) and we pick the positive one
/// by convention.  Similarly for polynomials with coefficients from a field,
/// the the gcd is not unique and any associate of the gcd is also a gcd.
/// Thus we restrict the gcd to be monic.
/// @param [out] d: polynomial in which to store the monic gcd, or store 0 if both f and n are 0
/// @param [in] f, g: polynomials to take the gcd of
/// @param [in] n: modulus for the coefficient ring
/// @param [in] tmps: polynomials for scratch work, cannot be NULL ({ @link nut_Poly_gcd_modn_tmptmp})
/// @return 1 on success, 0 on failure.  Can fail if a non-invertable leading coefficient is encountered, as well as on allocation failure.
[[gnu::nonnull(1, 2, 3, 5)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3) NUT_ATTR_ACCESS(read_write, 5)
bool nut_Poly_gcd_modn(nut_Poly *restrict d, const nut_Poly *restrict f, const nut_Poly *restrict g, int64_t n, nut_Poly tmps[restrict static 3]);

/// Find the gcd of two polynomials.
///
/// Like {@link nut_Poly_gcd_modn} except temporary polynomials are allocated and freed internally instead of using
/// supplied temporaries.  Less efficient because early checking is not done.
/// @param [out] d: polynomial in which to store the monic gcd, or store 0 if both f and n are 0
/// @param [in] f, g: polynomials to take the gcd of
/// @param [in] n: modulus for the coefficient ring
/// @return 1 on success, 0 on failure.  Can fail if a non-invertable leading coefficient is encountered, as well as on allocation failure.
[[gnu::nonnull(1, 2, 3)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 3)
bool nut_Poly_gcd_modn_tmptmp(nut_Poly *restrict d, const nut_Poly *restrict f, const nut_Poly *restrict g, int64_t n);



/// Compute the power of a polynomial mod another polynomial.
///
/// Uses binary exponentiation.  The coefficients themselves are in the cyclic group Z_n as usual.
/// @param [out] h: polynomial in which to store the result
/// @param [in] f: base
/// @param [in] e: exponent
/// @param [in] g: modulus
/// @param [in] n: coefficient modulus
/// @param [in] tmps: polynomials for scratch work, must have 3 initialized (or at least zeroed out) polynomials ({@link nut_Poly_powmod_modn_tmptmp})
/// @return 1 on success, 0 on failure.  Can fail if a non-invertable leading coefficient is encountered, as well as on allocation failure.
[[gnu::nonnull(1, 2, 4, 6)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 4) NUT_ATTR_ACCESS(read_write, 6)
bool nut_Poly_powmod_modn(nut_Poly *restrict h, const nut_Poly *restrict f, uint64_t e, const nut_Poly *restrict g, int64_t n, nut_Poly tmps[restrict static 3]);

/// Compute the power of a polynomial mod another polynomial.
///
/// Like {@link nut_Poly_powmod_modn} except temporary polynomials are allocated and freed internally instead of using
/// supplied temporaries.  Less efficient because early checking is not done.
/// @param [out] h: polynomial in which to store the result
/// @param [in] f: base
/// @param [in] e: exponent
/// @param [in] g: modulus
/// @param [in] n: coefficient modulus
/// @return 1 on success, 0 on failure.  Can fail if a non-invertable leading coefficient is encountered, as well as on allocation failure.
[[gnu::nonnull(1, 2, 4)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_only, 4)
bool nut_Poly_powmod_modn_tmptmp(nut_Poly *restrict h, const nut_Poly *restrict f, uint64_t e, const nut_Poly *restrict g, int64_t n);

/// Compute the product of all irreducible factors of degree d.
///
/// This should be called repeatedly on squarefree factors, or else it will extract repeat factors for d > 1.
/// Thus if we only want roots in the field Z_p, or equivalently we only want linear factors, we can call
/// this function without bothering to make f squarefree first.
/// In the general case, getting the squarefree factorization of f is easy: any repeated root of f (in Z_p or indeed the splitting field of f)
/// is also a root of f' the derivative of f, so we can say the squarefree part of f is f / gcd(f, f').
/// We can repeat this to find polynomials representing roots which occur exactly once (f / gcd(f, f')), roots which occur exactly twice (sqrt(gcd(f, f') / gcd(f, f''))),
/// roots which occur exactly thrice (cbrt(gcd(f, f'') / gcd(f, f'''))), etc.  Eventually, these quotients will saturate and just become 1 forever.
/// If f is already squarefree, let f_i be the product of all irreducible factors of f with degree i.
/// Then we have f_1 = gcd(x^p-x, f), f_2 = gcd(x^(p^2), f/f_1), f_3 = gcd(x^(p^3), f/f_1/f_2), etc.
/// Once the expression f/f_1/f_2/.../f_i becomes 1, we are done.  This clearly takes at most k steps,
/// where k is the degree of f.  It is possible to show an irreducible f is irreducible by only showing
/// it has no factors of degree k/q for any prime divisor q, ie f | x^(p^k) but gcd(f, x^(p^(k/q))) = 1 for all q.
/// This function only takes care of computing gcd(x^(p^d)-x, f/f_1/.../f_(d-1)) given p, d, and f/f_1/.../f_(d-1).
/// @param [out] f_d: the product of all degree d irreducible factors of f
/// @param [in] f: the polynomial to factor, which should be squarefree and have no irreducible factors with degree less than d
/// @param [in] d: the degree of irreducible factors to extract.  Note that f should be free from irreducible factors with lower degree,
/// so once d exceeds half the degree of f, f must be irreducible.
/// @param [in] n: coefficient modulus.  MUST be prime.
/// @param [in] tmps: polynomials for scratch work, must have 4 initialized (or at least zeroed out) polynomials
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 5)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_write, 5)
bool nut_Poly_factors_d_modn(nut_Poly *restrict f_d, const nut_Poly *restrict f, uint64_t d, int64_t n, nut_Poly tmps[restrict static 4]);

/// Attempt to find a divisor of f_d.
///
/// Given a squarefree polynomial f whose irreducible factors all have degree d,
/// if f has degree k which is greater than d (it must be a multiple of d) it must split.
/// Mod an odd prime p, we can find a nontrivial factor of f by 1) choosing a random r with
/// degree less than k (if gcd(r, f) isn't 1 we are done), 2) computing s = r^((p^d-1)/2) mod f,
/// 3) gcd(s + 1, f) is a nontrivial factor with probability greater than 1/2.
/// This works because not only does x^(p^d)-x factor as x(x^((p^d-1)/2)-1)(x^((p^d-1)/2)+1),
/// this holds if we replace x with any polynomial.  If a factor isn't contained in the first term, it has
/// about a 1/2 chance of being contained in the second term, hence this algorithm.
/// This function returns a single nontrivial factor g, and both g and f/g may need to be factored further
/// if they are not degree d.  Note that if f is degree 2 with irreducible factors of degree 1, completing the
/// square/using the quadratic formula is faster, although for such small polynomials this is unimportant.
/// @param [out] g: a nontrivial factor of f.  Both g and f/g could potentially need to be factored further with more calls to this method.
/// @param [in] f: the polynomial to factor, which should be squarefree and have multiple irreducible factors of degree d
/// @param [in] d: the degree of irreducible factors of f
/// @param [in] n: coefficient modulus.  MUST be prime.
/// @param [in] tmps: polynomials for scratch work, must have 4 initialized (or at least zeroed out) polynomials
/// @return 1 on success, 0 on failure
[[gnu::nonnull(1, 2, 5)]]
NUT_ATTR_ACCESS(read_write, 1) NUT_ATTR_ACCESS(read_only, 2) NUT_ATTR_ACCESS(read_write, 5)
bool nut_Poly_factor1_modn(nut_Poly *restrict g, const nut_Poly *restrict f, uint64_t d, int64_t n, nut_Poly tmps[restrict static 4]);

/// Find all roots of a polynomial in an odd prime field.
///
/// Uses the Cantor-Zassenhaus algorithm, specialized for linear factors (so we can skip squarefree factorization).
/// This is implemented in the functions {@link nut_Poly_factors_d_modn} and {@link nut_Poly_factor1_modn}.
/// After the first function, we know how many roots there are and allocate a buffer.  The second function is called
/// as many times as needed to find all the roots.  TODO: optimize factoring quadratics using quadratic formula.
/// @param [in] f: polynomial to find roots of
/// @param [in] n: modulus of the coefficient field.  MUST be prime.
/// @param [out] roots: store roots here.
/// @param [in] tmps: polynomials for scratch work, must have 6 initialized (or at least zeroed out) polynomials ({@link nut_Poly_roots_modn_tmptmp})
/// @return 1 on success (including if the polynomial has 0 roots), 0 on failure.
[[gnu::nonnull(1, 3, 4)]]
NUT_ATTR_ACCESS(read_only, 1) NUT_ATTR_ACCESS(read_write, 3) NUT_ATTR_ACCESS(read_write, 4)
bool nut_Poly_roots_modn(const nut_Poly *restrict f, int64_t n, nut_Roots *restrict roots, nut_Poly tmps[restrict static 6]);

/// Find all roots of a polynomial in an odd prime field.
///
/// Like {@link nut_Poly_roots_modn} except temporary polynomials are allocated and freed internally instead of using
/// supplied temporaries.
/// @param [in] f: polynomial to find roots of
/// @param [in] n: modulus of the coefficient field.  MUST be prime.
/// @param [out] roots: store roots here.
/// @return 1 on success (including if the polynomial has 0 roots), 0 on failure.
[[gnu::nonnull(1, 3)]]
NUT_ATTR_ACCESS(read_only, 1) NUT_ATTR_ACCESS(read_write, 3)
bool nut_Poly_roots_modn_tmptmp(const nut_Poly *restrict f, int64_t n, nut_Roots *restrict roots);


