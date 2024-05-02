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

#include <nut/debug.h>

typedef struct{
	int64_t rows, cols;
	int64_t *buf;
} nut_i64_Matrix;

/// Allocate internal buffers for a matrix
/// @param [in] rows, cols: The number of rows and columns.  Must not be negative.
/// If either is zero, most operations will become no-ops.
/// Many operations are specific to square matrices
/// @param [out] self: matrix to initialize
/// @return true on success, false on failure (allocation failure, size overflow, negative size)
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(write_only, 1)
bool nut_i64_Matrix_init(nut_i64_Matrix *self, int64_t rows, int64_t cols);

/// Copy the values from one matrix to another, which must be initialized
NUT_ATTR_NONNULL(1, 2)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_only, 2)
void nut_i64_Matrix_copy(nut_i64_Matrix *restrict dest, const nut_i64_Matrix *restrict src);

/// Deallocate internal buffers for a matrix table
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
void nut_i64_Matrix_destroy(nut_i64_Matrix *self);

/// Print a matrix to a given file.
/// Format is '[[0, 1, 2],\n[3, 4, 5],\n[6, 7, 8]]
NUT_ATTR_NONNULL(1, 2)
NUT_ATTR_ACCESS(read_only, 1)
NUT_ATTR_ACCESS(read_write, 2)
void nut_i64_Matrix_fprint(const nut_i64_Matrix *restrict self, FILE *file);

/// Fill a matrix with zeros off the diagonal and ones on the diagonal
/// The matrix must be squared
NUT_ATTR_NONNULL(1)
NUT_ATTR_ACCESS(read_write, 1)
bool nut_i64_Matrix_fill_I(nut_i64_Matrix *self);

/// Multiply a matrix by a vector, returning a vector
/// The vectors are stored as one dimensional arrays
/// This is a special case of matrix-matrix multiplication, since a row/column
/// vec is simply a matrix with 1 row/column respectively.
/// The advantage is that this can interprest more things (any buffer) as a vector.
/// @param [in] self: The matrix to multiply the input vector by
/// @param [in] vec: The vector to multiply.  The number of cols in the matrix must match the length of this vector.
/// If at least self->cols elements of vec are not accessible, it is undefined behavior.
/// @param [out] out: The vector to store the result in.  The number of rows in the matrix must match the length of this vector.
/// If at least self->rows elements of out are not accessible, it is undefined behavior
NUT_ATTR_NONNULL(1, 2, 3)
NUT_ATTR_ACCESS(read_only, 1)
NUT_ATTR_ACCESS(read_only, 2)
NUT_ATTR_ACCESS(write_only, 3)
void nut_i64_Matrix_mul_vec(const nut_i64_Matrix *restrict self, const int64_t vec[restrict static self->cols], int64_t out[restrict static self->rows]);

/// Fill a matrix with part of the Pascal Triangle
/// This is mostly internally used together with { @link nut_i64_Matrix_invert_ltr } to get a matrix of Faulhaber polynomial
/// coefficients to quickly sum powers of integer (see { @link nut_Diri_compute_Nk }).
/// In particular, the ith row of the matrix will be filled by the ith row of pascal's triangle without the final 1,
/// trimmed/zero padded as needed.  The first (0th) row is a 1 followed by all zeroes, the second (1st) row is a 1,
/// then a 2, then all zeroes, etc.
void nut_i64_Matrix_fill_short_pascal(nut_i64_Matrix *self);

/// Invert a lower triangular matrix
/// Integer matrices in general are not invertable because the integers are not a field, so this function also returns a denominator so
/// the result can be interpreted as a rational matrix.  This denominator may not be reduced in general.
/// Uses gaussian elimination specialized to lower triangular matrices.
/// @param [in] self: The matrix to invert, which is clobbered in place.  If the original matrix will still be needed, copy it beforehand
/// @param [out] tmp: The matrix to store the result in
/// @return 0 if self/tmp have different sizes or are not square or self is not invertable, otherwise the denominator (all entries should be understood as numeratos)
NUT_ATTR_NONNULL(1, 2)
NUT_ATTR_ACCESS(read_write, 1)
NUT_ATTR_ACCESS(read_write, 2)
int64_t nut_i64_Matrix_invert_ltr(nut_i64_Matrix *restrict self, nut_i64_Matrix *restrict tmp);

/// Fill a vector with increasing powers of x
/// @param [in] x: The number to take powers of
/// @param [in] k: The max exponent (inclusive)
/// @param [in] m: The modulus to reduce powers by, or 0 to skip reducing
/// @param [out] out: Buffer to store resulting vector in.  The first entry will be set to x%m, then x*x%m then x*x*x%m, etc
void nut_i64_Matrix_fill_vandemond_vec(uint64_t x, uint64_t k, uint64_t m, int64_t out[static k + 1]);

