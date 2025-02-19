#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Sieve based functions that work on segments, so that you can use all
/// the cores your little heart desires

#include <inttypes.h>
#include <stdlib.h>

#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/sieves.h>

typedef struct{
	uint64_t max;
	uint64_t sqrt_max;
	uint64_t num_primes;
	uint64_t *primes;
	uint64_t preferred_bucket_size;
} nut_Segsieve;

/// Set up the sieving primes header for a segmented sieve.
/// The user can then divide the range into segments as desired and process them in multiple threads,
/// by calling `nut_Segsieve_*` functions with `[a, b)` intervals that partition the range.
/// @param [out] self: the segsieve header to initialize.  Must be freed with { @link nut_Segsieve_destroy }
/// @param [in] max: the inclusive upper bound of the range.
/// @param [in] preferred_bucket_size: designed to tell algorithms how big buckets should be for optimal cache use etc.
/// Currently not really used.  If 0, we just pick sqrt_max, and if nonzero, we forward it to self
/// @return true on success, in which case self contains a list of sieving primes and other information, or false on
/// (allocation) failure
bool nut_Segsieve_init(nut_Segsieve *self, uint64_t max, uint64_t preferred_bucket_size);

/// Free the resources held by a segmented sieve header.
/// Note that this does not free per-thread work buffers or other resources not directly managed by self.
void nut_Segsieve_destroy(nut_Segsieve *self);

/// Sieve all factorizations in the range [a, b) using a modified, in place largest factor sieve
/// This uses a pitched array of factorization structs { @link nut_Factors }, so it uses a lot of memory per element
/// in the segment, and the segment size should be lowered accordingly.
/// In particular, we use 8 + 16*(omega+1) bytes PER segment entry, where omega is the max number of distinct prime divisors,
/// which can be up to 15, ie up to 264 bytes per entry.  So for example, if using a CPU with 512k L1d cache per core and smt,
/// and sieving numbers over 614889782588491410, omega can be 15 so we would use a preferred bucket size of 992.
/// If the upper bound were 10^12 instead, then omega is at most 11, so we use a preferred bucket size of 1310,
/// Or you can just try a larger bucket size, since ~1k buckets are insanely small.  My computer also has 16m L2 cache per core,
/// which would lead to a bucket size of 41943 for omega 11 or 31775 for omega 15.
///
/// @param [out] buffer: the factorizations are stored here, but for performance reasons, their first prime power
/// will have the form p^1, where p is either the LARGEST prime divisor OR 1, so code consuming this data must be aware of that.
void nut_Segsieve_factorizations(const nut_Segsieve *restrict self, uint64_t a, uint64_t b, size_t pitch, void *buffer);

/// Allocate a buffer for a thread to pass to { @link nut_Segsieve_factorizations } on its intervals.
/// @param [in,out] pitch: if 0, compute max omega for self->max and compute the related pitch, then store it here.
/// if nonzero, we assume you are passing the pitch calculated from a previous call.
/// This pointer must not be null.
void *nut_Segsieve_factorizations_mkbuffer(const nut_Segsieve *self, size_t *pitch);