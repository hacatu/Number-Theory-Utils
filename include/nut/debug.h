#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Functions for internal use in the tests

#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>

#include <nut/modular_math.h>

[[gnu::artificial]]
NUT_ATTR_ACCESS(read_write, 1)
static inline void cleanup_free(void *_p){
	free(*(void**)_p);
	*(void**)_p = NULL;
}

[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_only, 1) NUT_ATTR_ACCESS(none, 2)
static inline void check_alloc(const char *restrict what, const void *restrict buf){
	if(!buf){
		fprintf(stderr, "\e[1;31mAllocation failed for %s!\e[0m\n", what);
		exit(0);
	}
}

[[gnu::nonnull(1)]]
NUT_ATTR_ACCESS(read_only, 1)
static inline void print_summary(const char *what, uint64_t correct, uint64_t total){
	fprintf(stderr, "%s (found %s for %"PRIu64"/%"PRIu64" numbers correctly)\e[0m\n", correct == total ? "\e[1;32mPASSED" : "\e[1;31mFAILED", what, correct, total);
}

