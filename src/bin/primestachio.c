#define _POSIX_C_SOURCE 202208L
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>

#include "nut/sieves.h"

int main(int argc, char **argv){
	if(argc != 2){
		fprintf(stderr, "\e[1;31mNo upper bound specified.  Please use like\e[0m\n\e[1;31m%s <MAX>\e[0m\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char *str_end = NULL;
	uint64_t sieve_max = strtoull(argv[1], &str_end, 10);
	if(!str_end || str_end == argv[1]){
		fprintf(stderr, "\e[1;31mCould not parse upper bound\e[0m\n");
		exit(EXIT_FAILURE);
	}
	fprintf(stderr, "\e[1;34mSieving primes up to %"PRIu64"...\e[0m\n", sieve_max);
	struct timespec start_time, end_time;
	uint64_t num_primes;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	uint64_t *primes = sieve_primes(sieve_max, &num_primes);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
	if(!primes){
		fprintf(stderr, "\e[1;31mCould not allocate memory!\e[0m\n");
		exit(EXIT_FAILURE);
	}
	double tot_time = (end_time.tv_sec + end_time.tv_nsec*1e-9 - start_time.tv_sec - start_time.tv_nsec*1e-9);
	printf("Found %"PRId64" primes in %.9fs\n", num_primes, tot_time);
	free(primes);
}

