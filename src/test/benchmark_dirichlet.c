#define _POSIX_C_SOURCE 202305L
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <nut/modular_math.h>
#include <nut/dirichlet.h>
#include <nut/sieves.h>
#include <nut/factorization.h>
#include <nut/debug.h>

static const uint64_t sieve_max = 1'000'000;

int main(){
	int64_t *e_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	int64_t *o_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	int64_t *d2_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	int64_t *d3_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	check_alloc("d5 vals (euler sieve)", e_vals);
	check_alloc("d5 vals (grouped sieve)", o_vals);
	check_alloc("d2 vals", d2_vals);
	check_alloc("d3 vals", d3_vals);
	e_vals[0] = 0;
	for(uint64_t i = 1; i <= sieve_max; ++i){
		e_vals[i] = 1;
	}
	nut_euler_sieve_conv_u(sieve_max, e_vals, d2_vals);
	nut_euler_sieve_conv_u(sieve_max, d2_vals, d3_vals);
	struct timespec time_0, time_1, time_2;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_0);
	nut_euler_sieve_conv(sieve_max, d2_vals, d3_vals, e_vals);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_1);
	nut_grouped_sieve_conv(sieve_max, d2_vals, d3_vals, o_vals);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_2);
	double euler_sieve_secs = time_1.tv_sec - time_0.tv_sec + (time_1.tv_nsec - time_0.tv_nsec)*1e-9;
	double grouped_sieve_secs = time_2.tv_sec - time_1.tv_sec + (time_2.tv_nsec - time_1.tv_nsec)*1e-9;
	fprintf(stderr, "Euler sieve ran in %fs\n", euler_sieve_secs);
	fprintf(stderr, "Grouped sieve ran in %fs\n", grouped_sieve_secs);
	fprintf(stderr, memcmp(e_vals + 1, o_vals + 1, sieve_max*sizeof(int64_t)) ? "\e[1;31mOutputs do not match!\e[0m\n" : "\e[1;32mOutputs match!\e[0m\n");
	free(e_vals);
	free(o_vals);
	free(d2_vals);
	free(d3_vals);
}

