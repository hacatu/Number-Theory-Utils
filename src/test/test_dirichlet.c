#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nut/dirichlet.h>
#include <nut/sieves.h>

const uint64_t sieve_max = 1000;

static void check_alloc(const char *what, const void *buf){
	if(!buf){
		fprintf(stderr, "\e[1;31mAllocation failed for %s!\e[0m\n", what);
		exit(0);
	}
}

static void print_summary(const char *what, uint64_t correct){
	fprintf(stderr, "%s (found %s for %"PRIu64"/%"PRIu64" numbers correctly)\e[0m\n", correct == sieve_max ? "\e[1;32mPASSED" : "\e[1;31mFAILED", what, correct, sieve_max);
}

int main(){
	uint64_t *sigma_vals = sieve_sigma_0(sieve_max);
	int64_t *f_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	int64_t *h_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	check_alloc("Sigma vals", sigma_vals);
	check_alloc("f vals", f_vals);
	check_alloc("h vals", h_vals);
	uint64_t correct = 0;
	uint64_t acc = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		acc += sigma_vals[n];
		if(dirichlet_D(n) == acc){
			++correct;
		}else{
			fprintf(stderr, "\e[1;31mdirichlet_D(%"PRIu64") should be %"PRIu64"\e[0m\n", n, acc);
		}
		f_vals[n] = 1;
	}
	print_summary("dirichlet_D", correct);
	if(!euler_sieve_conv_u(sieve_max, f_vals, h_vals)){
		fprintf(stderr, "\e[1;31mAllocation failed for linear_sieve_conv_u!\e[0m\n");
	}
	if(memcmp(sigma_vals + 1, h_vals + 1, sieve_max*sizeof(int64_t))){
		fprintf(stderr, "\e[1;31mlinear sieve failed to compute u <*> u!\e[0m\n");
	}
	free(sigma_vals);
	memset(f_vals + 2, 0, (sieve_max - 2)*sizeof(int64_t));
	if(!euler_sieve_conv_u(sieve_max, f_vals, h_vals)){
		fprintf(stderr, "\e[1;31mAllocation failed for linear_sieve_conv_u!\e[0m\n");
	}
	for(uint64_t n = 1; n <= sieve_max; ++n){
		if(h_vals[n] != 1){
			fprintf(stderr, "\e[1;31mlinear sieve failed to compute I <*> u!\e[0m\n");
			break;
		}
	}
	free(f_vals);
	free(h_vals);
}

