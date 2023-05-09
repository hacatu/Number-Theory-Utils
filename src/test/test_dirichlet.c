#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nut/dirichlet.h>
#include <nut/sieves.h>
#include <nut/factorization.h>

static const uint64_t sieve_max = 1000;
static uint64_t *sigma_vals;
static int64_t *f_vals;
static int64_t *h_vals;

static void check_alloc(const char *what, const void *buf){
	if(!buf){
		fprintf(stderr, "\e[1;31mAllocation failed for %s!\e[0m\n", what);
		exit(0);
	}
}

static void print_summary(const char *what, uint64_t correct){
	fprintf(stderr, "%s (found %s for %"PRIu64"/%"PRIu64" numbers correctly)\e[0m\n", correct == sieve_max ? "\e[1;32mPASSED" : "\e[1;31mFAILED", what, correct, sieve_max);
}

static void test_dirichlet_D(){
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
}

static void test_euler_sieve_conv_u(){
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
}

static void test_compute_conv_u_diri(){
	diri_table dkp_table = {}, dk_table = {};
	diri_table_init(&dkp_table, sieve_max, 0);
	diri_table_init(&dk_table, sieve_max, 0);
	check_alloc("u table", dkp_table.buf);
	check_alloc("dk table", dk_table.buf);
	compute_u_diri_table(&dkp_table);
	for(uint64_t k = 2; k < 7; ++k){
		fprintf(stderr, "\e[1;34mFinding d_%"PRIu64"...\e[0m\n", k);
		compute_conv_u_diri_table(&dk_table, &dkp_table);
		uint64_t *dk_vals = sieve_dk(sieve_max, k);
		check_alloc("dk sieve", dk_vals);
		for(int64_t i = 1; i <= dk_table.y; ++i){
			if(dk_vals[i] != (uint64_t)diri_table_get_dense(&dk_table, i)){
				fprintf(stderr, "\e[1;31mu <*> u table was wrong at dense %"PRIi64"\e[0m\n", i);
			}
		}
		for(uint64_t i = 2; i <= sieve_max; ++i){
			dk_vals[i] += dk_vals[i - 1];
		}
		for(int64_t i = 1; i < dk_table.yinv; ++i){
			if(dk_vals[dk_table.x/i] != (uint64_t)diri_table_get_sparse(&dk_table, i)){
				fprintf(stderr, "\e[1;31mu <*> u table was wrong at sparse %"PRIi64"\e[0m\n", i);
			}
		}
		free(dk_vals);
		int64_t *tmp = dkp_table.buf;
		dkp_table.buf = dk_table.buf;
		dk_table.buf = tmp;
	}
	diri_table_destroy(&dkp_table);
	diri_table_destroy(&dk_table);
}

int main(){
	fprintf(stderr, "\e[1;34mTesting dirichlet functions...\e[0m\n");
	sigma_vals = sieve_sigma_0(sieve_max);
	f_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	h_vals = malloc((sieve_max + 1)*sizeof(int64_t));
	check_alloc("Sigma vals", sigma_vals);
	check_alloc("f vals", f_vals);
	check_alloc("h vals", h_vals);
	test_dirichlet_D();
	test_euler_sieve_conv_u();
	test_compute_conv_u_diri();
	free(f_vals);
	free(h_vals);
}

