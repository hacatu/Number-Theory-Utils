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

// oeis.org/A084237
static const int64_t M10[] = {1, -1, 1, 2, -23, -48, 212, 1037, 1928, -222, -33722, -87856, 62366, 599582, -875575, -3216373, -3195437, -21830254, -46758740, 899990187, 461113106, -3395895277, -2061910120, 62467771689};

static void test_mertens(uint64_t t){
	uint64_t x = pow_u64(10, t);
	uint64_t y = 0.25*pow(x, 2./3);
	struct timespec start, sieve_done, end;
	diri_table mertens_table = {};
	diri_table_init(&mertens_table, x, y);
	check_alloc("Mertens table", mertens_table.buf);
	y = mertens_table.y;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	uint8_t *mobius = sieve_mobius(y);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &sieve_done);
	check_alloc("Mobius sieve", mobius);
	compute_mertens_diri_table(&mertens_table, mobius);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
	double sieve_secs = sieve_done.tv_sec - start.tv_sec + (sieve_done.tv_nsec - start.tv_nsec)*1e-9;
	double diri_secs = end.tv_sec - sieve_done.tv_sec + (end.tv_nsec - sieve_done.tv_nsec)*1e-9;
	int64_t M10_12 = diri_table_get_sparse(&mertens_table, 1);
	diri_table_destroy(&mertens_table);
	free(mobius);
	fprintf(stderr, "%s: sieved up to %"PRIu64" in %.3fs; got Mertens(10^%"PRIu64") = %"PRIi64" in %.3fs\e[0m\n", M10_12 == M10[t] ? "\e[1;32mSUCCESS" : "\e[1;31mERROR", y, sieve_secs, t, M10_12, diri_secs);
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
	test_mertens(12);
	free(f_vals);
	free(h_vals);
}

