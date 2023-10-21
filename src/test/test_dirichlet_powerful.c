#define _POSIX_C_SOURCE 202305L
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <nut/debug.h>
#include <nut/modular_math.h>
#include <nut/dirichlet.h>
#include <nut/dirichlet_powerful.h>
#include <nut/sieves.h>
#include <nut/factorization.h>
#include <nut/polynomial.h>

static int64_t h_powerfulpart(uint64_t p, uint64_t pp, uint64_t e, int64_t m){
	uint64_t res;
	if(e < 2){
		return 1 - e;
	}else if(e == 2){
		res = pp - 1;
	}else{
		res = pp - pp/p;
	}
	return m ? res%m : res;
}

static void test_adjust_powerfulpart(uint64_t sieve_max){
	fprintf(stderr, "\e[1;34mTesting powerful part sums using powerful number adjustments...\e[0m\n");
	nut_Diri u_table = {};
	nut_Diri_init(&u_table, sieve_max, 0);
	check_alloc("u table", u_table.buf);
	nut_Diri_compute_u(&u_table, 0);
	int64_t res;
	check_alloc("pp iterator primes", (void*)nut_Diri_sum_adjusted_hpe(&res, 0, &u_table, h_powerfulpart));
	uint64_t w = nut_max_prime_divs(sieve_max);
	uint64_t *largest_factors = nut_sieve_largest_factors(sieve_max);
	nut_Factors *factors = nut_make_Factors_w(w);
	int64_t correct = 3;
	for(uint64_t n = 4; n <= sieve_max; ++n){
		nut_fill_factors_from_largest(factors, n, largest_factors);
		uint64_t pfp = 1;
		for(uint64_t i = 0; i < factors->num_primes; ++i){
			uint64_t p = factors->factors[i].prime;
			uint64_t e = factors->factors[i].power;
			if(e >= 2){
				pfp *= nut_u64_pow(p, e);
			}
		}
		correct += pfp;
	}
	nut_Diri_destroy(&u_table);
	free(largest_factors);
	free(factors);
	if(correct == res){
		fprintf(stderr, "\e[1;32mGot sum of powerful parts of numbers up to %"PRIu64" = %"PRIi64"\e[0m\n", sieve_max, res);
	}else{
		fprintf(stderr, "\e[1;31mGot sum of powerful parts of numbers up to %"PRIu64" = %"PRIi64"; expected %"PRIi64"\e[0m\n", sieve_max, res, correct);
	}
}

static void test_countprimes(uint64_t sieve_max){
	fprintf(stderr, "\e[1;34mTesting prime counting via powerful number trick\e[0m\n");
	int64_t x = nut_u64_next_prime_ge(nut_max_primes_le(sieve_max));
	fprintf(stderr, "\e[1;34mChose %"PRIi64" as representative for x\e[0m\n", x);
	uint64_t bits = 64 - __builtin_clzll(sieve_max);
	uint64_t num_primes = 0, num_semilike = 0, num_spheniclike = 0;
	uint64_t *sigma_0_vals = nut_sieve_sigma_0(sieve_max);
	for(uint64_t n = 2; n <= sieve_max; ++n){
		switch(sigma_0_vals[n]){
			case 2: ++num_primes; break;
			case 4: ++num_semilike; break;
			case 8: ++num_spheniclike;
		}
	}
	free(sigma_0_vals);
	fprintf(stderr, "\e[1;34mCorrect counts are:\n#{n:d(n)=2}: %"PRIu64"\n#{n:d(n)=4}: %"PRIu64"\n#{n:d(n)=8}: %"PRIu64"\e[0m\n", num_primes, num_semilike, num_spheniclike);
	int64_t *h_vals = malloc((bits + 1)*sizeof(uint64_t));
	check_alloc("perturbative bell series", h_vals);
	h_vals[0] = 1;
	h_vals[1] = 0;
	for(uint64_t e = 2; e <= bits; ++e){
		h_vals[e] = (x - nut_i64_modinv(e, x))*x;
	}
	nut_Diri dx_tbl, f_tbl, g_tbl;
	nut_Diri_init(&dx_tbl, sieve_max, 0);
	check_alloc("dx table", dx_tbl.buf);
	nut_Diri_init(&f_tbl, sieve_max, 0);
	check_alloc("f table", f_tbl.buf);
	nut_Diri_init(&g_tbl, sieve_max, 0);
	check_alloc("g table", g_tbl.buf);
	nut_Diri_compute_dk(&dx_tbl, x, x*x, &f_tbl, &g_tbl);
	int64_t res;
	nut_Diri_sum_adjusted_he(&res, x*x, &dx_tbl, h_vals);
	if(res == (int64_t)num_primes*x + 1){
		fprintf(stderr, "\e[1;32mGot sum of P_1(n)(x) = %"PRIi64"x + 1\e[0m\n", res/x);
	}else{
		fprintf(stderr, "\e[1;31mGot sum of P_1(n)(x) = %"PRIi64"x + %"PRIi64", expected %"PRIu64"x + 1\e[0m\n", res/x, res%x, num_primes);
	}
	nut_Diri_destroy(&dx_tbl);
	nut_Diri_destroy(&f_tbl);
	nut_Diri_destroy(&g_tbl);
	free(h_vals);
}

int main(int argc, char **argv){
	if(argc != 2){
		fprintf(stderr, "\e[1;31mIncorrect arguments, please specify the bound to sieve up to\e[0m\n");
		exit(EXIT_FAILURE);
	}
	char *end = NULL;
	uint64_t sieve_max = strtoull(argv[1], &end, 10);
	if(!end || *end){
		fprintf(stderr, "\e[1;31mCould not parse sieve upper bound\e[0m\n");
		exit(EXIT_FAILURE);
	}
	test_adjust_powerfulpart(sieve_max);
	test_countprimes(sieve_max);
}

