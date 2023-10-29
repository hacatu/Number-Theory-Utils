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

static int64_t h_powerfulpart(uint64_t p, uint64_t pp, uint64_t e, uint64_t m){
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
	nut_PfIt pf_it = {};
	check_alloc("powerful iterator", (void*)NUT_PfIt_INIT(&pf_it, sieve_max, 0, h_powerfulpart));
	nut_Diri_init(&u_table, sieve_max, 0);
	check_alloc("u table", u_table.buf);
	nut_Diri_compute_u(&u_table, 0);
	int64_t res;
	check_alloc("pp iterator primes", (void*)nut_Diri_sum_adjusted(&res, &u_table, &pf_it));
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

static void verify_dk(const nut_Diri *dx_tbl, int64_t x, int64_t modulus){
	uint64_t *dx_vals [[gnu::cleanup(cleanup_free)]] = nut_sieve_dk(dx_tbl->x, x, modulus);
	int64_t Dn = 0;
	for(int64_t n = 1; n <= dx_tbl->x; ++n){
		int64_t dn = dx_vals[n];
		Dn = nut_i64_mod(Dn + dn, modulus);
		if(n <= dx_tbl->y){
			if(dn != nut_Diri_get_dense(dx_tbl, n)){
				fprintf(stderr, "\e[1;31mdx_tbl.dense[%"PRIi64"] is wrong\e[0m\n", n);
			}
		}else if((dx_tbl->x/(n + 1) != dx_tbl->x/n) || n == dx_tbl->x){
			if(Dn != nut_Diri_get_sparse(dx_tbl, dx_tbl->x/n)){
				fprintf(stderr, "\e[1;31mdx_tbl.sparse[%"PRIi64"] is wrong\e[0m\n", dx_tbl->x/n);
			}
		}
	}
}

static void verify_series_product(uint64_t n, int64_t m, int64_t f[restrict static n], int64_t g[restrict static n], int64_t h[restrict static n]){
	for(uint64_t e = 0; e < n; ++e){
		int64_t term = 0;
		for(uint64_t k = 0; k <= e; ++k){
			term += g[k]*h[e - k];
			if(m){
				term %= m;
			}
		}
		if(m && term < 0){
			term += m;
		}
		if(term != (m ? nut_i64_mod(f[e], m) : f[e])){
			fprintf(stderr, "\e[1;31mf_%"PRIu64" does not match g <*> h\e[0m\n", e);
		}
	}
}

static void test_countprimes(uint64_t sieve_max){
	fprintf(stderr, "\e[1;34mTesting prime counting via powerful number trick\e[0m\n");
	// These constants were chosen as the set of primes 61 < p < 256 such that their product
	// is maximized without exceeding 2^64.  This allows us to compute any number n mod all of
	// these primes and then use the chinese remainder theorem to reconstruct n mod their product.
	// The lower bound 61 comes because the highest power below 2^64 is 2^63, so we only need 63
	// terms in the bell series, meaning that we need the numbers from 1 to 63 to be invertable,
	// so we need a prime larger than 61, the largest prime up to 63.
	// The upper bound 256 comes because we will work mod x^4 for some parts of the problem, so to avoid
	// overflow when multiplying we need x^8 < 2^64.
	// By a bit of luck, these also won't cause overflow on signed 64 bit integers.
	// This set was found using a meet in the middle approach in python:
	// - generate all primes on (61, 256)
	// - take their logarithms
	// - use a greedy recursive algorithm to find all possible subset sums <= log(2^64) of the second half of the logarithms
	// - sort the list of all possible subset sums of the second half
	// - now we can use a similar greedy recursive algorithm to build subset sums for each subset of the first half of the logarithms
	// -  and for each subset sum that doesn't go over the cap, we find the largest subset sum in the list of subset sums for the second
	//    half using binary search
	// this runs in a couple seconds in python.  Decimal must be used to get enough precision to reconstruct the max product from
	// the max sum without having to do extra sieving and searching at the end
	static const int64_t xs[] = {71, 73, 127, 137, 149, 163, 179, 211, 223};
	uint64_t ld_counts[4] = {[0]=1};
	uint64_t *sigma_0_vals [[gnu::cleanup(cleanup_free)]] = nut_sieve_sigma_0(sieve_max);
	check_alloc("divisor count sieve", sigma_0_vals);
	for(uint64_t n = 2; n <= sieve_max; ++n){
		switch(sigma_0_vals[n]){
			case 2: ld_counts[1]++; break;
			case 4: ld_counts[2]++; break;
			case 8: ld_counts[3]++;
		}
	}
	fprintf(stderr, "\e[1;34mCorrect counts are:\n#{n:d(n)=2}: %"PRIu64"\n#{n:d(n)=4}: %"PRIu64"\n#{n:d(n)=8}: %"PRIu64"\e[0m\n", ld_counts[1], ld_counts[2], ld_counts[3]);
	uint64_t bits = 64 - __builtin_clzll(sieve_max);
	uint64_t max_denom = 223 - 1;
	int64_t *h_vals [[gnu::cleanup(cleanup_free)]] = malloc((bits + 1)*sizeof(int64_t));
	int64_t *f_vals [[gnu::cleanup(cleanup_free)]] = malloc((bits + 1)*sizeof(int64_t));
	int64_t *g_vals [[gnu::cleanup(cleanup_free)]] = malloc((bits + 1)*sizeof(int64_t));
	uint64_t *factorials [[gnu::cleanup(cleanup_free)]] = malloc((bits + 223)*sizeof(uint64_t));
	uint64_t *inv_factorials [[gnu::cleanup(cleanup_free)]] = malloc((max_denom + 1)*sizeof(uint64_t));
	check_alloc("perturbative bell series", h_vals);
	check_alloc("target bell series", f_vals);
	check_alloc("dk bell series", g_vals);
	check_alloc("factorial table", factorials);
	check_alloc("inv factorial table", inv_factorials);
	h_vals[0] = 1;
	h_vals[1] = 0;
	nut_Diri dx_tbl [[gnu::cleanup(nut_Diri_destroy)]], f_tbl [[gnu::cleanup(nut_Diri_destroy)]], g_tbl [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri_init(&dx_tbl, sieve_max, 0);
	nut_Diri_init(&f_tbl, sieve_max, 0);
	nut_Diri_init(&g_tbl, sieve_max, 0);
	check_alloc("dx table", dx_tbl.buf);
	check_alloc("f table", f_tbl.buf);
	check_alloc("g table", g_tbl.buf);
	for(uint64_t i = 0; i < 9; ++i){
		int64_t x = xs[i];
		int64_t modulus = x*x*x*x;
		fprintf(stderr, "\e[1;34mx = %"PRIi64"\e[0m\n", x);
		nut_Diri_compute_dk(&dx_tbl, x, modulus, &f_tbl, &g_tbl);
		verify_dk(&dx_tbl, x, modulus);
		nut_u64_make_factorial_tbl(223, modulus, bits, 223 - 1, factorials, inv_factorials);
		memset(g_vals, 0, (bits + 1)*sizeof(int64_t));
		for(int64_t j = 0, t = 1; t <= dx_tbl.x; ++j, t <<= 1){
			g_vals[j] = factorials[j + x - 1]*inv_factorials[x - 1]%(uint64_t)modulus*inv_factorials[j]%(uint64_t)modulus;
		}
		memset(f_vals, 0, (bits + 1)*sizeof(int64_t));
		f_vals[0] = 1;
		for(uint64_t i = 1, px = x; i <= 3; ++i, px *= x){
			f_vals[(1ull << i) - 1] = px;
			nut_series_div(bits + 1, modulus, h_vals, f_vals, g_vals);
			verify_series_product(bits + 1, modulus, f_vals, g_vals, h_vals);
			nut_PfIt pf_it;
			NUT_PfIt_INIT(&pf_it, sieve_max, modulus, h_vals);
			check_alloc("powerful iterator", pf_it.entries);
			int64_t res;
			nut_Diri_sum_adjusted(&res, &dx_tbl, &pf_it);
			res = nut_i64_mod(res, px*x);
			int64_t res_xdigits[4] = {nut_i64_mod(res, x)};
			for(int64_t px = x, j = 1; px != modulus; ++j){
				px *= x;
				res_xdigits[j] = nut_i64_mod(res, px)/(px/x);
			}
			int64_t expect_xdigits[4] = {[0] = 1};
			for(uint64_t j = 1; j <= i; ++j){
				expect_xdigits[j] = nut_i64_mod(ld_counts[j], x);
			}
			bool status = !memcmp(res_xdigits, expect_xdigits, 4*sizeof(int64_t));
			fprintf(stderr, "%sGot sum of P_1(n)(x) = ", status ? "\e[1;32m" : "\e[1;31m");
			status = false;
			for(uint64_t j = 4; j > 0;){
				--j;
				if(res_xdigits[j]){
					if(j > 1){
						fprintf(stderr, "%s%"PRIi64"x^%"PRIu64, status ? " + " : "", res_xdigits[j], j);
					}else if(j == 1){
						fprintf(stderr, "%s%"PRIi64"x", status ? " + " : "", res_xdigits[j]);
					}else{
						fprintf(stderr, "%s%"PRIi64, status ? " + " : "", res_xdigits[j]);
					}
					status = true;
				}
			}
			fprintf(stderr, "%s\e[0m\n", status ? "" : "0");
			nut_PfIt_destroy(&pf_it);
		}
	}
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

