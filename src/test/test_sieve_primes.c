#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#include <nut/sieves.h>
#include <nut/debug.h>

static const uint64_t sieve_max = 1000000;
static uint64_t w;
static void *fzn_buf;
size_t fzn_pitch;

static void test_prime_sieve(){
	uint64_t num_primes;
	bool succeeded = true;
	fprintf(stderr, "\e[1;34mVerifying sieved primes up to %"PRIu64" using dmr...\e[0m\n", sieve_max);
	uint64_t *primes = nut_sieve_primes(sieve_max, &num_primes);
	check_alloc("prime sieve", primes);
	if(!num_primes){
		fprintf(stderr, "\e[1;31mSieving failed to generate any primes!\e[0m\n");
		exit(0);
	}else if(primes[0] != 2){
		fprintf(stderr, "\e[1;31mprimes[0] == %"PRIu64"\e[0m\n", primes[0]);
		succeeded = false;
	}
	uint64_t n = 3, i = 1;
	while(1){
		while(n < primes[i]){
			if(nut_u64_is_prime_dmr(n)){
				fprintf(stderr, "\e[1;31mSieve missed prime %"PRIu64"\e[0m\n", n);
				succeeded = false;
			}
			++n;
		}
		if(!nut_u64_is_prime_dmr(n)){
			fprintf(stderr, "\e[1;31mSieve output %"PRIu64" is not prime\e[0m\n", n);
			succeeded = false;
		}
		++n;
		++i;
		if(i == num_primes){
			break;
		}else if(primes[i] <= primes[i - 1]){
			succeeded = false;
			do{
				fprintf(stderr, "\e[1;31mSieve output is not strictly increasing (i=%"PRIu64"-%"PRIu64")\e[0m\n", i - 1, i);
				++i;
			}while(i < num_primes && primes[i] <= primes[i - 1]);
			if(i == num_primes){
				break;
			}
		}
		if(primes[i] > sieve_max){
			succeeded = false;
			do{
				fprintf(stderr, "\e[1;31mSieve output %"PRIu64" is over the upper bound\e[0m\n", primes[i]);
				++i;
			}while(i < num_primes);
			break;
		}
	}
	while(n <= sieve_max){
		if(nut_u64_is_prime_dmr(n)){
			fprintf(stderr, "\e[1;31mSieve missed prime %"PRIu64"\e[0m\n", n);
			succeeded = false;
		}
		++n;
	}
	free(primes);
	if(succeeded){
		fprintf(stderr, "\e[1;32mSieve found primes up to %"PRIu64" correctly\e[0m\n", sieve_max);
	}
}

static bool check_factorization(const nut_Factors *factors, uint64_t n){
	if(nut_Factors_prod(factors) != n){
		return false;
	}
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		if(!nut_u64_is_prime_dmr(factors->factors[i].prime)){
			return false;
		}
	}
	return true;
}

static void test_factorization_sieve(){
	fprintf(stderr, "\e[1;34mVerifying factorization sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	fzn_buf = nut_sieve_factorizations(sieve_max, &w);
	check_alloc("factorization sieve", fzn_buf);
	static const nut_Factors dummy;
	fzn_pitch = offsetof(nut_Factors, factors) + w*sizeof(dummy.factors[0]);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		const nut_Factors *factors = nut_Pitcharr_get(fzn_buf, fzn_pitch, n);
		if(check_factorization(factors, n)){
			++correct;
		}
	}
	print_summary("factorizations", correct, sieve_max);
}

static void test_factor_sieve(){
	fprintf(stderr, "\e[1;34mVerifying factor sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	void *fs_buf = nut_sieve_factors(sieve_max, &w);
	check_alloc("nut_sieve_factors", fs_buf);
	size_t fs_pitch = offsetof(nut_u64_Pitcharr, elems) + w*sizeof(uint64_t);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		const nut_Factors *factors = nut_Pitcharr_get(fzn_buf, fzn_pitch, n);
		const nut_u64_Pitcharr *arr = nut_Pitcharr_get(fs_buf, fs_pitch, n);
		if(arr->len != factors->num_primes){
			continue;
		}
		int all_factors_match = 1;
		for(uint64_t i = 0; i < factors->num_primes; ++i){
			if(factors->factors[i].prime != arr->elems[i]){
				all_factors_match = 0;
				break;
			}
		}
		correct += !!all_factors_match;
	}
	print_summary("factors", correct, sieve_max);
	free(fs_buf);
}

#define TEST_FUNCTION_SIEVE(sieve, single, sieve_name, single_name, plural_name, ...) do{\
	fprintf(stderr, "\e[1;34mVerifying %s up to %"PRIu64"...\e[0m\n", sieve_name, sieve_max);\
	uint64_t *tmp_buf [[gnu::cleanup(cleanup_free)]] = sieve(sieve_max __VA_OPT__(,) __VA_ARGS__);\
	check_alloc(sieve_name, tmp_buf);\
	uint64_t correct = 0;\
	for(uint64_t n = 1; n <= sieve_max; ++n){\
		const nut_Factors *factors = nut_Pitcharr_get(fzn_buf, fzn_pitch, n);\
		uint64_t tmp = single(factors __VA_OPT__(,) __VA_ARGS__);\
		if(tmp != tmp_buf[n]){\
			fprintf(stderr, "\e[1;31m%s mismatch at %"PRIu64"\e[0m\n", single_name, n);\
			continue;\
		}\
		++correct;\
	}\
	print_summary(plural_name, correct, sieve_max);\
}while(0)

static void test_largest_factor_sieve(){
	fprintf(stderr, "\e[1;34mVerifying largest factor sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	uint64_t *largest_factors = nut_sieve_largest_factors(sieve_max);
	check_alloc("largest_factors", largest_factors);
	nut_Factors *factors = nut_make_Factors_ub(sieve_max, 25, nut_small_primes);
	check_alloc("largest_factors(init_factors_t_ub)", factors);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		nut_fill_factors_from_largest(factors, n, largest_factors);
		if(check_factorization(factors, n)){
			++correct;
		}
	}
	free(largest_factors);
	free(factors);
	print_summary("largest factors", correct, sieve_max);
}

int main(){
	test_prime_sieve();
	test_factorization_sieve();
	test_factor_sieve();
	TEST_FUNCTION_SIEVE(nut_sieve_sigma_0, nut_Factor_divcount, "nut_sieve_sigma_0", "divisor count", "divisor counts");
	TEST_FUNCTION_SIEVE(nut_sieve_sigma_1, nut_Factor_divsum, "nut_sieve_sigma_1", "divisor sum", "divisor sums");
	TEST_FUNCTION_SIEVE(nut_sieve_sigma_e, nut_Factor_divpowsum, "nut_sieve_sigma_e(2)", "divisor square sum", "divisor square sums", 2);
	TEST_FUNCTION_SIEVE(nut_sieve_dk, nut_Factor_divtupcount, "nut_sieve_dk(3)", "divisor triple counts", "divisor square sums", 3, 0);
	TEST_FUNCTION_SIEVE(nut_sieve_dk, nut_Factor_divtupcount, "nut_sieve_dk(3) mod 65521", "divisor triple counts", "divisor square sums", 3, 65521);
	TEST_FUNCTION_SIEVE(nut_sieve_phi, nut_Factor_phi, "nut_sieve_phi", "euler phi", "euler phi");
	TEST_FUNCTION_SIEVE(nut_sieve_carmichael, nut_Factor_carmichael, "nut_sieve_carmichael", "carmichael lambda", "carmichael lambda");
	free(fzn_buf);
	test_largest_factor_sieve();
}

