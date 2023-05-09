#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#include <nut/sieves.h>

static const uint64_t sieve_max = 1000000;
static uint64_t w;
static void *fzn_buf;
size_t fzn_pitch;

static void check_alloc(const char *what, const void *buf){
	if(!buf){
		fprintf(stderr, "\e[1;31mAllocation failed for %s!\e[0m\n", what);
		exit(0);
	}
}

static void print_summary(const char *what, uint64_t correct){
	fprintf(stderr, "%s (found %s for %"PRIu64"/%"PRIu64" numbers correctly)\e[0m\n", correct == sieve_max ? "\e[1;32mPASSED" : "\e[1;31mFAILED", what, correct, sieve_max);
}

static void test_prime_sieve(){
	uint64_t num_primes;
	bool succeeded = true;
	fprintf(stderr, "\e[1;34mVerifying sieved primes up to %"PRIu64" using dmr...\e[0m\n", sieve_max);
	uint64_t *primes = sieve_primes(sieve_max, &num_primes);
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
			if(is_prime_dmr(n)){
				fprintf(stderr, "\e[1;31mSieve missed prime %"PRIu64"\e[0m\n", n);
				succeeded = false;
			}
			++n;
		}
		if(!is_prime_dmr(n)){
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
		if(is_prime_dmr(n)){
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

static bool check_factorization(const factors_t *factors, uint64_t n){
	if(factors_product(factors) != n){
		return false;
	}
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		if(!is_prime_dmr(factors->factors[i].prime)){
			return false;
		}
	}
	return true;
}

static void test_factorization_sieve(){
	fprintf(stderr, "\e[1;34mVerifying factorization sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	fzn_buf = sieve_factorizations(sieve_max, &w);
	check_alloc("factorization sieve", fzn_buf);
	static const factors_t dummy;
	fzn_pitch = offsetof(factors_t, factors) + w*sizeof(dummy.factors[0]);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		const factors_t *factors = pitch_arr_get(fzn_buf, fzn_pitch, n);
		if(check_factorization(factors, n)){
			++correct;
		}
	}
	print_summary("factorizations", correct);
}

static void test_factor_sieve(){
	fprintf(stderr, "\e[1;34mVerifying factor sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	void *fs_buf = sieve_factors(sieve_max, &w);
	check_alloc("sieve_factors", fs_buf);
	size_t fs_pitch = offsetof(fw_u64arr_t, elems) + w*sizeof(uint64_t);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		const factors_t *factors = pitch_arr_get(fzn_buf, fzn_pitch, n);
		const fw_u64arr_t *arr = pitch_arr_get(fs_buf, fs_pitch, n);
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
	print_summary("factors", correct);
	free(fs_buf);
}

static void test_function_sieve(uint64_t *(*sieve)(uint64_t max), uint64_t (*single)(const factors_t *factors), const char *sieve_name, const char *single_name, const char *plural_name){
	fprintf(stderr, "\e[1;34mVerifying %s up to %"PRIu64"...\e[0m\n", sieve_name, sieve_max);
	uint64_t *tmp_buf = sieve(sieve_max);
	check_alloc(sieve_name, tmp_buf);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		const factors_t *factors = pitch_arr_get(fzn_buf, fzn_pitch, n);
		uint64_t tmp = single(factors);
		if(tmp != tmp_buf[n]){
			fprintf(stderr, "\e[1;31m%s mismatch at %"PRIu64"\e[0m\n", single_name, n);
			continue;
		}
		++correct;
	}
	free(tmp_buf);
	print_summary(plural_name, correct);
}

static void test_function1_sieve(uint64_t *(*sieve)(uint64_t max, uint64_t x), uint64_t (*single)(const factors_t *factors, uint64_t x), uint64_t x, const char *sieve_name, const char *single_name, const char *plural_name){
	fprintf(stderr, "\e[1;34mVerifying %s up to %"PRIu64"...\e[0m\n", sieve_name, sieve_max);
	uint64_t *tmp_buf = sieve(sieve_max, x);
	check_alloc(sieve_name, tmp_buf);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		const factors_t *factors = pitch_arr_get(fzn_buf, fzn_pitch, n);
		uint64_t tmp = single(factors, x);
		if(tmp != tmp_buf[n]){
			fprintf(stderr, "\e[1;31m%s mismatch at %"PRIu64"\e[0m\n", single_name, n);
			continue;
		}
		++correct;
	}
	free(tmp_buf);
	print_summary(plural_name, correct);
}

static void test_largest_factor_sieve(){
	fprintf(stderr, "\e[1;34mVerifying largest factor sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	uint64_t *largest_factors = sieve_largest_factors(sieve_max);
	check_alloc("largest_factors", largest_factors);
	factors_t *factors = init_factors_t_ub(sieve_max, 25, small_primes);
	check_alloc("largest_factors(init_factors_t_ub)", factors);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		fill_factors_from_largest(factors, n, largest_factors);
		if(check_factorization(factors, n)){
			++correct;
		}
	}
	free(largest_factors);
	free(factors);
	print_summary("largest factors", correct);
}

int main(){
	test_prime_sieve();
	test_factorization_sieve();
	test_factor_sieve();
	test_function_sieve(sieve_sigma_0, divisor_count, "sieve_sigma_0", "divisor count", "divisor counts");
	test_function_sieve(sieve_sigma_1, divisor_sum, "sieve_sigma_1", "divisor sum", "divisor sums");
	test_function1_sieve(sieve_sigma_e, divisor_power_sum, 2, "sieve_sigma_e(2)", "divisor square sum", "divisor square sums");
	test_function1_sieve(sieve_dk, divisor_tuple_count, 3, "sieve_dk(3)", "divisor triple counts", "divisor square sums");
	test_function_sieve(sieve_phi, euler_phi, "sieve_phi", "euler phi", "euler phi");
	test_function_sieve(sieve_carmichael, carmichael_lambda, "sieve_carmichael", "carmichael lambda", "carmichael lambda");
	free(fzn_buf);
	test_largest_factor_sieve();
}

