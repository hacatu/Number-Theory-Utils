#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>

#include "sieves.h"

int main(){
	uint64_t sieve_max = 1000000;
	fprintf(stderr, "\e[1;34mVerifying sieved primes up to %"PRIu64" using dmr...\e[0m\n", sieve_max);
	uint64_t num_primes;
	uint64_t *primes = sieve_primes(sieve_max, &num_primes);
	int succeeded = 1;
	if(!primes){
		fprintf(stderr, "\e[1;31mCould not sieve primes up to %"PRIu64"!\e[0m\n", sieve_max);
		return 0;
	}else if(!num_primes){
		fprintf(stderr, "\e[1;31mSieving failed to generate any primes!\e[0m\n");
		return 0;
	}else if(primes[0] != 2){
		fprintf(stderr, "\e[1;31mprimes[0] == %"PRIu64"\e[0m\n", primes[0]);
		succeeded = 0;
	}
	uint64_t n = 3, i = 1;
	while(1){
		while(n < primes[i]){
			if(is_prime_dmr(n)){
				fprintf(stderr, "\e[1;31mSieve missed prime %"PRIu64"\e[0m\n", n);
				succeeded = 0;
			}
			++n;
		}
		if(!is_prime_dmr(n)){
			fprintf(stderr, "\e[1;31mSieve output %"PRIu64" is not prime\e[0m\n", n);
			succeeded = 0;
		}
		++n;
		++i;
		if(i == num_primes){
			break;
		}else if(primes[i] <= primes[i - 1]){
			uint64_t j = i - 1;
			do{
				fprintf(stderr, "\e[1;31mSieve output is not strictly increasing (i=%"PRIu64"-%"PRIu64")\e[0m\n", j, i);
				succeeded = 0;
				++i;
			}while(i < num_primes && primes[i] <= primes[j]);
			if(i == num_primes){
				break;
			}
		}
		if(primes[i] > sieve_max){
			do{
				fprintf(stderr, "\e[1;31mSieve output %"PRIu64" is over the upper bound\e[0m\n", primes[i]);
				succeeded = 0;
				++i;
			}while(i < num_primes);
			break;
		}
	}
	while(n <= sieve_max){
		if(is_prime_dmr(n)){
			fprintf(stderr, "\e[1;31mSieve missed prime %"PRIu64"\e[0m\n", n);
			succeeded = 0;
		}
		++n;
	}
	free(primes);
	if(succeeded){
		fprintf(stderr, "\e[1;32mSieve found primes up to %"PRIu64" correctly\e[0m\n", sieve_max);
	}
	fprintf(stderr, "\e[1;34mVerifying factorization sieve up to %"PRIu64"...\e[0m\n", sieve_max);
	uint64_t w;
	void *fzn_buf = sieve_factorizations(sieve_max, &w);
	if(!fzn_buf){
		fprintf(stderr, "\e[1;31mCould not sieve factorizatons up to %"PRIu64"!\e[0m\n", sieve_max);
		return 0;
	}
	static const factors_t dummy;
	size_t fzn_pitch = offsetof(factors_t, factors) + w*sizeof(dummy.factors[0]);
	uint64_t correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		factors_t *factors = pitch_arr_get(fzn_buf, fzn_pitch, n);
		if(factors_product(factors) != n){
			continue;
		}
		int all_factors_prime = 1;
		for(uint64_t i = 0; i < factors->num_primes; ++i){
			if(!is_prime_dmr(factors->factors[i].prime)){
				all_factors_prime = 0;
				break;
			}
		}
		correct += !!all_factors_prime;
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" numbers factored correctly)\e[0m\n", correct == sieve_max ? "\e[1;32mPASSED" : "\e[1;31mFAILED", correct, sieve_max);
	void *fs_buf = sieve_factors(sieve_max, &w);
	if(!fs_buf){
		fprintf(stderr, "\e[1;31mCould not sieve factors up to %"PRIu64"!\e[0m\n", sieve_max);
		return 0;
	}
	size_t fs_pitch = offsetof(fw_u64arr_t, elems) + w*sizeof(uint64_t);
	correct = 0;
	for(uint64_t n = 1; n <= sieve_max; ++n){
		factors_t *factors = pitch_arr_get(fzn_buf, fzn_pitch, n);
		fw_u64arr_t *arr = pitch_arr_get(fs_buf, fs_pitch, n);
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
	fprintf(stderr, "%s (found factors for %"PRIu64"/%"PRIu64" numbers correctly)\e[0m\n", correct == sieve_max ? "\e[1;32mPASSED" : "\e[1;31mFAILED", correct, sieve_max);
	free(fs_buf);
	free(fzn_buf);
	//TODO: test sigma_0, sigma_1, sigma_e, phi, and carmichael
}

