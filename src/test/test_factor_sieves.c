#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

#include <nut/factorization.h>
#include <nut/sieves.h>
#include <nut/debug.h>

int main(){
	fprintf(stderr, "\e[1;34mSieving Factorizations  up to 1000...\e[0m\n");
	uint64_t w = nut_max_prime_divs(1000);
	uint64_t *largest_factors [[gnu::cleanup(cleanup_free)]] = nut_sieve_largest_factors(1000);
	uint32_t *smallest_factors [[gnu::cleanup(cleanup_free)]] = nut_sieve_smallest_factors(1000);
	uint8_t *pdiv_counts [[gnu::cleanup(cleanup_free)]] = nut_sieve_omega(1000);
	nut_Factors *factors [[gnu::cleanup(cleanup_free)]] = nut_make_Factors_w(w);
	bool passed = true;
	for(uint64_t i = 1; i <= 1000; ++i){
		nut_fill_factors_from_largest(factors, i, largest_factors);
		if(nut_Factors_prod(factors) != i){
			fprintf(stderr, "\e[1;31mFactorization from lfs for %"PRIu64" doesn't match!\e[0m\n", i);
			passed = false;
			break;
		}
		nut_fill_factors_from_smallest(factors, i, smallest_factors);
		if(nut_Factors_prod(factors) != i){
			fprintf(stderr, "\e[1;31mFactorization from sfs for %"PRIu64" doesn't match!\e[0m\n", i);
			passed = false;
			break;
		}
		if(pdiv_counts[i] != factors->num_primes){
			fprintf(stderr, "\e[1;31mDistinct prime divisor count for %"PRIu64" doesn't match!\e[0m\n", i);
		}
	}
	if(passed){
		fprintf(stderr, "\e[1;32mPASSED\e[0m\n");
	}else{
		fprintf(stderr, "\e[1;31mFAILED\e[0m\n");
	}
}

