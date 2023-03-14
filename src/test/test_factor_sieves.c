#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>

#include <nut/factorization.h>
#include <nut/sieves.h>

int main(){
	fprintf(stderr, "\e[1;34mSieving Factorizations  up to 1000...\e[0m\n");
	uint64_t w = max_prime_divs(1000);
	uint64_t *largest_factors = sieve_largest_factors(1000);
	factors_t *factors = init_factors_t_w(w);
	bool passed = true;
	for(uint64_t i = 2; i <= 1000; ++i){
		fill_factors_from_largest(factors, i, largest_factors);
		if(factors_product(factors) != i){
			fprintf(stderr, "\e[1;31mFactorization for %"PRIu64" doesn't match!\e[0m\n", i);
			passed = false;
			break;
		}
	}
	free(largest_factors);
	free(factors);
	if(passed){
		fprintf(stderr, "\e[1;32mPASSED\e[0m\n");
	}else{
		fprintf(stderr, "\e[1;31mFAILED\e[0m\n");
	}
}

