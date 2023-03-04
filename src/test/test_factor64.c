#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>

#include <nut/factorization.h>

int main(){
	factors_t *factors = init_factors_t_w(MAX_PRIMES_64);
	uint64_t trials = 1000, passed = 0;
	fprintf(stderr, "\e[1;34mFactoring %"PRIu64" random numbers...\e[0m\n", trials);
	for(uint64_t i = 0; i < trials; ++i){
		uint64_t n = rand_u64(2, 1ull << 30);
		//factor_heuristic will check if the input and intermediate factors are prime
		//it returns the input n divided by all factors found, ie if n is factored completely it returns 1,
		//and it should not return > 1 unless n is greater than all the *_max fields in factor_conf and composite.
		//the primes and num_primes arguments here are primes_2_5 and 2 respectively.  this is required if
		//n could be a multiple of 4 or 25 and pollard_max > 3, since otherwise pollard rho will never find a factor
		//and factor_heuristic will loop indefinitely
		if(factor_heuristic(n, 25, small_primes, &default_factor_conf, factors) != 1){
			fprintf(stderr, "\e[1;31mFailed to factor %"PRIu64"!\e[0m\n", n);
			continue;
		}else if(factors_product(factors) != n){
			fprintf(stderr, "\e[1;31mProduct of factorization doesn't match for %"PRIu64"\e[0m\n", n);
			continue;
		}
		int all_factors_prime = 1;
		for(uint64_t i = 0; i < factors->num_primes; ++i){
			if(!is_prime_dmr(factors->factors[i].prime)){
				all_factors_prime = 0;
				break;
			}
		}
		if(!all_factors_prime){
			fprintf(stderr, "\e[1;31mNot all factors are prime for %"PRIu64"\e[0m\n", n);
			continue;
		}
		++passed;
	}
	free(factors);
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	//TODO: add tests for Pollard Rho, Montgomery ECF, and trial division
}

