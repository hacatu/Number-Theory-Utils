#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>

#include <nut/factorization.h>

int main(){
	nut_Factors *factors_n = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	nut_Factors *factors_m = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	uint64_t trials = 500, passed = 0;
	fprintf(stderr, "\e[1;34mTesting combine factors for %"PRIu64" pairs of random numbers...\e[0m\n", trials);
	for(uint64_t i = 0; i < trials; ++i){
		uint64_t n = nut_u64_rand(1ull << 29, 1ull << 30);
		uint64_t m = nut_u64_rand(1ull << 29, 1ull << 30);
		if(nut_u64_factor_heuristic(n, 25, nut_small_primes, &nut_default_factor_conf, factors_n) != 1){
			fprintf(stderr, "\e[1;31mFailed to factor %"PRIu64"!\e[0m\n", n);
			continue;
		}else if(nut_u64_factor_heuristic(m, 25, nut_small_primes, &nut_default_factor_conf, factors_m) != 1){
			fprintf(stderr, "\e[1;31mFailed to factor %"PRIu64"!\e[0m\n", m);
			continue;
		}else if(nut_Factors_prod(factors_n) != n){
			fprintf(stderr, "\e[1;31mProduct of factorization doesn't match for %"PRIu64"\e[0m\n", n);
			continue;
		}else if(nut_Factors_prod(factors_n) != n){
			fprintf(stderr, "\e[1;31mProduct of factorization doesn't match for %"PRIu64"\e[0m\n", m);
			continue;
		}
		int all_factors_prime = 1;
		for(uint64_t i = 0; i < factors_n->num_primes; ++i){
			if(!nut_u64_is_prime_dmr(factors_n->factors[i].prime)){
				all_factors_prime = 0;
				break;
			}
		}
		if(!all_factors_prime){
			fprintf(stderr, "\e[1;31mNot all factors are prime for %"PRIu64"\e[0m\n", n);
			continue;
		}
		for(uint64_t i = 0; i < factors_m->num_primes; ++i){
			if(!nut_u64_is_prime_dmr(factors_m->factors[i].prime)){
				all_factors_prime = 0;
				break;
			}
		}
		if(!all_factors_prime){
			fprintf(stderr, "\e[1;31mNot all factors are prime for %"PRIu64"\e[0m\n", m);
			continue;
		}
		nut_Factor_combine(factors_n, factors_m, 1);
		if(nut_Factors_prod(factors_n) != n*m){
			fprintf(stderr, "\e[1;31mProduct of combined factorization doesn't match for %"PRIu64"*%"PRIu64"\e[0m\n", n, m);
			continue;
		}
		nut_Factor_ipow(factors_m, 2);
		if(nut_Factors_prod(factors_m) != m*m){
			fprintf(stderr, "\e[1;31mPower of factorization doesn't match for %"PRIu64"^2\e[0m\n", m);
			continue;
		}
		++passed;
	}
	free(factors_n);
	free(factors_m);
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
}

