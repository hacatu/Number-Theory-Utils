#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

#include <nut/factorization.h>

static int dvisitor_coll_moments(const nut_Factors *dfactors, uint64_t d, void *data){
	uint64_t *moments_out = data;
	moments_out[0]++;
	moments_out[1] += d;
	moments_out[2] += d*d;
	return 0;
}

static void naive_coll_moments(uint64_t n, uint64_t d_max, uint64_t moments_out[static 3]){
	memset(moments_out, 0, 3*sizeof(uint64_t));
	for(uint64_t d = 1; d <= d_max; ++d){
		if(n%d == 0){
			moments_out[0]++;
			moments_out[1] += d;
			moments_out[2] += d*d;
		}
	}
}

int main(){
	uint64_t brute_moments[3] = {};
	uint64_t visitor_moments[3] = {};
	uint64_t trials = 1000, passed = 0;
	nut_Factors *factors = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	fprintf(stderr, "\e[1;34mGenerating bounded divisors for %"PRIu64" random numbers...\e[0m\n", trials);
	for(uint64_t i = 0; i < trials; ++i){
		uint64_t n = nut_u64_rand(2, 10000);
		uint64_t d_max = nut_u64_rand(1, n + 1);
		if(nut_u64_factor_heuristic(n, 25, nut_small_primes, &nut_default_factor_conf, factors) != 1){
			fprintf(stderr, "\e[1;31mFailed to factor %"PRIu64"!\e[0m\n", n);
			continue;
		}else if(nut_Factors_prod(factors) != n){
			fprintf(stderr, "\e[1;31mProduct of factorization doesn't match for %"PRIu64"\e[0m\n", n);
			continue;
		}
		bool status = true;
		for(uint64_t i = 0; i < factors->num_primes; ++i){
			if(!nut_u64_is_prime_dmr(factors->factors[i].prime)){
				status = 0;
				break;
			}
		}
		if(!status){
			fprintf(stderr, "\e[1;31mNot all factors are prime for %"PRIu64"\e[0m\n", n);
			continue;
		}
		memset(brute_moments, 0, 3*sizeof(uint64_t));
		memset(visitor_moments, 0, 3*sizeof(uint64_t));
		naive_coll_moments(n, d_max, brute_moments);
		nut_Factor_forall_divs_le_tmptmp(factors, d_max, dvisitor_coll_moments, visitor_moments);
		if(memcmp(brute_moments, visitor_moments, 3*sizeof(uint64_t))){
			fprintf(stderr, "\e[1;31mDivisor mismatch for d | %"PRIu64", d < %"PRIu64"\e[0m\n", n, d_max);
		}
		++passed;
	}
	free(factors);
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
}

