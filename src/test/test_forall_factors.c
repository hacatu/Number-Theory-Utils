#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

#include <nut/factorization.h>

static int dvisitor_print(const factors_t *dfactors, uint64_t d, void *data){
	factors_fprint(stdout, dfactors);
	putc('\n', stdout);
	return 0;
}

int main(){
	factors_t *factors = init_factors_t_w(MAX_PRIMES_64);
	uint64_t n = 360;
	if(factor_heuristic(n, 25, small_primes, &default_factor_conf, factors) != 1){
		fprintf(stderr, "\e[1;31mFailed to factor %"PRIu64"!\e[0m\n", n);
		exit(1);
	}else if(factors_product(factors) != n){
		fprintf(stderr, "\e[1;31mProduct of factorization doesn't match for %"PRIu64"\e[0m\n", n);
		exit(1);
	}
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		if(!is_prime_dmr(factors->factors[i].prime)){
			fprintf(stderr, "\e[1;31mNot all factors are prime for %"PRIu64"\e[0m\n", n);
			exit(1);
		}
	}
	forall_divisors_tmptmp(factors, dvisitor_print, NULL);
	free(factors);
}

