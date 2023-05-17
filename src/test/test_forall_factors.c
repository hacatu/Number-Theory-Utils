#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

#include <nut/factorization.h>

static int dvisitor_print(const nut_Factors *dfactors, uint64_t d, void *data){
	nut_Factor_fprint(stdout, dfactors);
	putc('\n', stdout);
	return 0;
}

int main(){
	nut_Factors *factors = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	uint64_t n = 360;
	if(nut_u64_factor_heuristic(n, 25, nut_small_primes, &nut_default_factor_conf, factors) != 1){
		fprintf(stderr, "\e[1;31mFailed to factor %"PRIu64"!\e[0m\n", n);
		exit(1);
	}else if(nut_Factors_prod(factors) != n){
		fprintf(stderr, "\e[1;31mProduct of factorization doesn't match for %"PRIu64"\e[0m\n", n);
		exit(1);
	}
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		if(!nut_u64_is_prime_dmr(factors->factors[i].prime)){
			fprintf(stderr, "\e[1;31mNot all factors are prime for %"PRIu64"\e[0m\n", n);
			exit(1);
		}
	}
	nut_Factor_forall_divs_tmptmp(factors, dvisitor_print, NULL);
	free(factors);
}

