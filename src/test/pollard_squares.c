#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "factorization.h"

int main(){
	uint64_t num_primes = 0;
	uint64_t *primes = malloc(1000*sizeof(uint64_t));
	if(!primes){
		fprintf(stderr, "\e[1;31mCould not allocate primes\e[0m\n");
		exit(EXIT_FAILURE);
	}
	for(uint64_t i = 2; i < 1000; ++i){
		if(is_prime_dmr(i)){
			primes[num_primes++] = i;
		}
	}
	uint64_t *tmp = realloc(primes, num_primes*sizeof(uint64_t));
	if(tmp){
		primes = tmp;
	}
	fprintf(stderr, "\e[1;34mFound %"PRIu64" primes below 1000\e[0m\n", num_primes);
	for(uint64_t i = 0; i < num_primes; ++i){
		uint64_t p = primes[i];
		int found_factor = 0;
		for(uint64_t x = 0; x < p; ++x){
			if(factor1_pollard_rho(p*p, x) == p){
				found_factor = 1;
				if(x){
					fprintf(stderr, "\e[1;32mSmallest pollard witness for %"PRIu64"**2 is %"PRIu64"\e[0m\n", p, x);
				}
				break;
			}
		}
		if(!found_factor){
			fprintf(stderr, "\e[1;33mNo pollard witness for %"PRIu64"**2\e[0m\n", p);
		}
	}
	free(primes);
}

