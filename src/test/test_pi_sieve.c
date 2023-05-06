#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#include <nut/sieves.h>

static const uint64_t sieve_max = 1000000;

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
	uint64_t correct = 0, count = 0;
	fprintf(stderr, "\e[1;34mVerifying is_composite and compute_pi_from_tables up to %"PRIu64" using dmr...\e[0m\n", sieve_max);
	uint8_t *buf = sieve_is_composite(sieve_max);
	uint64_t *pi_table = compute_pi_range(sieve_max, buf);
	check_alloc("prime sieve", buf);
	check_alloc("pi range", pi_table);
	if(!is_composite(1, buf)){
		fprintf(stderr, "\e[1;31m1 should be composite!\e[0m\n");
	}else if(compute_pi_from_tables(1, pi_table, buf)){
		fprintf(stderr, "\e[1;31mpi(1) should be 0!\e[0m\n");
	}else{
		++correct;
	}
	if(is_composite(2, buf)){
		fprintf(stderr, "\e[1;31m2 should be prime!\e[0m\n");
	}else if(compute_pi_from_tables(2, pi_table, buf) != 1){
		fprintf(stderr, "\e[1;31mpi(2) should be 1!\e[0m\n");
	}else{
		++correct;
		++count;
	}
	for(uint64_t n = 3; n <= sieve_max; ++n){
		bool is_prime = is_prime_dmr(n);
		if(is_prime){
			++count;
		}
		if(is_prime == is_composite(n, buf)){
			fprintf(stderr, "\e[1;31mis_composite(%"PRIu64") is wrong\e[0m\n", n);
		}else if(compute_pi_from_tables(n, pi_table, buf) != count){
			fprintf(stderr, "\e[1;31mpi(%"PRIu64") should be %"PRIu64"!\e[0m\n", n, count);
		}else{
			++correct;
		}
	}
	free(pi_table);
	free(buf);
	print_summary("is_composite", correct);
}

int main(){
	test_prime_sieve();
}

