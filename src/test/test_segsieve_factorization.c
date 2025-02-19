#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <pthread.h>

#include <nut/sieves.h>
#include <nut/segsieves.h>
#include <nut/debug.h>

static bool check_factorization(const nut_Factors *factors, uint64_t n){
	if(nut_Factors_prod(factors) != n){
		return false;
	}
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		uint64_t p = factors->factors[i].prime;
		if(!((i == 0 && p == 1) || nut_u64_is_prime_dmr(p))){
			return false;
		}
	}
	return true;
}

#define NUM_THREADS 2

typedef struct{
	const nut_Segsieve *header;
	uint64_t tidx;
} worker_arg_s;

void *worker_fn(void *_args){
	worker_arg_s *args = _args;
	size_t pitch = 0;
	void *buffer [[gnu::cleanup(cleanup_free)]] = nut_Segsieve_factorizations_mkbuffer(args->header, &pitch);
	check_alloc("Thread work buffer", buffer);
	uint64_t bucket_size = args->header->preferred_bucket_size;
	uint64_t correct = 0;
	for(uint64_t a = args->tidx*bucket_size, b = a + bucket_size; a <= args->header->max; a += NUM_THREADS*bucket_size, b = a + bucket_size){
		if(b > args->header->max){
			b = args->header->max + 1;
		}
		if(a == 0){
			a = 2;
		}
		nut_Segsieve_factorizations(args->header, a, b, pitch, buffer);
		for(uint64_t n = a; n < b; ++n){
			nut_Factors *fxn = nut_Pitcharr_get(buffer, pitch, n - a);
			if(!check_factorization(fxn, n)){
				fprintf(stderr, "\e[1;31mERROR: factorization of %"PRIu64" is wrong!\e[0m\n", n);
				exit(1);
			}
			++correct;
		}
		if(a == 2){
			a = 0;
		}
	}
	fprintf(stderr, "\e[1;32mWorker thread checked %"PRIu64" factorizations\e[0m\n", correct);
	return NULL;
}

int main(){
	nut_Segsieve ssv [[gnu::cleanup(nut_Segsieve_destroy)]] = {};
	nut_Segsieve_init(&ssv, 1000000, 0);
	check_alloc("Segsieve", ssv.primes);
	pthread_t pids[NUM_THREADS] = {};
	worker_arg_s worker_args[NUM_THREADS] = {};
	for(uint64_t i = 0; i < NUM_THREADS; ++i){
		worker_args[i].header = &ssv;
		worker_args[i].tidx = i;
		pthread_create(pids + i, NULL, worker_fn, worker_args + i);
	}
	for(uint64_t i = 0; i < NUM_THREADS; ++i){
		pthread_join(pids[i], NULL);
	}
}

