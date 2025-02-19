#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>

#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/sieves.h>
#include <nut/segsieves.h>

bool nut_Segsieve_init(nut_Segsieve *self, uint64_t max, uint64_t preferred_bucket_size){
	self->max = max;
	self->sqrt_max = nut_u64_nth_root(max, 2);
	self->preferred_bucket_size = preferred_bucket_size ?: self->sqrt_max;
	return (self->primes = nut_sieve_primes(self->sqrt_max, &self->num_primes));
}

void nut_Segsieve_destroy(nut_Segsieve *self){
	free(self->primes);
}

void nut_Segsieve_factorizations(const nut_Segsieve *restrict self, uint64_t a, uint64_t b, size_t pitch, void *buffer){
	for(uint64_t n = a; n < b; ++n){
		nut_Factors *fxn = nut_Pitcharr_get(buffer, pitch, n - a);
		fxn->num_primes = 1;
		fxn->factors[0].prime = n;
		fxn->factors[0].power = 1;
	}
	for(uint64_t i = 0; i < self->num_primes; ++i){
		uint64_t p = self->primes[i];
		for(uint64_t m = (a + p - 1)/p*p; m < b; m += p){
			nut_Factors *fxn = nut_Pitcharr_get(buffer, pitch, m - a);
			fxn->factors[fxn->num_primes].prime = p;
			fxn->factors[fxn->num_primes].power = 1;
			fxn->factors[0].prime /= p;
			while(fxn->factors[0].prime%p == 0){
				fxn->factors[0].prime /= p;
				fxn->factors[fxn->num_primes].power++;
			}
			fxn->num_primes++;
		}
	}
}

void *nut_Segsieve_factorizations_mkbuffer(const nut_Segsieve *self, size_t *pitch){
	if(!*pitch){
		*pitch = offsetof(nut_Factors, factors) + (nut_max_prime_divs(self->max) + 1)*sizeof(uint64_t)*2;
	}
	return malloc(*pitch*self->preferred_bucket_size);
}

