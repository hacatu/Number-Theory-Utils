#include <stddef.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "factorization.h"
#include "sieves.h"

uint64_t max_prime_divs(uint64_t max){
	//2*3*5*7*11*13*17*19*23*29*31*37*41*43*47
	if(max < 2ull*3){//the compiler can optimize this to a binary search if it is faster
		return 1;
	}else if(max < 2ull*3*5){
		return 2;
	}else if(max < 2ull*3*5*7){
		return 3;
	}else if(max < 2ull*3*5*7*11){
		return 4;
	}else if(max < 2ull*3*5*7*11*13){
		return 5;
	}else if(max < 2ull*3*5*7*11*13*17){
		return 6;
	}else if(max < 2ull*3*5*7*11*13*17*19){
		return 7;
	}else if(max < 2ull*3*5*7*11*13*17*19*23){
		return 8;
	}else if(max < 2ull*3*5*7*11*13*17*19*23*29){
		return 9;
	}else if(max < 2ull*3*5*7*11*13*17*19*23*29*31){
		return 10;
	}else if(max < 2ull*3*5*7*11*13*17*19*23*29*31*37){
		return 11;
	}else if(max < 2ull*3*5*7*11*13*17*19*23*29*31*37*41){
		return 12;
	}else if(max < 2ull*3*5*7*11*13*17*19*23*29*31*37*41*43){
		return 13;
	}else if(max < 2ull*3*5*7*11*13*17*19*23*29*31*37*41*43*47){
		return 14;
	}
	return 15;
}

void *sieve_factorizations(uint64_t max, uint64_t *_w){
	static const factors_t dummy;
	uint64_t w = max_prime_divs(max);
	size_t pitch = offsetof(factors_t, factors) + w*sizeof(dummy.factors[0]);
	void *buf = calloc(pitch, max + 1);
	if(!buf){
		return NULL;
	}
	*_w = w;
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		factors_t *factors = pitch_arr_get(buf, pitch, n);
		if(factors->num_primes){
			continue;
		}
		for(uint64_t p = n, e = 1; p <= max; ++e){
			for(uint64_t m = p; m <= max;){
				factors_t *factors = pitch_arr_get(buf, pitch, m);
				uint64_t i = factors->num_primes;
				if(e == 1){
					++factors->num_primes;
				}else{
					--i;
				}
				factors->factors[i].prime = n;
				factors->factors[i].power = e;
				if(__builtin_add_overflow(m, p, &m)){
					break;
				}
			}
			if(__builtin_mul_overflow(p, n, &p)){
				break;
			}
		}
	}
	return buf;
}

void *sieve_factors(uint64_t max, uint64_t *_w){
	uint64_t w = max_prime_divs(max);
	size_t pitch = offsetof(fw_u64arr_t, elems) + w*sizeof(uint64_t);
	void *buf = calloc(pitch, max + 1);
	if(!buf){
		return NULL;
	}
	*_w = w;
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		fw_u64arr_t *factors = pitch_arr_get(buf, pitch, n);
		if(factors->len){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			fw_u64arr_t *factors = pitch_arr_get(buf, pitch, m);
			factors->elems[factors->len++] = n;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

// Works by multiplying power + 1 for all prime factors
uint64_t *sieve_sigma_0(uint64_t max){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(buf[n] != 1){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			uint64_t a = m/n, e = 1;
			while(a%n == 0){
				a /= n;
				++e;
			}
			buf[m] *= e + 1;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

// Works by multiplying (prime**(power+1)-1)/(prime-1) for all prime factors
uint64_t *sieve_sigma_1(uint64_t max){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(buf[n] != 1){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			uint64_t a = n*n;
			while(m%a == 0){
				a *= n;
			}
			buf[m] *= (a - 1)/(n - 1);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

// Works by multiplying (prime**((power+1)*e)-1)/(prime**e-1) for all prime factors, where power is the
// power of each prime and e is the power of divisors to sum
uint64_t *sieve_sigma_e(uint64_t max, uint64_t e){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(buf[n] != 1){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			uint64_t a = n*n;
			while(m%a == 0){
				a *= n;
			}
			buf[m] *= (pow_u64(a, e) - 1)/(pow_u64(n, e) - 1);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint64_t *sieve_phi(uint64_t max){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(buf[n] != 1){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			uint64_t a = n*n;
			while(m%a == 0){
				a *= n;
			}
			buf[m] *= a/n - a/n/n;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint64_t *sieve_carmichael(uint64_t max){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(buf[n] != 1){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			uint64_t a = n*n;
			while(m%a == 0){
				a *= n;
			}
			if(n == 2){
				if(a > 8){
					buf[m] = a >> 3;
				}else{
					buf[m] = a >> 2;
				}
			}else{
				buf[m] = lcm(buf[m], a/n - a/n/n);
			}
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint64_t *sieve_primes(uint64_t max, uint64_t *_num_primes){
	uint64_t is_composite_len = max/64 + 1;
	uint64_t *is_composite = calloc(is_composite_len, sizeof(uint64_t));
	if(!is_composite){
		return NULL;
	}
	uint64_t *primes = malloc((size_t)(1.25506*max/log(max))*sizeof(uint64_t));
	if(!primes){
		free(is_composite);
		return NULL;
	}
	uint64_t num_primes = 0;
	uint64_t n = 2;
	for(uint64_t n2; !__builtin_mul_overflow(n, n, &n2) && n2 <= max; ++n){
		if(is_composite[n/64] & (1ull << (n%64))){
			continue;
		}
		primes[num_primes++] = n;
		for(uint64_t m = n2; m <= max;){
			is_composite[m/64] |= 1ull << (m%64);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	for(; n <= max; ++n){
		if(!(is_composite[n/64] & (1ull << (n%64)))){
			primes[num_primes++] = n;
		}
	}
	free(is_composite);
	*_num_primes = num_primes;
	return realloc(primes, num_primes*sizeof(uint64_t)) ?: primes;
}

