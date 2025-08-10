#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>

#include <nut/debug.h>
#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/sieves.h>

uint64_t nut_max_prime_divs(uint64_t max){
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

uint64_t nut_max_primes_le(uint64_t max){
	return 1.25506*max/log(max);
}

void *nut_sieve_factorizations(uint64_t max, uint64_t *_w){
	uint64_t w = nut_max_prime_divs(max);
	size_t pitch = nut_get_factorizations_pitch(w);
	void *buf = calloc(pitch, max + 1);
	if(!buf){
		return NULL;
	}
	*_w = w;
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		nut_Factors *factors = nut_Pitcharr_get(buf, pitch, n);
		if(factors->num_primes){
			continue;
		}
		for(uint64_t p = n, e = 1; p <= max; ++e){
			for(uint64_t m = p; m <= max;){
				nut_Factors *factors = nut_Pitcharr_get(buf, pitch, m);
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

uint64_t nut_get_factorizations_pitch(uint64_t w){
	static const nut_Factors dummy;
	return offsetof(nut_Factors, factors) + w*sizeof(dummy.factors[0]);
}

void *nut_sieve_factors(uint64_t max, uint64_t *_w){
	uint64_t w = nut_max_prime_divs(max);
	size_t pitch = nut_get_factors_pitch(w);
	void *buf = calloc(max + 1, pitch);
	if(!buf){
		return NULL;
	}
	*_w = w;
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		nut_u64_Pitcharr *factors = nut_Pitcharr_get(buf, pitch, n);
		if(factors->len){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			nut_u64_Pitcharr *factors = nut_Pitcharr_get(buf, pitch, m);
			factors->elems[factors->len++] = n;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint8_t *nut_sieve_omega(uint64_t max){
	uint8_t *buf = calloc(max + 1, sizeof(uint8_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t n = 2; n <= max && n; ++n){
		if(buf[n]){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			buf[m]++;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint64_t *nut_sieve_largest_factors(uint64_t max){
	uint64_t *buf = calloc(max + 1, sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t n = 2; n <= max && n; ++n){
		if(buf[n]){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			buf[m] = n;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint32_t *nut_sieve_smallest_factors(uint64_t max){
	uint32_t *buf = calloc(max + 1, sizeof(uint32_t));
	if(!buf){
		return NULL;
	}
	uint64_t rmax = nut_u64_nth_root(max, 2);
	for(uint64_t n = 2; n <= rmax; ++n){
		if(buf[n]){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			if(buf[m] == 0 || buf[m] > n){
				buf[m] = n;
			}
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
		buf[n] = 1; // this means the first store in the above loop is dead
	}
	for(uint64_t n = rmax + 1; n <= max; ++n){
		if(!buf[n]){
			buf[n] = 1;
		}
	}
	return buf;
}

void nut_fill_factors_from_largest(nut_Factors *restrict out, uint64_t n, const uint64_t largest_factors[restrict static n + 1]){
	out->num_primes = 0;
	for(uint64_t p = largest_factors[n], k = 1; p;){
		n /= p;
		uint64_t q = largest_factors[n];
		if(p == q){
			++k;
		}else{
			nut_Factor_append(out, p, k);
			p = q;
			k = 1;
		}
	}
}

void nut_fill_factors_from_smallest(nut_Factors *restrict out, uint64_t n, const uint32_t smallest_factors[restrict static n + 1]){
	out->num_primes = 0;
	for(uint64_t p = smallest_factors[n], k = 1; p;){
		if(p == 1){
			p = n;
		}
		n /= p;
		uint64_t q = smallest_factors[n];
		if(q == 1){
			q = n;
		}
		if(p == q){
			++k;
		}else{
			nut_Factor_append(out, p, k);
			p = q;
			k = 1;
		}
	}
}

uint64_t nut_get_factors_pitch(uint64_t w){
	return offsetof(nut_u64_Pitcharr, elems) + w*sizeof(uint64_t);
}

// Works by multiplying power + 1 for all prime factors
uint64_t *nut_sieve_sigma_0(uint64_t max){
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
uint64_t *nut_sieve_sigma_1(uint64_t max){
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
uint64_t *nut_sieve_sigma_e(uint64_t max, uint64_t e){
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
			buf[m] *= (nut_u64_pow(a, e) - 1)/(nut_u64_pow(n, e) - 1);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

bool nut_u64_make_factorial_tbl(uint64_t k, uint64_t modulus, uint64_t bits, uint64_t max_denom, uint64_t factorials[static bits + k], uint64_t inv_factorials[static max_denom + 1]){
	factorials[0] = 1;
	for(uint64_t i = 1; i < bits + k; ++i){
		factorials[i] = factorials[i - 1]*i%modulus;
	}
	memcpy(inv_factorials, factorials, 2*sizeof(uint64_t));
	for(uint64_t i = 2; i <= max_denom; ++i){
		int64_t inv;
		if(nut_i64_egcd(factorials[i], modulus, &inv, NULL) != 1){
			return false;
		}else if(inv < 0){
			inv_factorials[i] = inv + (int64_t)modulus;
		}else{
			inv_factorials[i] = inv;
		}
	}
	return true;
}

static uint64_t *sieve_dk_mod(uint64_t max, uint64_t k, uint64_t modulus){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	uint8_t *is_c_buf [[gnu::cleanup(cleanup_free)]] = calloc(max/8 + 1, sizeof(uint8_t));
	uint64_t bits = 64 - __builtin_clzll(max);
	uint64_t max_denom = bits > k - 1 ? bits : k - 1;
	uint64_t *factorials [[gnu::cleanup(cleanup_free)]] = malloc((bits + k)*sizeof(uint64_t));
	uint64_t *inv_factorials [[gnu::cleanup(cleanup_free)]] = malloc((max_denom + 1)*sizeof(uint64_t));
	if(!buf || !is_c_buf || !factorials || !inv_factorials || modulus > INT64_MAX){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	if(!nut_u64_make_factorial_tbl(k, modulus, bits, max_denom, factorials, inv_factorials)){
		return NULL;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(nut_Bitarray_get(is_c_buf, n)){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			nut_Bitarray_set(is_c_buf, m, 1);
			uint64_t m1 = m, a = 0;
			while(m1%n == 0){
				m1 /= n;
				++a;
			}
			uint64_t term = buf[m]*factorials[a + k - 1]%modulus;
			term = term*inv_factorials[k - 1]%modulus;
			buf[m] = term*inv_factorials[a]%modulus;
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

static uint64_t *sieve_dk_u64(uint64_t max, uint64_t k){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	uint8_t *is_c_buf [[gnu::cleanup(cleanup_free)]] = calloc(max/8 + 1, sizeof(uint8_t));
	if(!buf || !is_c_buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	for(uint64_t n = 2; n <= max && n; ++n){//TODO: optimize loop
		if(nut_Bitarray_get(is_c_buf, n)){
			continue;
		}
		for(uint64_t m = n; m <= max;){
			nut_Bitarray_set(is_c_buf, m, 1);
			uint64_t m1 = m, a = 0;
			while(m1%n == 0){
				m1 /= n;
				++a;
			}
			buf[m] *= nut_u64_binom(a + k - 1, k - 1);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return buf;
}

uint64_t *nut_sieve_dk(uint64_t max, uint64_t k, uint64_t modulus){
	return modulus ? sieve_dk_mod(max, k, modulus) : sieve_dk_u64(max, k);
}

uint64_t *nut_sieve_phi(uint64_t max){
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

uint64_t *nut_sieve_carmichael(uint64_t max){
	uint64_t *buf = malloc((max + 1)*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	for(uint64_t i = 1; i <= max && i; ++i){
		buf[i] = 1;
	}
	uint64_t t = 4, b = 1;
	while((t >> 1) <= max){
		for(uint64_t m = t >> 1; m <= max; m += t){
			buf[m] = b;
		}
		t <<= 1;
		if(t != 16){
			b <<= 1;
		}
	}
	for(uint64_t n = 1; !__builtin_add_overflow(n, 2, &n) && n <= max;){//TODO: optimize loop
		if(buf[n] != 1){
			continue;
		}
		for(uint64_t k = 1, m = n; k <= max/n; ++k, m += n){
			uint64_t a = n - 1;
			for(uint64_t j = k; j%n == 0; j /= n){
				a *= n;
			}
			buf[m] = nut_i64_lcm(buf[m], a);
		}
	}
	return buf;
}

uint8_t *nut_sieve_mobius(uint64_t max){
	uint64_t buf_len = max/4 + 1;
	uint8_t *buf = malloc(buf_len*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	// Initialize buf to be 0 for 0, 1 for 1, and 2 for every other n <= max.
	// 2 is not a valid output value, so we use it here to mark a number as prime
	// until they are marked off by being the factor of some other number.
	buf[0] = 0xA4u;
	memset(buf + 1, 0xAA, (buf_len - 1)*sizeof(uint8_t));
	for(uint64_t n = 2; n <= max && n; ++n){
		if(((buf[n/4] >> (n%4*2)) & 3) != 2){
			continue;
		}
		buf[n/4] ^= 1ull << (n%4*2);
		uint64_t m;
		if(__builtin_mul_overflow(n, 2, &m)){
			continue;
		}
		// First we visit all multiples of n and flip their mobius value
		// We need to send 2 to 3, 1 to 3, 3 to 1, and 0 to 0.
		// We can do this by using an xor mask so we only have to
		// get x from the array and then xor the array with this mask.
		while(m <= max){
			uint8_t x = (buf[m/4] >> (m%4*2)) & 3;
			x = (0x98u >> (x*2)) & 3;// this constant acts as a packed lookup table for the xor mask for x
			buf[m/4] ^= x << (m%4*2);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
		uint64_t n2;
		// Here we visit all multiples of n**2 and set their mobius value to 0,
		// since they are not square free.
		if(!__builtin_mul_overflow(n, n, &n2)){
			for(uint64_t m = n2; m <= max;){
				buf[m/4] &= ~(3u << (m%4*2));
				if(__builtin_add_overflow(m, n2, &m)){
					break;
				}
			}
		}
	}
	return buf;
}

int64_t *nut_compute_mertens_range(uint64_t max, const uint8_t mobius[static max/4 + 1]){
	int64_t *buf = malloc((max + 1)*sizeof(int64_t));
	if(!buf){
		return NULL;
	}
	buf[0] = 0;
	for(uint64_t n = 1; n <= max; ++n){
		int64_t v = nut_Bitfield2_arr_get(mobius, n);
		if(v == 3){
			v = -1;
		}
		buf[n] = buf[n - 1] + v;
	}
	return buf;
}

static void mark_is_composite(uint8_t *is_composite, uint64_t q, uint64_t r, uint64_t q_ub){
	/* We need to find the first multiple of the prime 30*q + r which could possibly be prime and is at least the prime squared.
	 * Since the prime squared is still coprime to the wheel size, we start at the prime squared.
	 * (30*q + r)**2 = 30*30*q**2 + 30*2*q*r + r**2 = 30*(30*q**2 + 2*q*r) + r**2
	 */
	uint64_t mq;
	switch(r){
		case 1:
			mq = q*(30*q + 2) + 0;// mr = 1
			for(; mq + 28*q + 0 < q_ub; mq += 30*q + 1){
				is_composite[mq] |= 0x1;// mr = 1
				is_composite[mq + 6*q + 0] |= 0x2;// mr = 7
				is_composite[mq + 10*q + 0] |= 0x4;// mr = 11
				is_composite[mq + 12*q + 0] |= 0x8;// mr = 13
				is_composite[mq + 16*q + 0] |= 0x10;// mr = 17
				is_composite[mq + 18*q + 0] |= 0x20;// mr = 19
				is_composite[mq + 22*q + 0] |= 0x40;// mr = 23
				is_composite[mq + 28*q + 0] |= 0x80;// mr = 29
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			if((mq += 6*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x2;// mr = 7
			if((mq += 4*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			if((mq += 2*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 4*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x10;// mr = 17
			if((mq += 2*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 4*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x40;// mr = 23
			mq += 6*q + 0; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 7:
			mq = q*(30*q + 14) + 1;// mr = 19
			for(; mq + 24*q + 6 < q_ub; mq += 30*q + 7){
				is_composite[mq] |= 0x20;// mr = 19
				is_composite[mq + 4*q + 1] |= 0x10;// mr = 17
				is_composite[mq + 6*q + 2] |= 0x1;// mr = 1
				is_composite[mq + 10*q + 2] |= 0x80;// mr = 29
				is_composite[mq + 12*q + 3] |= 0x8;// mr = 13
				is_composite[mq + 16*q + 4] |= 0x4;// mr = 11
				is_composite[mq + 22*q + 5] |= 0x40;// mr = 23
				is_composite[mq + 24*q + 6] |= 0x2;// mr = 7
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 4*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x10;// mr = 17
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			if((mq += 4*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 4*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			if((mq += 6*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x40;// mr = 23
			mq += 2*q + 1; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 11:
			mq = q*(30*q + 22) + 4;// mr = 1
			for(; mq + 26*q + 9 < q_ub; mq += 30*q + 11){
				is_composite[mq] |= 0x1;// mr = 1
				is_composite[mq + 2*q + 0] |= 0x40;// mr = 23
				is_composite[mq + 6*q + 2] |= 0x2;// mr = 7
				is_composite[mq + 8*q + 2] |= 0x80;// mr = 29
				is_composite[mq + 12*q + 4] |= 0x8;// mr = 13
				is_composite[mq + 18*q + 6] |= 0x20;// mr = 19
				is_composite[mq + 20*q + 7] |= 0x4;// mr = 11
				is_composite[mq + 26*q + 9] |= 0x10;// mr = 17
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			if((mq += 2*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x40;// mr = 23
			if((mq += 4*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x2;// mr = 7
			if((mq += 2*q + 0) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 4*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 6*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			mq += 6*q + 2; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 13:
			mq = q*(30*q + 26) + 5;// mr = 19
			for(; mq + 28*q + 12 < q_ub; mq += 30*q + 13){
				is_composite[mq] |= 0x20;// mr = 19
				is_composite[mq + 4*q + 2] |= 0x4;// mr = 11
				is_composite[mq + 6*q + 3] |= 0x2;// mr = 7
				is_composite[mq + 10*q + 4] |= 0x80;// mr = 29
				is_composite[mq + 16*q + 7] |= 0x10;// mr = 17
				is_composite[mq + 18*q + 8] |= 0x8;// mr = 13
				is_composite[mq + 24*q + 11] |= 0x1;// mr = 1
				is_composite[mq + 28*q + 12] |= 0x40;// mr = 23
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 4*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x2;// mr = 7
			if((mq += 4*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 6*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x10;// mr = 17
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 6*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			mq += 4*q + 1; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 17:
			mq = q*(30*q + 34) + 9;// mr = 19
			for(; mq + 26*q + 15 < q_ub; mq += 30*q + 17){
				is_composite[mq] |= 0x20;// mr = 19
				is_composite[mq + 2*q + 1] |= 0x40;// mr = 23
				is_composite[mq + 6*q + 4] |= 0x1;// mr = 1
				is_composite[mq + 12*q + 7] |= 0x8;// mr = 13
				is_composite[mq + 14*q + 8] |= 0x10;// mr = 17
				is_composite[mq + 20*q + 11] |= 0x80;// mr = 29
				is_composite[mq + 24*q + 14] |= 0x2;// mr = 7
				is_composite[mq + 26*q + 15] |= 0x4;// mr = 11
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x40;// mr = 23
			if((mq += 4*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			if((mq += 6*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x10;// mr = 17
			if((mq += 6*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 4*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x2;// mr = 7
			mq += 2*q + 1; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 19:
			mq = q*(30*q + 38) + 12;// mr = 1
			for(; mq + 28*q + 17 < q_ub; mq += 30*q + 19){
				is_composite[mq] |= 0x1;// mr = 1
				is_composite[mq + 4*q + 2] |= 0x10;// mr = 17
				is_composite[mq + 10*q + 6] |= 0x4;// mr = 11
				is_composite[mq + 12*q + 7] |= 0x20;// mr = 19
				is_composite[mq + 18*q + 11] |= 0x8;// mr = 13
				is_composite[mq + 22*q + 13] |= 0x80;// mr = 29
				is_composite[mq + 24*q + 15] |= 0x2;// mr = 7
				is_composite[mq + 28*q + 17] |= 0x40;// mr = 23
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			if((mq += 4*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x10;// mr = 17
			if((mq += 6*q + 4) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 6*q + 4) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 4*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 2*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x2;// mr = 7
			mq += 4*q + 2; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 23:
			mq = q*(30*q + 46) + 17;// mr = 19
			for(; mq + 26*q + 20 < q_ub; mq += 30*q + 23){
				is_composite[mq] |= 0x20;// mr = 19
				is_composite[mq + 6*q + 5] |= 0x2;// mr = 7
				is_composite[mq + 8*q + 6] |= 0x40;// mr = 23
				is_composite[mq + 14*q + 11] |= 0x4;// mr = 11
				is_composite[mq + 18*q + 14] |= 0x8;// mr = 13
				is_composite[mq + 20*q + 15] |= 0x80;// mr = 29
				is_composite[mq + 24*q + 19] |= 0x1;// mr = 1
				is_composite[mq + 26*q + 20] |= 0x10;// mr = 17
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 6*q + 5) >= q_ub){break;}
			is_composite[mq] |= 0x2;// mr = 7
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x40;// mr = 23
			if((mq += 6*q + 5) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			if((mq += 4*q + 3) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 4*q + 4) >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			mq += 2*q + 1; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		case 29:
			mq = q*(30*q + 58) + 28;// mr = 1
			for(; mq + 24*q + 23 < q_ub; mq += 30*q + 29){
				is_composite[mq] |= 0x1;// mr = 1
				is_composite[mq + 2*q + 1] |= 0x80;// mr = 29
				is_composite[mq + 8*q + 7] |= 0x40;// mr = 23
				is_composite[mq + 12*q + 11] |= 0x20;// mr = 19
				is_composite[mq + 14*q + 13] |= 0x10;// mr = 17
				is_composite[mq + 18*q + 17] |= 0x8;// mr = 13
				is_composite[mq + 20*q + 19] |= 0x4;// mr = 11
				is_composite[mq + 24*q + 23] |= 0x2;// mr = 7
			}
			if(mq >= q_ub){break;}
			is_composite[mq] |= 0x1;// mr = 1
			if((mq += 2*q + 1) >= q_ub){break;}
			is_composite[mq] |= 0x80;// mr = 29
			if((mq += 6*q + 6) >= q_ub){break;}
			is_composite[mq] |= 0x40;// mr = 23
			if((mq += 4*q + 4) >= q_ub){break;}
			is_composite[mq] |= 0x20;// mr = 19
			if((mq += 2*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x10;// mr = 17
			if((mq += 4*q + 4) >= q_ub){break;}
			is_composite[mq] |= 0x8;// mr = 13
			if((mq += 2*q + 2) >= q_ub){break;}
			is_composite[mq] |= 0x4;// mr = 11
			mq += 4*q + 4; break; // If we get here, we have mq >= q_ub by the fact that we made it out of the for loop
		default:
			__builtin_unreachable();
	}
}

uint8_t *nut_sieve_is_composite(uint64_t max){
	uint64_t is_composite_len = max/30 + 1;
	uint8_t *is_composite = calloc(is_composite_len, sizeof(uint64_t));
	if(!is_composite){
		return NULL;
	}
	uint64_t q_ub = is_composite_len;
	uint64_t p_max = nut_u64_nth_root(max, 2);
	uint64_t q_max = (p_max - 1)/30;
	uint64_t q = 0, r = 7;
	while(q <= q_max){
		uint8_t flags = is_composite[q];
		switch(r){
			case 1:
				if(!(flags & 0x1)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 7; [[fallthrough]];
			case 7:
				if(!(flags & 0x2)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 11; [[fallthrough]];
			case 11:
				if(!(flags & 0x4)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 13; [[fallthrough]];
			case 13:
				if(!(flags & 0x8)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 17; [[fallthrough]];
			case 17:
				if(!(flags & 0x10)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 19; [[fallthrough]];
			case 19:
				if(!(flags & 0x20)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 23; [[fallthrough]];
			case 23:
				if(!(flags & 0x40)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 29; [[fallthrough]];
			case 29:
				if(!(flags & 0x80)){mark_is_composite(is_composite, q, r, q_ub);}
				r = 1;
				++q;
				break;
			default:
				__builtin_unreachable();
		}
	}
	return is_composite;
}

bool nut_is_composite(uint64_t n, const uint8_t buf[static n/30 + 1]){
	if(n < 6){
		return n != 2 && n != 3 && n != 5;
	}
	uint64_t r = n%30;
	uint64_t q = n/30;
	switch(r){
		case 1: return buf[q] & 0x1;
		case 7: return buf[q] & 0x2;
		case 11: return buf[q] & 0x4;
		case 13: return buf[q] & 0x8;
		case 17: return buf[q] & 0x10;
		case 19: return buf[q] & 0x20;
		case 23: return buf[q] & 0x40;
		case 29: return buf[q] & 0x80;
		default: return true;
	}
}

uint64_t *nut_compute_pi_range(uint64_t max, const uint8_t buf[static max/30 + 1]){
	/* buf[i] tells us the primality of 30*i + 1, 30*i + 7, 30*i + 11, ..., 30*i + 29.
	 * Thus, we want res[i] to tell us the number of prime numbers from 1 up to 30*(i + 1)
	 * and so in general res[i] = popcount(buf[i]) + res[i - 1].
	 * However, res[0], the number of primes up to 30, is instead 10.
	 * Finally, since max is inclusive, we must have res[max/30 - 1].
	 * The -1 is because since res[0] is the number of primes up to 30, for arguments i < 30 we have to
	 * special case their lookup, so the domain essentially starts at 30.
	 */
	uint64_t res_len = max/30;
	uint64_t *res = calloc(res_len, sizeof(uint64_t));
	if(!res){
		return NULL;
	}
	res[0] = 10;
	for(uint64_t i = 1, acc = 10; i < res_len; ++i){
		res[i] = acc += __builtin_popcount(0xFF&~buf[i]);
	}
	return res;
}

uint64_t nut_compute_pi_from_tables(uint64_t n, const uint64_t pi_table[restrict static n/30], const uint8_t buf[restrict static n/30 + 1]){
	if(n < 30){
		uint64_t res;
		for(res = 0; nut_small_primes[res] <= n; ++res);
		return res;
	}
	uint64_t q = n/30;
	uint64_t r = n%30;
	uint8_t mask;
	if(!r){mask = 0x00;}
	else if(r < 7){mask = 0x01;}
	else if(r < 11){mask = 0x03;}
	else if(r < 13){mask = 0x07;}
	else if(r < 17){mask = 0x0f;}
	else if(r < 19){mask = 0x1f;}
	else if(r < 23){mask = 0x3f;}
	else if(r < 29){mask = 0x7f;}
	else{mask = 0xff;}
	return pi_table[q-1] + __builtin_popcount(mask&~buf[q]);
}

static uint64_t *copy_small_primes(uint64_t max, uint64_t *_num_primes){
	uint64_t n = 0;
	while(n < 25 && nut_small_primes[n] <= max){
		++n;
	}
	uint64_t *primes = malloc(n*sizeof(uint64_t));
	if(!primes){
		return NULL;
	}
	memcpy(primes, nut_small_primes, n*sizeof(uint64_t));
	*_num_primes = n;
	return primes;
}

uint64_t *nut_sieve_primes(uint64_t max, uint64_t *_num_primes){
	if(max <= 100){
		return copy_small_primes(max, _num_primes);
	}
	uint64_t *primes = malloc((size_t)nut_max_primes_le(max)*sizeof(uint64_t));
	if(!primes){
		return NULL;
	}
	uint8_t *is_composite = nut_sieve_is_composite(max);
	if(!is_composite){
		free(primes);
		return NULL;
	}
	memcpy(primes, nut_small_primes, 3*sizeof(uint64_t));
	uint64_t num_primes = 3;
	uint64_t q = 0;
	uint64_t r = 7;
	while(30*q + 29 <= max){
		uint8_t flags = is_composite[q];
		switch(r){
			case 1:
				if(!(flags & 0x1)){primes[num_primes++] = 30*q + r;}
				r = 7; [[fallthrough]];
			case 7:
				if(!(flags & 0x2)){primes[num_primes++] = 30*q + r;}
				r = 11; [[fallthrough]];
			case 11:
				if(!(flags & 0x4)){primes[num_primes++] = 30*q + r;}
				r = 13; [[fallthrough]];
			case 13:
				if(!(flags & 0x8)){primes[num_primes++] = 30*q + r;}
				r = 17; [[fallthrough]];
			case 17:
				if(!(flags & 0x10)){primes[num_primes++] = 30*q + r;}
				r = 19; [[fallthrough]];
			case 19:
				if(!(flags & 0x20)){primes[num_primes++] = 30*q + r;}
				r = 23; [[fallthrough]];
			case 23:
				if(!(flags & 0x40)){primes[num_primes++] = 30*q + r;}
				r = 29; [[fallthrough]];
			case 29:
				if(!(flags & 0x80)){primes[num_primes++] = 30*q + r;}
				r = 1;
				++q;
				break;
			default:
				__builtin_unreachable();
		}
	}
	do{
		if(30*q + 1 > max){break;}
		uint8_t flags = is_composite[q];
		if(!(flags & 0x1)){primes[num_primes++] = 30*q + 1;}
		if(30*q + 7 > max){break;}
		if(!(flags & 0x2)){primes[num_primes++] = 30*q + 7;}
		if(30*q + 11 > max){break;}
		if(!(flags & 0x4)){primes[num_primes++] = 30*q + 11;}
		if(30*q + 13 > max){break;}
		if(!(flags & 0x8)){primes[num_primes++] = 30*q + 13;}
		if(30*q + 17 > max){break;}
		if(!(flags & 0x10)){primes[num_primes++] = 30*q + 17;}
		if(30*q + 19 > max){break;}
		if(!(flags & 0x20)){primes[num_primes++] = 30*q + 19;}
		if(30*q + 23 > max){break;}
		if(!(flags & 0x40)){primes[num_primes++] = 30*q + 23;}
		if(30*q + 29 > max){break;}
		if(!(flags & 0x80)){primes[num_primes++] = 30*q + 29;}
	}while(0);
	*_num_primes = num_primes;
	free(is_composite);
	return realloc(primes, num_primes*sizeof(uint64_t)) ?: primes;
}

