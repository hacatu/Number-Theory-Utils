#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>

#include <nut/factorization.h>
#include <nut/sieves.h>

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

uint64_t max_primes_le(uint64_t max){
	return 1.25506*max/log(max);
}

void *sieve_factorizations(uint64_t max, uint64_t *_w){
	uint64_t w = max_prime_divs(max);
	size_t pitch = get_factorizations_pitch(w);
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

uint64_t get_factorizations_pitch(uint64_t w){
	static const factors_t dummy;
	return offsetof(factors_t, factors) + w*sizeof(dummy.factors[0]);
}

void *sieve_factors(uint64_t max, uint64_t *_w){
	uint64_t w = max_prime_divs(max);
	size_t pitch = get_factors_pitch(w);
	void *buf = calloc(max + 1, pitch);
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

uint64_t *sieve_largest_factors(uint64_t max){
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

void fill_factors_from_largest(factors_t *out, uint64_t n, const uint64_t largest_factors[static n + 1]){
	out->num_primes = 0;
	for(uint64_t p = largest_factors[n], k = 1; p;){
		n /= p;
		uint64_t q = largest_factors[n];
		if(p == q){
			++k;
		}else{
			factors_append(out, p, k);
			p = q;
			k = 1;
		}
	}
}

uint64_t get_factors_pitch(uint64_t w){
	return offsetof(fw_u64arr_t, elems) + w*sizeof(uint64_t);
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

uint64_t *sieve_mobius(uint64_t max){
	uint64_t buf_len = max/32 + 1;
	uint64_t *buf = malloc(buf_len*sizeof(uint64_t));
	if(!buf){
		return NULL;
	}
	// Initialize buf to be 0 for 0, 1 for 1, and 2 for every other n <= max.
	// 2 is not a valid output value, so we use it here to mark a number as prime
	// until they are marked off by being the factor of some other number.
	buf[0] = 0xAAAAAAAAAAAAAAA4ull;
	memset(buf + 1, 0xAA, (buf_len - 1)*sizeof(uint64_t));
	for(uint64_t n = 2; n <= max && n; ++n){
		if(((buf[n/32] >> (n%32*2)) & 3) != 2){
			continue;
		}
		buf[n/32] ^= 1ull << (n%32*2);
		uint64_t m;
		if(__builtin_mul_overflow(n, 2, &m)){
			continue;
		}
		// First we visit all multiples of n and flip their mobius value
		// We need to send 2 to 3, 1 to 3, 3 to 1, and 0 to 0.
		// We can do this by using an xor mask so we only have to
		// get x from the array and then xor the array with this mask.
		while(m <= max){
			// fprintf(stderr, "[%"PRIu64"]: %"PRIu64, );
			uint64_t x = (buf[m/32] >> (m%32*2)) & 3;
			x = (0x98 >> (x*2)) & 3;// this constant acts as a packed lookup table for the xor mask for x
			buf[m/32] ^= x << (m%32*2);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
		uint64_t n2;
		// Here we visit all multiples of n**2 and set their mobius value to 0,
		// since they are not square free.
		if(!__builtin_mul_overflow(n, n, &n2)){
			for(uint64_t m = n2; m <= max;){
				buf[m/32] &= ~(3ull << (m%32*2));
				if(__builtin_add_overflow(m, n2, &m)){
					break;
				}
			}
		}
	}
	return buf;
}

int64_t *compute_mertens_range(uint64_t max, const uint64_t mobius[static max/32 + 1]){
	int64_t *buf = malloc((max + 1)*sizeof(int64_t));
	if(!buf){
		return NULL;
	}
	buf[0] = 0;
	typedef struct{
		int64_t x:2;
	} int2_t;// used to signed extend 2 bit integers to 64 bit
	for(uint64_t n = 1; n <= max; ++n){
		int2_t mu = {bitfield2_arr_get(mobius, n)};
		buf[n] = buf[n - 1] + mu.x;
	}
	return buf;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverride-init"
static const uint64_t wheel_30_offsets[30] = {
	[0 ... 29] = 8, // error
	[1] = 0,
	[7] = 1,
	[11] = 2,
	[13] = 3,
	[17] = 4,
	[19] = 5,
	[23] = 6,
	[29] = 7
};
#pragma GCC diagnostic pop

static const int64_t wheel_30_lb_offsets[30] = {
	[0] = -1,
	[1 ... 6] = 0,
	[7 ... 10] = 1,
	[11 ... 12] = 2,
	[13 ... 16] = 3,
	[17 ... 18] = 4,
	[19 ... 22] = 5,
	[23 ... 28] = 6,
	[29] = 7
};

uint64_t *sieve_is_composite(uint64_t max){
	uint64_t is_composite_len = max/(16*15) + 1;
	uint64_t *is_composite = calloc(is_composite_len, sizeof(uint64_t));
	if(!is_composite){
		return NULL;
	}
	uint64_t num_primes = 0;
	uint64_t n = 7;
	for(uint64_t n_max = u64_nth_root(max, 2); n <= n_max; ++n){
		uint64_t i = n/30 + wheel_30_offsets[n%30];
		if(is_composite[i >> 6] & (1ull << (i&0x7F))){
			continue;
		}
		for(uint64_t m = n*n; m <= max;){
			uint64_t j = m/30 + wheel_30_offsets[m%30]; // TODO: try looping over residues mod 30 in outer loop and other things.
			is_composite[j >> 6] |= 1ull << (j&0x7F);
			if(__builtin_add_overflow(m, n, &m)){
				break;
			}
		}
	}
	return is_composite;
}

bool is_composite(uint64_t n, const uint64_t buf[static n/(16*15) + 1]){
	uint64_t offset = wheel_30_offsets[n%30];
	if(offset == 8){
		switch(n){
			case 2:
			case 3:
			case 5:
			return false;
		}
		return true;
	}
	uint64_t i = n/30 + offset;
	return buf[i >> 6] & (1ull << (i&0x7F));
}

uint64_t *compute_pi_range(uint64_t max, const uint64_t buf[static max/(16*15) + 1]){
	uint64_t res_len = max/(16*15) + 1;
	uint64_t *res = calloc(res_len, sizeof(uint64_t));
	if(!res){
		return NULL;
	}
	for(uint64_t i = 0; i < res_len; ++i){
		res[i] += __builtin_popcountll(~buf[i]);
	}
	return res;
}

uint64_t compute_pi_from_tables(uint64_t n, const uint64_t pi_table[static n/(16*15) + 1], const uint64_t buf[static n/(16*15) + 1]){
	int64_t i = (int64_t)(n/30) + wheel_30_lb_offsets[n%30];
	if(i < 64){
		
	}
	uint64_t mask = ~buf[i >> 6];
	mask &= ~(~0ull << (i&0x7));
	return pi_table[i >> 6] + __builtin_popcountll(mask);
}

uint64_t *sieve_primes(uint64_t max, uint64_t *_num_primes){
	uint64_t is_composite_len = max/64 + 1;
	uint64_t *is_composite = calloc(is_composite_len, sizeof(uint64_t));
	if(!is_composite){
		return NULL;
	}
	uint64_t *primes = malloc((size_t)max_primes_le(max)*sizeof(uint64_t));
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
	*_num_primes = num_primes;
	free(is_composite);
	return realloc(primes, num_primes*sizeof(uint64_t)) ?: primes;
}

typedef struct{
	uint64_t num_primes;
	uint64_t size;
	uint64_t num_residues;
	uint64_t *residues;
	uint64_t *inv_residues;
	uint64_t *residue_idxs;
	uint64_t *product_idxs;
	uint64_t *product_offs;
} wheel_info_t;

static inline bool init_wheel(wheel_info_t *self, uint64_t num_primes){
	if(num_primes > 15){
		return false;
	}
	uint64_t size = 1;
	uint64_t num_residues = 1;
	for(uint64_t i = 0; i < num_primes; ++i){
		size *= small_primes[i];
		num_residues *= small_primes[i] - 1;
	}
	if(!(self->residues = malloc(num_residues*sizeof(uint64_t)))){
		goto CLEANUP_0;
	}else if(!(self->inv_residues = malloc(num_residues*sizeof(uint64_t)))){
		goto CLEANUP_1;
	}else if(!(self->residue_idxs = malloc(size*sizeof(uint64_t)))){
		goto CLEANUP_2;
	}else if(!(self->product_idxs = malloc(num_residues*num_residues*sizeof(uint64_t)))){
		goto CLEANUP_3;
	}else if(!(self->product_offs = malloc(num_residues*num_residues*sizeof(uint64_t)))){
		goto CLEANUP_4;
	}
	self->num_primes = num_primes;
	self->size = size;
	self->num_residues = num_residues;
	for(uint64_t j = 0, r = 1; r < size; ++ r){
		int64_t r_inv;
		if(egcd(r, size, &r_inv, NULL) == 1){
			self->residues[j] = r;
			self->inv_residues[j] = r_inv < 0 ? r_inv + (int64_t)size : r_inv;
			++j;
		}
	}
	for(uint64_t i = 0, j = 0; i < size; ++i){
		self->residue_idxs[i] = i == self->residues[j] ? j++ : num_residues;
	}
	for(uint64_t i = 0; i < num_residues; ++i){
		for(uint64_t j = 0; j <= i; ++j){
			uint64_t rs = self->residues[i]*self->residues[j];
			uint64_t k = self->residue_idxs[rs%size];
			uint64_t offset = rs / size;
			self->product_idxs[i*num_residues + j] = k;
			self->product_offs[i*num_residues + j] = offset;
			if(i != j){
				self->product_idxs[j*num_residues + i] = k;
				self->product_offs[j*num_residues + i] = offset;
			}
		}
	}
	return true;
	CLEANUP_4:
	free(self->product_idxs);
	CLEANUP_3:
	free(self->residue_idxs);
	CLEANUP_2:
	free(self->inv_residues);
	CLEANUP_1:
	free(self->residues);
	CLEANUP_0:
	return false;
}

static inline void delete_wheel(wheel_info_t *self){
	free(self->residues);
	free(self->inv_residues);
	free(self->residue_idxs);
	free(self->product_idxs);
	free(self->product_offs);
}

typedef struct{
	uint64_t num_primes;
	uint64_t *buf;
} spoke_info_t;

static inline void delete_spokes(uint64_t num_residues, spoke_info_t spokes[static num_residues]){
	for(uint64_t i = 0; i < num_residues; ++i){
		free(spokes[i].buf);
	}
}

static inline spoke_info_t *init_spokes(const wheel_info_t *wheel_info, uint64_t spoke_prime_ub){
	uint64_t max_num_spoke_primes = 2*spoke_prime_ub/(wheel_info->num_residues*log((double)spoke_prime_ub/wheel_info->size));
	spoke_info_t *sieve_spokes = malloc(wheel_info->num_residues*sizeof(spoke_info_t));
	if(!sieve_spokes){
		return NULL;
	}
	{
		uint64_t i;
		for(i = 0; i < wheel_info->num_residues; ++i){
			sieve_spokes[i].num_primes = 0;
			sieve_spokes[i].buf = malloc((size_t)max_num_spoke_primes*sizeof(uint64_t));
			if(!sieve_spokes[i].buf){
				break;
			}
		}
		if(i != wheel_info->num_residues){
			delete_spokes(i, sieve_spokes);
			free(sieve_spokes);
			return NULL;
		}
	}
	return sieve_spokes;
}

typedef struct{
	uint64_t a, b;
	uint8_t *is_composite;
} bucket_info_t;

static inline void delete_buckets(uint64_t num_buckets, bucket_info_t buckets[static num_buckets]){
	for(uint64_t i = 0; i < num_buckets; ++i){
		free(buckets[i].is_composite);
	}
}

static inline bucket_info_t *init_buckets(uint64_t num_buckets, uint64_t bucket_size){
	bucket_info_t *buckets = malloc(num_buckets*sizeof(bucket_info_t));
	if(!buckets){
		return NULL;
	}
	{
		uint64_t i;
		for(i = 0; i < num_buckets; ++i){
			buckets[i].is_composite = malloc(bucket_size);
			if(!buckets[i].is_composite){
				break;
			}
		}
		if(i != num_buckets){
			delete_buckets(i, buckets);
			free(buckets);
			return NULL;
		}
	}
	return buckets;
}	

static inline uint64_t setup_sieving_primes(uint64_t *out, const wheel_info_t *wheel_info, spoke_info_t sieve_spokes[static wheel_info->num_residues], bucket_info_t *bucket){
	//printf("Generating sieving primes (%"PRIu64" to %"PRIu64")...\n", bucket->a, bucket->b);
	uint64_t bucket_num_wheels = bucket->b/wheel_info->size;
	uint64_t b_idx = wheel_info->residue_idxs[bucket->b%wheel_info->size];
	uint64_t bucket_num_residues = bucket_num_wheels*wheel_info->num_residues + b_idx + 1;
	uint64_t num_primes = 0;
	memset(bucket->is_composite, 0, bucket_num_residues);
	// (kw + r)(lw + s) = l(kw + r)w + ksw + rs <= B
	// first we sieve the first copy of the wheel
	//printf("First copy of the wheel:\n");
	for(uint64_t r_idx = 1; r_idx < wheel_info->num_residues; ++r_idx){
		if(bucket->is_composite[r_idx]){
			continue;
		}
		uint64_t r = wheel_info->residues[r_idx];
		//printf("    Marking off multiples of %"PRIu64"\n", r);
		out[num_primes++] = r;
		sieve_spokes[r_idx].buf[sieve_spokes[r_idx].num_primes++] = 0;
		for(uint64_t s_idx = 0; s_idx < wheel_info->num_residues; ++s_idx){
			//uint64_t s = wheel_info->residues[s_idx];
			uint64_t offset = wheel_info->product_offs[r_idx*wheel_info->num_residues + s_idx];
			uint64_t rs_idx = wheel_info->product_idxs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t coeff = offset + !!(rs_idx > b_idx);
				if(coeff >= bucket->b/wheel_info->size){
					break;
				}
			uint64_t l_lb = (s_idx < r_idx);
			uint64_t l_ub = (bucket->b/wheel_info->size - coeff)/r;
			for(uint64_t l = l_lb; l <= l_ub; ++l){
				bucket->is_composite[(l*r + offset)*wheel_info->num_residues + rs_idx] = true;
			}
		}
	}
	// now we sieve the remaining copies of the wheel.  TODO special case the last copy if it is incomplete
	uint64_t k;
	// this loop condition is currently that (k*w + 1)^2 <= bucket_max, but if we make it
	// (k*w + w - 1)^2 then we can eliminate some checks and only do them in an optimized loop
	// for the interveining k
	//printf("Remaining copies of the wheel with at least one multiple in the sieving prime range:\n");
	for(k = 1; (k*wheel_info->size + 2)*k*wheel_info->size + 1 <= bucket->b; ++k){
		for(uint64_t r_idx = 0; r_idx < wheel_info->num_residues; ++r_idx){
			if(bucket->is_composite[k*wheel_info->num_residues + r_idx]){
				continue;
			}
			uint64_t r = wheel_info->residues[r_idx];
			uint64_t n = k*wheel_info->size + r;
			//printf("    Marking off multiples of %"PRIu64"\n", n);
			out[num_primes++] = n;
			sieve_spokes[r_idx].buf[sieve_spokes[r_idx].num_primes++] = k;
			for(uint64_t s_idx = 0; s_idx < wheel_info->num_residues; ++s_idx){
				uint64_t s = wheel_info->residues[s_idx];
				uint64_t offset = wheel_info->product_offs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t rs_idx = wheel_info->product_idxs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t coeff = k*s + offset + !!(rs_idx > b_idx);
				if(coeff >= bucket->b/wheel_info->size){
					break;
				}
				uint64_t l_lb = k + (s_idx < r_idx);
				uint64_t l_ub = (bucket->b/wheel_info->size - coeff)/n;
				for(uint64_t l = l_lb; l <= l_ub; ++l){
					bucket->is_composite[(l*n + k*s + offset)*wheel_info->num_residues + rs_idx] = true;
				}
			}
		}
	}
	//printf("Finding numbers still not marked as composite in the sieving prime range:\n");
	for(; k <= (bucket->b + 1)/wheel_info->size - 1; ++k){
		for(uint64_t r_idx = 0; r_idx < wheel_info->num_residues; ++r_idx){
			if(!bucket->is_composite[k*wheel_info->num_residues + r_idx]){
				uint64_t r = wheel_info->residues[r_idx];
				uint64_t n = k*wheel_info->size + r;
				//printf("    found %"PRIu64"\n", n);
				out[num_primes++] = n;
				sieve_spokes[r_idx].buf[sieve_spokes[r_idx].num_primes++] = k;
			}
		}
	}
	if((bucket->b + 1)%wheel_info->size == 0){
		return num_primes;
	}
	// primes between floor(B/w)*w and B will be found by a sieve_one_bucket step because we want these to always start on residue 0
	// so even though we may have found some, we don't add them to the output here, just the spokes
	for(uint64_t r_idx = 0; r_idx <= b_idx; ++r_idx){
		if(!bucket->is_composite[k*wheel_info->num_residues + r_idx]){
			//uint64_t r = wheel_info->residues[r_idx];
			//uint64_t n = k*wheel_info->size + r;
			sieve_spokes[r_idx].buf[sieve_spokes[r_idx].num_primes++] = k;
		}
	}
	return num_primes;
}

static inline void sieve_one_bucket(const wheel_info_t *wheel_info, const spoke_info_t sieve_spokes[static wheel_info->num_residues], bucket_info_t *bucket, bool is_full){
	// full buckets are characterized by bucket->a being the first possible residue (1) and bucket->b being the last possible residue (w - 1) (mod w), so
	// to find the number of wheels spanned from bucket->a to bucket->b we have to add 2 to get to a multiple of w (bucket size)
	uint64_t bucket_num_wheels;
	uint64_t bucket_num_residues;
	uint64_t b_idx;
	if(is_full){
		bucket_num_wheels = (bucket->b - bucket->a + 2)/wheel_info->size;
		bucket_num_residues = bucket_num_wheels*wheel_info->num_residues;
		b_idx = wheel_info->num_residues - 1;
	}else{
		bucket_num_wheels = (bucket->b - bucket->a + 1)/wheel_info->size;
		b_idx = wheel_info->residue_idxs[bucket->b%wheel_info->size];
		bucket_num_residues = bucket_num_wheels*wheel_info->num_residues + b_idx + 1;
	}
	uint64_t l0 = bucket->a/wheel_info->size;
	//printf(
	//	"Sieving %s bucket from %"PRIu64" to %"PRIu64
	//	" with %"PRIu64" copies of the wheel, %"PRIu64" residues of the wheel, and minimal k=%"PRIu64"...\n",
	//	is_full ? "a full" : "a non-full (final)", bucket->a, bucket->b, bucket_num_wheels, bucket_num_residues, l0);
	memset(bucket->is_composite, 0, bucket_num_residues);
	// A <= (l*w + s)*(k*w + r) <= B
	// For what k*w + r will all values of s be permissible for at least one l?
	// for full buckets, the following holds.  for non-full buckets, it is overly restrictive, but that's ok
	// A <= (l*w + 1)*(k*w + r), (l*w + w - 1)*(k*w + r) <= B
	// A/w <= k*l*w + k + l*r, k*(l+1)*w - k + (l+1)*r <= B/w
	// A/w + k*(l+1)*w - k + (l+1)*r <= B/w + k*l*w + k + l*r
	// k*w + r <= (B - A)/w + 2*k
	// k*(w - 2) <= (B - A)/w - r
	for(uint64_t r_idx = 0; r_idx < wheel_info->num_residues; ++r_idx){
		if(!sieve_spokes[r_idx].num_primes){
			continue;
		}
		uint64_t r = wheel_info->residues[r_idx];
		uint64_t k_ub = 0;
		if(is_full || bucket_num_wheels > r){
			k_ub = (bucket_num_wheels - r)/(wheel_info->size - 2);
		}
		//printf("    For residue %"PRIu64", the maximal k allowing all s is %"PRIu64"\n    Processing small k...\n", r, k_ub);
		uint64_t k_idx, k;
		// for this first loop, we only process k which are small enough so that all residues s will have some permissible l
		for(k_idx = 0; k_idx < sieve_spokes[r_idx].num_primes; ++k_idx){
			k = sieve_spokes[r_idx].buf[k_idx];
			if(k > k_ub){
				break;
			}
			uint64_t n = k*wheel_info->size + r;// TODO: add this to spoke_info
			//printf("        prime: %"PRIu64"\n", n);
			// TODO: we do two divisions inside this loop to find the bounds for l (l_lb and l_ub)
			// For small primes n this is faster than also keeping track of the actual values of the multiples,
			// but we should determine emperically how large n should be to make simple keeping track of the multiples
			// faster and then switch this loop.
			for(uint64_t s_idx = 0; s_idx < wheel_info->num_residues; ++s_idx){
				uint64_t s = wheel_info->residues[s_idx];
				uint64_t offset = wheel_info->product_offs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t rs_idx = wheel_info->product_idxs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t coeff;
				if(is_full || rs_idx <= b_idx){
					coeff = k*s + offset;
				}else{
					coeff = k*s + offset + 1;
				}
				if(coeff >= bucket->b/wheel_info->size){
					break;
				}
				// A <= (l*w + s)*(k*w + r)
				// A <= (k*l*w + k*s + l*r + offset(r*s))*w [+ residue(r*s)]
				// ((A/w - offset(r*s) - k*s) + n - 1)/n <= l
				uint64_t l_lb = (l0 - offset - k*s + n - 1)/n;
				// (l*w + s)*(k*w + r) <= B
				// (k*l*w + k*s + l*r + offset(r*s))*w [+ residue(r*s)] <= B
				// l <= (B/w - offset(r*s) - k*s)/n
				uint64_t l_ub = (bucket->b/wheel_info->size - coeff)/n;
				//printf("        %"PRIu64": [%"PRIu64", %"PRIu64"]\n", s, l_lb, l_ub);
				for(uint64_t l = l_lb; l <= l_ub; ++l){
					bucket->is_composite[(k*l*wheel_info->size + k*s + l*r + offset - l0)*wheel_info->num_residues + rs_idx] = true;
				}
			}
		}
		// now we have to process large k where not all s might have permissible l.  If the bucket size is at least the square root
		// of the max prime up to the square root of the sieve, there will always be at least s with a permissible l
		// Also, if we want to allow buckets smaller than this (ie the final bucket), sometimes some of the primes in sieve_spokes
		// won't have multiples in a given bucket, so we will need a special case for these :(
		//printf("    Processing large k >=%"PRIu64"...\n", k);
		while(1){
			uint64_t n = k*wheel_info->size + r;
			//printf("        prime: %"PRIu64"\n", n);
			for(uint64_t s_idx = 0; s_idx < wheel_info->num_residues; ++s_idx){
				uint64_t s = wheel_info->residues[s_idx];
				uint64_t offset = wheel_info->product_offs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t rs_idx = wheel_info->product_idxs[r_idx*wheel_info->num_residues + s_idx];
				uint64_t coeff;
				if(is_full || rs_idx <= b_idx){
					coeff = k*s + offset;
				}else{
					coeff = k*s + offset + 1;
				}
				if(coeff >= bucket->b/wheel_info->size){
					break;
				}
				uint64_t l_lb = (l0 - offset - k*s + n - 1)/n;
				uint64_t l_ub = (bucket->b/wheel_info->size - coeff)/n;
				//printf("        %"PRIu64": [%"PRIu64", %"PRIu64"]\n", s, l_lb, l_ub);
				for(uint64_t l = l_lb; l <= l_ub; ++l){
					bucket->is_composite[(k*l*wheel_info->size + k*s + l*r + offset - l0)*wheel_info->num_residues + rs_idx] = true;
				}
			}
			if(++k_idx >= sieve_spokes[r_idx].num_primes){
				break;
			}
			k = sieve_spokes[r_idx].buf[k_idx];
			// TODO: if we want to support small buckets, we need to stop and handle primes with one or fewer multiples per bucket
			// note that n is computed at the top of the loop; it's rolled a little weirdly due to having k_idx and k from the above loop
		}
	}
}

static inline uint64_t collect_bucket_primes(uint64_t *out, const wheel_info_t *wheel_info, const bucket_info_t *bucket, bool is_full){
	uint64_t num_primes = 0;
	uint64_t b_idx = wheel_info->residue_idxs[bucket->b%wheel_info->size];
	uint64_t l0 = bucket->a/wheel_info->size;
	uint64_t k, k_ub = (bucket->b - wheel_info->size + 1)/wheel_info->size;
	// TODO: try reading all the indices and only converting them to numbers if they are not marked as composite
	for(k = l0; k <= k_ub; ++k){
		for(uint64_t r_idx = 0; r_idx < wheel_info->num_residues; ++r_idx){
			uint64_t r = wheel_info->residues[r_idx];
			if(!bucket->is_composite[(k - l0)*wheel_info->num_residues + r_idx]){
				uint64_t n = k*wheel_info->size + r;
				out[num_primes++] = n;
				//printf("%"PRIu64"\n", n);
			}
		}
	}
	if(is_full){
		return num_primes;
	}
	for(uint64_t r_idx = 0; r_idx <= b_idx; ++r_idx){
		uint64_t r = wheel_info->residues[r_idx];
		if(!bucket->is_composite[(k - l0)*wheel_info->num_residues + r_idx]){
			uint64_t n = k*wheel_info->size + r;
			out[num_primes++] = n;
			//printf("%"PRIu64"\n", n);
		}
	}
	return num_primes;
}

static inline uint64_t candidate_le(uint64_t n, const wheel_info_t *wheel_info){
	uint64_t r = n % wheel_info->size;
	if(!r){
		return n - 1;
	}
	while(wheel_info->residue_idxs[r--] == wheel_info->num_residues){
		n--;
	}
	return n;
}

uint64_t *sieve_primes_wheel(uint64_t max, uint64_t *_num_primes){
	uint64_t cache_size_bytes = (uint64_t)sysconf(_SC_LEVEL1_DCACHE_SIZE);
	uint64_t num_cpus = (uint64_t)sysconf(_SC_NPROCESSORS_ONLN);
	wheel_info_t wheel_info;
	uint64_t spoke_prime_ub = sqrtl(max);
	uint64_t *primes = NULL;
	if(!init_wheel(&wheel_info, max_prime_divs(sqrt(spoke_prime_ub)))){
		goto CLEANUP_0;
	}
	if(wheel_info.size < 6){
		fprintf(stderr, "\e[1;31mERROR: segmeneted sieve cannot currently operate for max < 1296: wheel sizes < 6 not supported!\e[0m\n");
		goto CLEANUP_1;
	}
	spoke_prime_ub = candidate_le(spoke_prime_ub, &wheel_info);
	max = candidate_le(max, &wheel_info);
	printf(
		"Sieving primes up to %"PRIu64"\n"
		"Min divisor upper bound: %"PRIu64"\n"
		"Wheel: %"PRIu64" (%"PRIu64" residues)\n"
		"Cache size: %.1fKB\n"
		"Processors available: %"PRIu64"\n",
		max, spoke_prime_ub, wheel_info.size, wheel_info.num_residues, cache_size_bytes/1024., num_cpus);
	spoke_info_t *sieve_spokes = init_spokes(&wheel_info, spoke_prime_ub);
	if(!sieve_spokes){
		goto CLEANUP_1;
	}
	uint64_t num_buckets = num_cpus;
	uint64_t bucket_size = wheel_info.size*wheel_info.size;
	uint64_t bucket_num_residues;
	if(bucket_size >= spoke_prime_ub){
		bucket_num_residues = wheel_info.size*wheel_info.num_residues;
	}else{
		bucket_num_residues = spoke_prime_ub/wheel_info.size*wheel_info.num_residues + wheel_info.residue_idxs[spoke_prime_ub%wheel_info.size] + 1;
	}
	bucket_info_t *buckets = init_buckets(num_buckets, bucket_num_residues);
	if(!buckets){
		goto CLEANUP_2;
	}
	primes = malloc((size_t)max_primes_le(max)*sizeof(uint64_t));
	if(!primes){
		goto CLEANUP_3;
	}
	uint64_t num_primes = wheel_info.num_primes;
	memcpy(primes, small_primes, num_primes*sizeof(uint64_t));
	buckets->a = wheel_info.residues[1];
	buckets->b = spoke_prime_ub;
	num_primes += setup_sieving_primes(primes + num_primes, &wheel_info, sieve_spokes, buckets);
	buckets->a = buckets->b/wheel_info.size*wheel_info.size + 1;
	buckets->b = buckets->a + bucket_size - 2;
	while(1){
		if(buckets->b <= max){
			sieve_one_bucket(&wheel_info, sieve_spokes, buckets, true);
			num_primes += collect_bucket_primes(primes + num_primes, &wheel_info, buckets, true);
		}else if(buckets->a <= max){
			buckets->b = max;
			sieve_one_bucket(&wheel_info, sieve_spokes, buckets, false);
			num_primes += collect_bucket_primes(primes + num_primes, &wheel_info, buckets, false);
			break;
		}else{
			break;
		}
		buckets->a += bucket_size;
		buckets->b += bucket_size;
	}
	*_num_primes = num_primes;
	CLEANUP_3:
	delete_buckets(num_buckets, buckets);
	free(buckets);
	CLEANUP_2:
	delete_spokes(wheel_info.num_residues, sieve_spokes);
	free(sieve_spokes);
	CLEANUP_1:
	delete_wheel(&wheel_info);
	CLEANUP_0:
	return primes ? realloc(primes, num_primes*sizeof(uint64_t)) ?: primes : NULL;
}

