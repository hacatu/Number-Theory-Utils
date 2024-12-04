#include <stddef.h>
#if __has_include(<linux/version.h>)
#include <linux/version.h>
#if LINUX_VERSION_CODE >= KERNEL_VERSION(3,17,0)
#if __has_include(<gnu/libc-version.h>)
#include <gnu/libc-version.h>
#include <sys/random.h>
#if (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 25) || __GLIBC__ > 2
#define _NUT_RAND_USE_GETRANDOM
#else
#include <sys/syscall.h>
#include <unistd.h>
#define _NUT_RAND_USE_SYSCALL
#endif
#else
#include <stdio.h>
#define _NUT_RAND_USE_URANDOM
#endif
#endif
#elifdef __MINGW32__
#define _CRT_RAND_S
#include <stdlib.h>
#define _NUT_RAND_USE_RAND_S
#endif
#include <string.h>
#include <nut/debug.h>
#include <nut/modular_math.h>

uint64_t nut_u64_pow(uint64_t b, uint64_t e){
	if(!e){
		return 1;
	}
	uint64_t r = 1;
	while(1){
		if(e&1){
			r = r*b;
		}
		if(!(e >>= 1)){
			return r;
		}
		b *= b;
	}
}

uint128_t nut_u128_pow(uint128_t b, uint64_t e){
	if(!e){
		return 1;
	}
	uint128_t r = 1;
	while(1){
		if(e&1){
			r = r*b;
		}
		if(!(e >>= 1)){
			return r;
		}
		b *= b;
	}
}

uint64_t nut_u64_powmod(uint64_t b, uint64_t e, uint64_t n){
	if(!e){
		return 1;
	}
	uint64_t r = 1;
	b %= n;
	while(1){
		if(e&1){
			r = (uint128_t)r*b%n;
		}
		if(!(e >>= 1)){
			return r;
		}
		b = (uint128_t)b*b%n;
	}
}

uint64_t nut_u64_binom(uint64_t n, uint64_t k){
	uint64_t res = 1;
	if(n - k < k){
		k = n - k;
	}
	for(uint64_t i = 0; i < k; ++i){
		res = res*(n - i)/(1 + i);
	}
	return res;
}

uint64_t nut_u64_binom_next(uint64_t n, uint64_t k, uint64_t prev){
	return prev*(n - k + 1)/k;
}

#if defined(_NUT_RAND_USE_GETRANDOM) || defined(_NUT_RAND_USE_SYSCALL) || defined(_NUT_RAND_USE_URANDOM)
uint64_t nut_u64_rand(uint64_t a, uint64_t b){
#ifdef _NUT_RAND_USE_URANDOM
	__thread FILE *urandom = NULL;
	if(!urandom){
		fopen("/dev/urandom", "rb");
	}
#endif
	uint64_t l = b - a, r = 0, bytes = (71 - __builtin_clzll(l))/8;
	uint64_t ub;
	if(bytes == 8){
		ub = ~0ull%l + 1;
		ub = (ub == l) ? 0 : -ub;
	}else{
		ub = 1ull << (bytes*8);
		ub -= ub%l;
	}
	do{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#ifdef _NUT_RAND_USE_GETRANDOM
		getrandom(&r, bytes, 0);
#elifdef _NUT_RAND_USE_SYSCALL
		syscall(SYS_getrandom, &r, bytes, 0);
#else
		fread(&r, 8, 1, urandom);
#endif
#pragma GCC diagnostic pop
	}while(ub && r >= ub);
	return r%l + a;
}
#endif

#if defined(_NUT_RAND_USE_RAND_S)
uint64_t nut_u64_rand(uint64_t a, uint64_t b){
	uint64_t l = b - a, r = 0;
	uint64_t ub;
	if(l > (1ull << 32)){
		ub = ~0ull%l + 1;
		ub = (ub == l) ? 0 : -ub;
		do{
			rand_s((uint32_t*)&r);
			rand_s((uint32_t*)&r + 1);
		}while(ub && r >= ub);
	}else{
		ub = 1ull << 32;
		ub -= ub%l;
		do{
			rand_s((uint32_t*)&r);
		}while(ub && r >= ub);
	}
	return r%l + a;
}
#endif

uint64_t nut_u64_prand(uint64_t a, uint64_t b){
	return nut_u64_rand(a, b);//TODO: test if this is a bottleneck, we don't need calls to this function to be secure
}

int64_t nut_i64_egcd(int64_t a, int64_t b, int64_t *restrict _t, int64_t *restrict _s){
	int64_t r0 = b, r1 = a;
	int64_t s0 = 1, s1 = 0;
	int64_t t0 = 0, t1 = 1;
	while(r1){
		int64_t q = r0/r1, t;
		t = r1;
		r1 = r0 - q*r1;
		r0 = t;
		t = s1;
		s1 = s0 - q*s1;
		s0 = t;
		t = t1;
		t1 = t0 - q*t1;
		t0 = t;
	}
	if(_t){
		*_t = t0;
	}
	if(_s){
		*_s = s0;
	}
	return r0;
}

int64_t nut_i64_modinv(int64_t a, int64_t b){
	int64_t ainv;
	nut_i64_egcd(a, b, &ainv, NULL);
	return ainv < 0 ? ainv + b : ainv;
}

// This uses a hensel/newton like iterative algorithm, described here
// https://crypto.stackexchange.com/a/47496
// Basically, we use a lookup table to get the inverse mod 2**8, and then
// use the fact that ax = 1 mod 2**k --> ax(2-ax) = 1 mod 2**(2k) to lift the inverse to
// mod 2**16, mode 2**32, etc as needed.  This does cap out at 2**64.
NUT_ATTR_NO_SAN("unsigned-shift-base")
NUT_ATTR_NO_SAN("unsigned-integer-overflow")
uint64_t nut_u64_modinv_2t(uint64_t a, uint64_t t){
	static const uint8_t modinv_256_tbl[] = {
		1, 171, 205, 183, 57, 163, 197, 239, 241, 27, 61, 167, 41, 19, 53, 223,
		225, 139, 173, 151, 25, 131, 165, 207, 209, 251, 29, 135, 9, 243, 21, 191,
		193, 107, 141, 119, 249, 99, 133, 175, 177, 219, 253, 103, 233, 211, 245, 159,
		161, 75, 109, 87, 217, 67, 101, 143, 145, 187, 221, 71, 201, 179, 213, 127,
		129, 43, 77, 55, 185, 35, 69, 111, 113, 155, 189, 39, 169, 147, 181, 95,
		97, 11, 45, 23, 153, 3, 37, 79, 81, 123, 157, 7, 137, 115, 149, 63,
		65, 235, 13, 247, 121, 227, 5, 47, 49, 91, 125, 231, 105, 83, 117, 31,
		33, 203, 237, 215, 89, 195, 229, 15, 17, 59, 93, 199, 73, 51, 85, 255
	};
	assert(t <= 64 && a&1);
	uint64_t t1 = 8;
	uint64_t x = modinv_256_tbl[(a >> 1)&0x7F];
	while(t1 < t){
		t1 *= 2;
		x = (x*(2 - a*x));
		if(t1 != 64){
			x &= (1ull << t1) - 1;
		}
	}
	return t < 64 ? x&((1ull << t) - 1) : x;
}

int64_t nut_i64_mod(int64_t a, int64_t n){
	int64_t r = a%n;
	if(r < 0){
		r += n;
	}
	return r;
}

int64_t nut_i64_crt(int64_t a, int64_t p, int64_t b, int64_t q){
	int64_t x, y;
	nut_i64_egcd(p, q, &x, &y);
	return nut_i64_mod(b*p%(p*q)*x + a*q%(p*q)*y, p*q);
}

int128_t nut_i128_crt(int64_t a, int64_t p, int64_t b, int64_t q){
	int64_t _x, _y;
	nut_i64_egcd(p, q, &_x, &_y);
	int128_t x = _x, y = _y;
	x = b*p%(p*q)*x + a*q%(p*q)*y;
	x %= p*q;
	return x < 0 ? x + p*q : x;
}

int64_t nut_i64_lcm(int64_t a, int64_t b){
	return a*b/nut_i64_egcd(a, b, NULL, NULL);
}

NUT_ATTR_NO_SAN("unsigned-shift-base")
NUT_ATTR_NO_SAN("unsigned-integer-overflow")
uint64_t nut_u64_binom_next_mod_2t(uint64_t n, uint64_t k, uint64_t t, uint64_t *restrict v2, uint64_t *restrict p2){
	uint64_t num = n - k + 1;
	uint64_t num_v2 = __builtin_ctz(num);
	uint64_t denom_v2 = __builtin_ctz(k);
	*v2 = *v2 + num_v2 - denom_v2;
	uint64_t denom_pinv = nut_u64_modinv_2t(k >> denom_v2, t);
	uint64_t mask = t < 64 ? ((1ull << t) - 1) : ~0ull;
	*p2 = (*p2 * (num >> num_v2) * denom_pinv) & mask;
	return (*p2 << *v2) & mask;
}

int64_t nut_i64_jacobi(int64_t n, int64_t k){
	if(n%k == 0){
		return 0;
	}
	int64_t j = 1;
	while(1){
		int64_t s = __builtin_ctzll(n);
		int64_t q = n >> s;
		if((s&1) && ((k&7) == 3 || ((k&7) == 5))){
			j = -j;
		}
		if(q == 1){
			return j;
		}else if(q == k - 1){
			return (k&3)==1 ? j : -j;
		}else if((q&2) && (k&2)){
			j = -j;
		}
		n = k%q;
		k = q;
	}
}

uint64_t *nut_u64_make_jacobi_tbl(uint64_t p, int64_t partial_sums[restrict static p]){
	uint64_t *is_qr = calloc((p + 63)/64, sizeof(uint64_t));
	if(!is_qr){
		return NULL;
	}
	for(uint64_t n = 1, nn = 1; n <= p/2; ++n){
		is_qr[nn/64] |= 1ull << (nn%64);
		nn += 2*n + 1;
		if(nn >= p){
			nn -= p;
		}
	}
	partial_sums[0] = 0;
	int64_t acc = 0;
	for(uint64_t n = 1; n < p; ++n){
		bool qr = is_qr[n/64] & (1ull << (n%64));
		acc += qr ? 1 : -1;
		partial_sums[n] = acc;
	}
	return is_qr;
}

int64_t nut_u64_jacobi_tbl_get(uint64_t n, uint64_t p, const uint64_t is_qr[restrict static (p + 63)/64]){
	uint64_t r = n%p;
	if(!r){
		return 0;
	}
	return (is_qr[r/64] & (1ull << (r%64))) ? 1 : -1;
}

int64_t nut_i64_rand_nr(int64_t p){
	while(1){
		int64_t z = nut_u64_rand(2, p);
		if(nut_i64_jacobi(z, p) == -1){
			return z;
		}
	}
}

int64_t nut_i64_sqrt_shanks(int64_t n, int64_t p){
	int64_t s = __builtin_ctzll(p-1);
	int64_t q = p >> s;//p-1 = q*2^s
	int64_t z = nut_i64_rand_nr(p);
	//printf("trying \"nonresidue\" %"PRIu64"\n", z);
	int64_t m = s;
	int64_t c = nut_u64_powmod(z, q, p);
	int64_t t = nut_u64_powmod(n, q, p);
	int64_t r = nut_u64_powmod(n, (q + 1) >> 1, p);
	while(t != 1){
		int64_t i = 1;
		for(int64_t s = (int128_t)t*t%p; s != 1; s = (int128_t)s*s%p, ++i);
		int64_t b = c;
		for(int64_t j = 0; j < m - i - 1; ++j){
			b = (int128_t)b*b%p;
		}
		m = i;
		c = (int128_t)b*b%p;
		t = (int128_t)t*c%p;
		r = (int128_t)r*b%p;
	}
	return r;
}

int64_t nut_i64_sqrt_cipolla(int64_t n, int64_t p){
	int64_t a, w;
	do{
		a = nut_u64_rand(2, p);
		w = nut_i64_mod((int128_t)a*a%p - n, p);
	}while(nut_i64_jacobi(w, p) != -1);
	int64_t u_s = a, w_s = 1;
	int64_t u_r = 1, w_r = 0;
	for(int64_t k = (p + 1) >> 1; k; k >>= 1){
		if(k&1){
			int64_t _w_r = (int128_t)u_r*w_s%p;
			_w_r = (_w_r + (int128_t)w_r*u_s)%p;
			u_r = (int128_t)u_r*u_s%p;
			w_r = (int128_t)w_r*w_s%p;
			w_r = (int128_t)w_r*w%p;
			u_r = ((int128_t)u_r + w_r)%p;
			w_r = _w_r;
		}
		int64_t _w_s = (int128_t)2*u_s*w_s%p;
		u_s = (int128_t)u_s*u_s%p;
		w_s = (int128_t)w_s*w_s%p;
		w_s = (int128_t)w*w_s%p;
		u_s = ((int128_t)u_s + w_s)%p;
		w_s = _w_s;
	}
	return u_r;
}

int64_t nut_i64_sqrt_mod(int64_t n, int64_t p){
	int64_t r;
	if((p&3) == 3){
		r = nut_u64_powmod(n, (p + 1) >> 2, p);
	}else if((p&7) == 5){
		r = nut_u64_powmod(n, (p + 3) >> 3, p);
		if(r*r%p != n){
			r = (int128_t)r*nut_u64_powmod(2, (p - 1) >> 2, p)%p;
		}//can add 9 mod 16 case
	}else{
		int64_t m = 64 - __builtin_clzll(p);
		int64_t s = __builtin_ctzll(p - 1);
		if(8*m + 20 >= s*(s - 1)){
			r = nut_i64_sqrt_shanks(n, p);
		}else{
			r = nut_i64_sqrt_cipolla(n, p);
		}
	}
	return r;
}


/// See https://arxiv.org/pdf/1902.01961.pdf
uint64_t nut_i32_fastmod_init(uint32_t pd){
	return ~0ull/pd + 1 + (__builtin_popcount(pd) == 1);
}

uint64_t nut_u32_fastmod_init(uint32_t d){
	return ~0ull/d + 1;
}

uint128_t nut_i64_fastmod_init(uint64_t pd){
	return ~(uint128_t)0/pd + 1 + (__builtin_popcount(pd) == 1);
}

uint128_t nut_u64_fastmod_init(uint64_t d){
	return ~(uint128_t)0/d + 1;
}


int32_t nut_i32_fastmod_trunc(int32_t n, uint32_t pd, uint64_t c){
	uint64_t cn = c*n;
	int32_t cnd = ((uint128_t)cn*pd) >> 64;
	return cnd - ((pd - 1) & (n >> 31));
}

int32_t nut_i32_fastmod_floor(int32_t n, uint32_t pd, uint64_t c){
	uint64_t cn = c*n;
	int32_t cnd = ((uint128_t)cn*pd) >> 64;
	return cnd - ((pd - 1) && (n >> 31));
}

uint32_t nut_u32_fastmod(uint32_t n, uint32_t d, uint64_t c){
	uint64_t cn = c*n;
	return ((uint128_t)cn*d) >> 64;
}

int64_t nut_i64_fastmod_trunc(int64_t n, uint64_t pd, uint128_t c){
	uint128_t cn = c*n;
	uint64_t cn_hi = cn >> 64;
	int64_t cnd_hi = ((uint128_t)cn_hi*pd) >> 64;
	// c*n*pd >> 128 does not actually depend on the low bits of cn
	return cnd_hi - ((pd - 1) & (n >> 63));
}

int64_t nut_i64_fastmod_floor(int64_t n, uint64_t pd, uint128_t c){
	uint128_t cn = c*n;
	uint64_t cn_hi = cn >> 64;
	int64_t cnd_hi = ((uint128_t)cn_hi*pd) >> 64;
	// c*n*pd >> 128 does not actually depend on the low bits of cn
	return cnd_hi - ((pd - 1) && (n >> 63));
}

uint64_t nut_u64_fastmod(uint64_t n, uint64_t d, uint128_t c){
	uint128_t cn = c*n;
	uint64_t cn_hi = cn >> 64;
	return ((uint128_t)cn_hi*d) >> 64;
}


