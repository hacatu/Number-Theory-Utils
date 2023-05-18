#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include <nut/factorization.h>
#include <nut/dirichlet.h>
#include <nut/sieves.h>

uint64_t nut_dirichlet_D(uint64_t max){
	uint64_t y = nut_u64_nth_root(max, 2);
	uint64_t res = 0;
	for(uint64_t n = 1; n <= y; ++n){
		res += max/n;
	}
	return 2*res - y*y;
}

static inline bool is_composite_unpacked(uint64_t n, const uint8_t buf[static n/8 + 1]){
	return buf[n/8] & (1 << (n%8));
}

static inline void mark_composite_unpacked(uint64_t n, uint8_t buf[static n/8 + 1]){
	buf[n/8] |= (1 << (n%8));
}

// Implemented based on (https://gbroxey.github.io/blog/2023/04/30/mult-sum-1.html)
// in turn based on (https://codeforces.com/blog/entry/54090)
// which is just euler's sieve (https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes#Euler's_sieve)
// this is a modification of the sieve of eratosthenes which marks off each number only once.
// this makes strapping computation of an arbitrary multiplicative function on top of it easier,
// but even though the codeforces post claims it is linear, in practice it will often perform worse than
// a basic sieve of eratosthenes.  Unfortunately, while there are other ways to extend the sieve of eratosthenes
// to compute arbitrary multiplicative functions, they are nontrivial.  The approach I would like
// to try is a sieve of eratosthenes which keeps iterators to higher powers of the current prime p it is
// looking at ... basically, the normal sieve of eratosthenes doesn't have to worry about powers of p at all,
// and many specific function sieves are well behaved with respect to ascending prime powers, letting us just
// apply a quick modification to all the multiples of p**2 that just sits on top of the one we applied for multiples
// of p.  But in general, we need to make sure that we only apply the value to each output entry for the maximal
// power of each prime.
bool nut_euler_sieve_conv_u(int64_t n, const int64_t f_vals[static n+1], int64_t f_conv_u_vals[static n+1]){
	int64_t *smallest_ppow = malloc((n+1)*sizeof(int64_t));
	uint8_t *is_c_buf = calloc(n/8 + 1, sizeof(uint8_t));
	uint64_t num_primes = 0;
	int64_t *primes = malloc(nut_max_primes_le(n)*sizeof(int64_t));
	if(!smallest_ppow || !is_c_buf || !primes){
		free(smallest_ppow);
		free(is_c_buf);
		free(primes);
		return false;
	}
	f_conv_u_vals[1] = 1;
	for(int64_t i = 2; i <= n; ++i){
		// if i is prime, add it to the list of primes
		if(!is_composite_unpacked(i, is_c_buf)){
			primes[num_primes++] = i;
			f_conv_u_vals[i] = f_vals[i] + 1;
			smallest_ppow[i] = i;
		}
		// the key feature of euler's sieve is that it ONLY marks off multiples
		// of each number i times primes <= i.  This requires storing an array of
		// primes and reading it A LOT of times, but ensures that each composite number
		// is visited only once
		for(uint64_t j = 0; j < num_primes; ++j){
			int64_t p = primes[j], m;
			if(__builtin_mul_overflow(p, i, &m) || m > n){
				break;
			}
			mark_composite_unpacked(m, is_c_buf);
			if(i%p == 0){// i is a multiple of p, so in particular i and p are not coprime
				// also notice that because we break out of this loop at the bottom of this if statement,
				// p is in fact the smallest prime divisor of i
				int64_t ppow = smallest_ppow[i];
				int64_t B = ppow*p;
				smallest_ppow[m] = B;
				int64_t v = i/ppow;
				if(v != 1){// i is not a perfect power of p, that is, i = v*p**a with v != 1
					f_conv_u_vals[m] = f_conv_u_vals[v]*f_conv_u_vals[B];
				}else{// i is a perfect power of p (and so is m)
					// from the definition, (f <*> u)(p**a) = f(1)*u(p**a) + f(p)*u(p**(a-1)) + ... + f(p**a)*u(1)
					// = f(1) + f(p) + ... + f(p**a) = (f <*> u)(p**(a-1)) + f(p**a)
					f_conv_u_vals[m] = f_conv_u_vals[i] + f_vals[m];
				}
				break;
			}// otherwise, i and p are coprime, so (f <*> u)(i*p) = (f <*> u)(i)*(f <*> u)(p)
			f_conv_u_vals[m] = f_conv_u_vals[i]*f_conv_u_vals[p];
			smallest_ppow[m] = p;
		}
	}
	free(smallest_ppow);
	free(primes);
	free(is_c_buf);
	return true;
}

bool nut_euler_sieve_conv_N(int64_t n, const int64_t f_vals[static n+1], int64_t f_conv_N_vals[static n+1]){
	int64_t *smallest_ppow = malloc((n+1)*sizeof(int64_t));
	uint8_t *is_c_buf = calloc(n/8 + 1, sizeof(uint8_t));
	uint64_t num_primes = 0;
	int64_t *primes = malloc(nut_max_primes_le(n)*sizeof(int64_t));
	if(!smallest_ppow || !is_c_buf || !primes){
		free(smallest_ppow);
		free(is_c_buf);
		free(primes);
		return false;
	}
	f_conv_N_vals[1] = 1;
	for(int64_t i = 2; i <= n; ++i){
		// if i is prime, add it to the list of primes
		if(!is_composite_unpacked(i, is_c_buf)){
			primes[num_primes++] = i;
			f_conv_N_vals[i] = f_vals[i] + i;
			smallest_ppow[i] = i;
		}
		// the key feature of euler's sieve is that it ONLY marks off multiples
		// of each number i times primes <= i.  This requires storing an array of
		// primes and reading it A LOT of times, but ensures that each composite number
		// is visited only once
		for(uint64_t j = 0; j < num_primes; ++j){
			int64_t p = primes[j], m;
			if(__builtin_mul_overflow(p, i, &m) || m > n){
				break;
			}
			mark_composite_unpacked(m, is_c_buf);
			if(i%p == 0){// i is a multiple of p, so in particular i and p are not coprime
				// also notice that because we break out of this loop at the bottom of this if statement,
				// p is in fact the smallest prime divisor of i
				int64_t ppow = smallest_ppow[i];
				int64_t B = ppow*p;
				smallest_ppow[m] = B;
				int64_t v = i/ppow;
				if(v != 1){// i is not a perfect power of p, that is, i = v*p**a with v != 1
					f_conv_N_vals[m] = f_conv_N_vals[v]*f_conv_N_vals[B];
				}else{// i is a perfect power of p (and so is m)
					// from the definition, (f <*> N)(p**a) = f(1)*N(p**a) + f(p)*N(p**(a-1)) + ... + f(p**a)*N(1)
					// = f(1)*p**a + f(p)*p**(a-1) + ... + f(p**a)*1 = p*(f <*> N)(p**(a-1)) + f(p**a)
					f_conv_N_vals[m] = p*f_conv_N_vals[i] + f_vals[m];
				}
				break;
			}// otherwise, i and p are coprime, so (f <*> u)(i*p) = (f <*> u)(i)*(f <*> u)(p)
			f_conv_N_vals[m] = f_conv_N_vals[i]*f_conv_N_vals[p];
			smallest_ppow[m] = p;
		}
	}
	free(smallest_ppow);
	free(primes);
	free(is_c_buf);
	return true;
}

bool nut_euler_sieve_conv(int64_t n, const int64_t f_vals[static n+1], const int64_t g_vals[static n+1], int64_t f_conv_vals[static n+1]){
	int64_t *smallest_ppow = malloc((n+1)*sizeof(int64_t));
	uint8_t *is_c_buf = calloc(n/8 + 1, sizeof(uint8_t));
	uint64_t num_primes = 0;
	int64_t *primes = malloc(nut_max_primes_le(n)*sizeof(int64_t));
	if(!smallest_ppow || !is_c_buf || !primes){
		free(smallest_ppow);
		free(is_c_buf);
		free(primes);
		return false;
	}
	f_conv_vals[1] = 1;
	for(int64_t i = 2; i <= n; ++i){
		// if i is prime, add it to the list of primes
		if(!is_composite_unpacked(i, is_c_buf)){
			primes[num_primes++] = i;
			f_conv_vals[i] = f_vals[i] + g_vals[i];
			smallest_ppow[i] = i;
		}
		// the key feature of euler's sieve is that it ONLY marks off multiples
		// of each number i times primes <= i.  This requires storing an array of
		// primes and reading it A LOT of times, but ensures that each composite number
		// is visited only once
		for(uint64_t j = 0; j < num_primes; ++j){
			int64_t p = primes[j], m;
			if(__builtin_mul_overflow(p, i, &m) || m > n){
				break;
			}
			mark_composite_unpacked(m, is_c_buf);
			if(i%p == 0){// i is a multiple of p, so in particular i and p are not coprime
				// also notice that because we break out of this loop at the bottom of this if statement,
				// p is in fact the smallest prime divisor of i
				int64_t ppow = smallest_ppow[i];
				int64_t B = ppow*p;
				smallest_ppow[m] = B;
				int64_t v = i/ppow;
				if(v != 1){// i is not a perfect power of p, that is, i = v*p**a with v != 1
					f_conv_vals[m] = f_conv_vals[v]*f_conv_vals[B];
				}else{// i is a perfect power of p (and so is m)
					// from the definition, (f <*> g)(p**a) = f(1)*g(p**a) + f(p)*g(p**(a-1)) + ... + f(p**a)*g(1)
					int64_t c = f_vals[i] + g_vals[i];
					int64_t A = p;
					while(A < ppow){
						c += f_vals[A]*g_vals[ppow] + f_vals[ppow]*g_vals[A];
						A *= p;
						ppow /= p;
					}
					if(A == ppow){
						c += f_vals[A]*g_vals[A];
					}
					f_conv_vals[m] = c;
				}
				break;
			}// otherwise, i and p are coprime, so (f <*> u)(i*p) = (f <*> u)(i)*(f <*> u)(p)
			f_conv_vals[m] = f_conv_vals[i]*f_conv_vals[p];
			smallest_ppow[m] = p;
		}
	}
	free(smallest_ppow);
	free(primes);
	free(is_c_buf);
	return true;
}

bool nut_Diri_init(nut_Diri *self, int64_t x, int64_t y){
	int64_t ymin = nut_u64_nth_root(x, 2);
	if(y < ymin){
		y = ymin;
	}
	int64_t yinv = x/y;
	if(!(self->buf = malloc((y + yinv + 2)*sizeof(int64_t)))){
		return false;
	}
	self->x = x;
	self->y = y;
	self->yinv = yinv;
	return true;
}

void nut_Diri_destroy(nut_Diri *self){
	free(self->buf);
	*self = (nut_Diri){};
}

void nut_Diri_compute_I(nut_Diri *self){
	memset(self->buf, 0, (self->y + 1)*sizeof(int64_t));
	self->buf[1] = 1;
	for(int64_t i = 1; i <= self->yinv; ++i){
		self->buf[self->y + i] = 1;
	}
}

void nut_Diri_compute_u(nut_Diri *self){
	for(int64_t i = 0; i <= self->y; ++i){
		self->buf[i] = 1;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		self->buf[self->y + i] = self->x/i;
	}
}

void nut_Diri_compute_N(nut_Diri *self){
	for(int64_t i = 0; i <= self->y; ++i){
		self->buf[i] = i;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		self->buf[self->y + i] = v*(v + 1)/2;
	}
}

void nut_Diri_compute_mertens(nut_Diri *self, const uint8_t mobius[self->y/4 + 1]){
	nut_Diri_set_dense(self, 0, 0);
	nut_Diri_set_dense(self, 1, 1);
	for(int64_t i = 2, acc = 1, v; i <= self->y; ++i){
		if((v = nut_Bitfield2_arr_get(mobius, i)) == 3){
			v = -1;
		}
		nut_Diri_set_dense(self, i, acc += v);
	}
	for(int64_t i = self->yinv; i; --i){
		int64_t v = self->x/i;
		int64_t M = 1;
		int64_t vr = nut_u64_nth_root(v, 2);
		for(int64_t j = 1; j <= vr; ++j){
			M -= (nut_Diri_get_dense(self, j) - nut_Diri_get_dense(self, j - 1))*(v/j);
		}
		for(int64_t j = 2; j <= vr; ++j){
			if(i*j > self->yinv){
				M -= nut_Diri_get_dense(self, v/j);
			}else{
				M -= nut_Diri_get_sparse(self, i*j);
			}
		}
		if(vr <= self->y){
			M += nut_Diri_get_dense(self, vr)*vr;
		}else{
			M += nut_Diri_get_sparse(self, self->x/vr)*vr;
		}
		nut_Diri_set_sparse(self, i, M);
	}
	for(int64_t i = 2, v; i <= self->y; ++i){
		if((v = nut_Bitfield2_arr_get(mobius, i)) == 3){
			v = -1;
		}
		nut_Diri_set_dense(self, i, v);
	}
}

bool nut_Diri_compute_conv_u(nut_Diri *self, const nut_Diri *f_tbl){
	if(self->y != f_tbl->y || self->x != f_tbl->x){
		return false;
	}
	// use the dense part of self to temporarily store the sums of f for small n up to y
	self->buf[0] = 0;
	for(int64_t i = 1; i <= self->y; ++i){
		self->buf[i] = self->buf[i-1] + f_tbl->buf[i];
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);// TODO: v is close to the previous v, so only one newton step should be needed here
		int64_t h = 0;
		for(int64_t n = 1; n <= vr; ++n){
			if(v/n <= self->y){
				h += nut_Diri_get_dense(self, v/n);
			}else{
				h += nut_Diri_get_sparse(f_tbl, i*n);
			}
			h += nut_Diri_get_dense(f_tbl, n)*(v/n);
		}
		h -= nut_Diri_get_dense(self, vr)*vr;
		nut_Diri_set_sparse(self, i, h);
	}
	return nut_euler_sieve_conv_u(self->y, f_tbl->buf, self->buf);
}

bool nut_Diri_compute_conv_N(nut_Diri *self, const nut_Diri *f_tbl){
	if(self->y != f_tbl->y || self->x != f_tbl->x){
		return false;
	}
	// use the dense part of self to temporarily store the sums of f for small n up to y
	self->buf[0] = 0;
	for(int64_t i = 1; i <= self->y; ++i){
		self->buf[i] = self->buf[i-1] + f_tbl->buf[i];
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);
		int64_t h = 0;
		// by hyperbola formula:
		// (f <*> N)(v) = sum(n = 1 ... vr, f(n)*(v/n)*((v/n) + 1)/2) + sum(n = 1 ... vr, F(v/n)*n) - F(vr)*vr*(vr + 1)/2
		// both sums are folded together into the following loop
		for(int64_t n = 1; n <= vr; ++n){
			if(v/n <= self->y){
				h += nut_Diri_get_dense(self, v/n)*n;
			}else{
				h += nut_Diri_get_sparse(f_tbl, i*n)*n;
			}
			h += nut_Diri_get_dense(f_tbl, n)*((v/n)*(v/n + 1)/2);
		}
		// finally we apply the corrective term (remove double counted values)
		h -= nut_Diri_get_dense(self, vr)*(vr*(vr + 1)/2);
		nut_Diri_set_sparse(self, i, h);
	}
	return nut_euler_sieve_conv_N(self->y, f_tbl->buf, self->buf);
}

bool nut_Diri_compute_conv(nut_Diri *self, const nut_Diri *f_tbl, const nut_Diri *g_tbl){
	if(self->y != f_tbl->y || self->x != f_tbl->x){
		return false;
	}
	// use the dense part of self to temporarily store the sums of f for small n up to y
	self->buf[0] = 0;
	for(int64_t i = 1; i <= self->y; ++i){
		self->buf[i] = self->buf[i-1] + f_tbl->buf[i];
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);
		int64_t h = 0;
		for(int64_t n = 1; n <= vr; ++n){
			if(v/n <= self->y){
				h += nut_Diri_get_dense(self, v/n)*nut_Diri_get_dense(g_tbl, n);
			}else{
				h += nut_Diri_get_sparse(f_tbl, i*n)*nut_Diri_get_dense(g_tbl, n);
			}
		}
		h -= nut_Diri_get_dense(self, vr)*(vr*(vr + 1)/2);
		nut_Diri_set_sparse(self, i, h);
	}
	// use the dense part of self to temporarily store the sums of g for small n up to y
	for(int64_t i = 1; i <= self->y; ++i){
		self->buf[i] = self->buf[i-1] + g_tbl->buf[i];
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);
		int64_t h = nut_Diri_get_sparse(self, i);
		for(int64_t n = 1; n <= vr; ++n){
			if(v/n <= self->y){
				h += nut_Diri_get_dense(f_tbl, n)*nut_Diri_get_dense(self, v/n);
			}else{
				h += nut_Diri_get_dense(f_tbl, n)*nut_Diri_get_sparse(g_tbl, i*n);
			}
		}
		nut_Diri_set_sparse(self, i, h);
	}
	return nut_euler_sieve_conv(self->y, f_tbl->buf, g_tbl->buf, self->buf);
}

