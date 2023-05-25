#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/dirichlet.h>
#include <nut/sieves.h>

uint64_t nut_dirichlet_D(uint64_t max, uint64_t m){
	uint64_t y = nut_u64_nth_root(max, 2);
	uint64_t res = 0;
	for(uint64_t n = 1; n <= y; ++n){
		uint64_t term = max/n;
		res = m ? (res + term)%m : res + term;
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
// a basic sieve of eratosthenes.  However, preliminary testing shows that the sieve of Eratosthenes
// is only about 14% faster than Euler's sieve in the release build, so it isn't worth upgrading.
// Euler's sieve generally works as follows:
// 1. Allocate result[1, ..., n] (in this case, this is f_conv_u_vals and smallest_ppow)
// 2. Initialize is_c_buf[1, ..., n] = false (array of bit flags for each number from 1 to n, 0 for primes and 1 for composites)
// 3. Initialize primes = [] a list of primes
// 4. For i = 2, ..., n:
// 5.   If i is not marked as composite by now (checked using is_c_buf), add it to the list of primes primes
// 6.   Either way, for every prime p in primes:
// 7.     If p divides i: (p is the smallest prime divisor of i since in this case we will break out of this loop)
// 8.       If the largest power of the smallest prime divisor of i (smallest_ppow[i]) equals i, then i is a perfect power of p, and so is p*i
// 9.         In this case, we must compute (f <*> g)(p*i) using the definition of dirichlet convolution
// 10.        We can set smallest_ppow[i*p] = i*p
// 11.      Otherwise, p divides i but i is not a perfect power of p, so we can write i*p as smallest_ppow[i]*p * i/smallest_ppow[i]
// 12.        These two parts are coprime, so we can find (f <*> g)(i*p) = (f <*> g)(smallest_ppow[i]*p) * (f <*> g)(i/smallest_ppow[i])
// 13.        And smallest_ppow[i*p] = smallest_ppow[i]*p
// 14.      Either way, break out of the loop over p
// 15.    Otherwise, p does not divide i, so p and i are already coprime (and p is smaller than any prime divisor of i)
// 16.      So we can set (f <*> g)(i*p) = (f <*> g)(i) * (f <*> g)(p)
// 17.      And smallest_ppow[i*p] = p
// Notice that we don't really need the smallest_ppow array, since we only use it when we know p the smallest prime divisor of i,
// we could just use a while loop to count how many times p divides i or find the largest power of p which divides i, but this is a space/time tradeoff.
// Also notice that we can cache smallest_ppow, primes, and is_c_buf, although the current implementation does not do that.
// On the other hand, the sieve of Eratosthenes generally works as follows:
// 1. Initialize result[1, ..., n] = 1
// 2. For i = 2, ..., n:
// 3.   If result[i] != 1, i is not prime so we can continue
// 4.   Otherwise, i is prime (with a caveat noted below)
// 5.   For every power pp of i which is <= n:
// 6.     Compute (f <*> g)(pp) using the definition of dirichlet convolution, and store this as h
// 7.     For every multiple m of pp <= n starting at 2*pp WHICH IS NOT a multiple of pp*p:
// 8.       Multiply (f <*> g)(m) by h
// When i is prime, result[i] will remain set to 1 when we get to the ith loop.
// HOWEVER, when (f <*> g)(n) can be 1 for composite values of n, this is a problem because result[i] = 1 no longer means i is prime.
// For convolutions which can be 1 for composite inputs, we would need to allocate a separate array to store whether or not each i is prime.
bool nut_euler_sieve_conv_u(int64_t n, int64_t modulus, const int64_t f_vals[restrict static n+1], int64_t f_conv_u_vals[restrict static n+1]){
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
			int64_t term = f_vals[i] + 1;
			f_conv_u_vals[i] = modulus ? nut_i64_mod(term, modulus) : term;
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
					int64_t term = f_conv_u_vals[v]*f_conv_u_vals[B];
					f_conv_u_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
				}else{// i is a perfect power of p (and so is m)
					// from the definition, (f <*> u)(p**a) = f(1)*u(p**a) + f(p)*u(p**(a-1)) + ... + f(p**a)*u(1)
					// = f(1) + f(p) + ... + f(p**a) = (f <*> u)(p**(a-1)) + f(p**a)
					int64_t term = f_conv_u_vals[i] + f_vals[m];
					f_conv_u_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
				}
				break;
			}// otherwise, i and p are coprime, so (f <*> u)(i*p) = (f <*> u)(i)*(f <*> u)(p)
			int64_t term = f_conv_u_vals[i]*f_conv_u_vals[p];
			f_conv_u_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
			smallest_ppow[m] = p;
		}
	}
	free(smallest_ppow);
	free(primes);
	free(is_c_buf);
	return true;
}

bool nut_euler_sieve_conv_N(int64_t n, int64_t modulus, const int64_t f_vals[restrict static n+1], int64_t f_conv_N_vals[restrict static n+1]){
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
			int64_t term = f_vals[i] + i;
			f_conv_N_vals[i] = modulus ? nut_i64_mod(term, modulus) : term;
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
					int64_t term = f_conv_N_vals[v]*f_conv_N_vals[B];
					f_conv_N_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
				}else{// i is a perfect power of p (and so is m)
					// from the definition, (f <*> N)(p**a) = f(1)*N(p**a) + f(p)*N(p**(a-1)) + ... + f(p**a)*N(1)
					// = f(1)*p**a + f(p)*p**(a-1) + ... + f(p**a)*1 = p*(f <*> N)(p**(a-1)) + f(p**a)
					if(modulus){
						int64_t term = (p%modulus)*f_conv_N_vals[i] + f_vals[m];
						f_conv_N_vals[m] = nut_i64_mod(term, modulus);
					}else{
						f_conv_N_vals[m] = p*f_conv_N_vals[i] + f_vals[m];
					}
				}
				break;
			}// otherwise, i and p are coprime, so (f <*> u)(i*p) = (f <*> u)(i)*(f <*> u)(p)
			int64_t term = f_conv_N_vals[i]*f_conv_N_vals[p];
			f_conv_N_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
			smallest_ppow[m] = p;
		}
	}
	free(smallest_ppow);
	free(primes);
	free(is_c_buf);
	return true;
}

bool nut_euler_sieve_conv(int64_t n, int64_t modulus, const int64_t f_vals[static n+1], const int64_t g_vals[static n+1], int64_t f_conv_vals[restrict static n+1]){
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
			int64_t term = f_vals[i] + g_vals[i];
			f_conv_vals[i] = modulus ? nut_i64_mod(term, modulus) : term;
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
					int64_t term = f_conv_vals[v]*f_conv_vals[B];
					f_conv_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
				}else{// i is a perfect power of p (and so is m)
					// from the definition, (f <*> g)(p**a) = f(1)*g(p**a) + f(p)*g(p**(a-1)) + ... + f(p**a)*g(1)
					int64_t c = f_vals[m] + g_vals[m];
					if(modulus){
						c = nut_i64_mod(c, modulus);
					}
					int64_t A = p;
					while(A < ppow){
						if(modulus){
							c = nut_i64_mod(c + f_vals[A]*g_vals[ppow] + f_vals[ppow]*g_vals[A], modulus);
						}else{
							c += f_vals[A]*g_vals[ppow] + f_vals[ppow]*g_vals[A];
						}
						A *= p;
						ppow /= p;
					}
					if(A == ppow){
						if(modulus){
							c = nut_i64_mod(c + f_vals[A]*g_vals[A], modulus);
						}else{
							c += f_vals[A]*g_vals[A];
						}
					}
					f_conv_vals[m] = c;
				}
				break;
			}// otherwise, i and p are coprime, so (f <*> u)(i*p) = (f <*> u)(i)*(f <*> u)(p)
			int64_t term = f_conv_vals[i]*f_conv_vals[p];
			f_conv_vals[m] = modulus ? nut_i64_mod(term, modulus) : term;
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

void nut_Diri_copy(nut_Diri *restrict dest, const nut_Diri *restrict src){
	memcpy(dest->buf, src->buf, (src->y + src->yinv + 1)*sizeof(int64_t));
}

void nut_Diri_compute_I(nut_Diri *self){
	memset(self->buf, 0, (self->y + 1)*sizeof(int64_t));
	self->buf[1] = 1;
	for(int64_t i = 1; i <= self->yinv; ++i){
		self->buf[self->y + i] = 1;
	}
}

void nut_Diri_compute_u(nut_Diri *self, int64_t m){
	for(int64_t i = 0; i <= self->y; ++i){
		self->buf[i] = 1;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t term = self->x/i;
		self->buf[self->y + i] = m ? nut_i64_mod(term, m) : term;
	}
}

void nut_Diri_compute_N(nut_Diri *self, int64_t m){
	for(int64_t i = 0; i <= self->y; ++i){
		self->buf[i] = i;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t term = (v&1) ? v*((v + 1) >> 1) : (v + 1)*(v >> 1);
		if(m){
			self->buf[self->y + i] = term%m;
		}else{
			self->buf[self->y + i] = term;
		}
	}
}

void nut_Diri_compute_mertens(nut_Diri *restrict self, int64_t m, const uint8_t mobius[restrict static self->y/4 + 1]){
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
			int64_t term = (nut_Diri_get_dense(self, j) - nut_Diri_get_dense(self, j - 1))*(v/j);
			if(m){
				M = nut_i64_mod(M - term, m);
			}else{
				M -= term;
			}
		}
		for(int64_t j = 2, term; j <= vr; ++j){
			if(i*j > self->yinv){
				term = nut_Diri_get_dense(self, v/j);
			}else{
				term = nut_Diri_get_sparse(self, i*j);
			}
			if(m){
				M = nut_i64_mod(M - term, m);
			}else{
				M -= term;
			}
		}
		int64_t term;
		if(vr <= self->y){
			term = nut_Diri_get_dense(self, vr)*vr;
		}else{
			term = nut_Diri_get_sparse(self, self->x/vr)*vr;
		}
		if(m){
			M = nut_i64_mod(M + term, m);
		}else{
			M += term;
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

bool nut_Diri_compute_conv_u(nut_Diri *restrict self, int64_t m, const nut_Diri *restrict f_tbl){
	if(self->y != f_tbl->y || self->x != f_tbl->x){
		return false;
	}
	// use the dense part of self to temporarily store the sums of f for small n up to y
	self->buf[0] = 0;
	for(int64_t i = 1; i <= self->y; ++i){
		int64_t term = self->buf[i-1] + f_tbl->buf[i];
		self->buf[i] = m ? nut_i64_mod(term, m) : term;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);// TODO: v is close to the previous v, so only one newton step should be needed here
		int64_t h = 0;
		for(int64_t n = 1, term; n <= vr; ++n){
			if(v/n <= self->y){
				term = nut_Diri_get_dense(self, v/n);
			}else{
				term = nut_Diri_get_sparse(f_tbl, i*n);
			}
			if(m){
				h = nut_i64_mod(h + term, m);
			}else{
				h += term;
			}
			term = nut_Diri_get_dense(f_tbl, n)*(v/n);
			if(m){
				h = nut_i64_mod(h + term, m);
			}else{
				h += term;
			}
		}
		int64_t term = nut_Diri_get_dense(self, vr)*vr;
		if(m){
			h = nut_i64_mod(h - term, m);
		}else{
			h -= term;
		}
		nut_Diri_set_sparse(self, i, h);
	}
	return nut_euler_sieve_conv_u(self->y, m, f_tbl->buf, self->buf);
}

bool nut_Diri_compute_conv_N(nut_Diri *restrict self, int64_t m, const nut_Diri *restrict f_tbl){
	if(self->y != f_tbl->y || self->x != f_tbl->x){
		return false;
	}
	// use the dense part of self to temporarily store the sums of f for small n up to y
	self->buf[0] = 0;
	for(int64_t i = 1; i <= self->y; ++i){
		int64_t term = self->buf[i-1] + f_tbl->buf[i];
		self->buf[i] = m ? nut_i64_mod(term, m) : term;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);
		int64_t h = 0;
		// by hyperbola formula:
		// (f <*> N)(v) = sum(n = 1 ... vr, f(n)*(v/n)*((v/n) + 1)/2) + sum(n = 1 ... vr, F(v/n)*n) - F(vr)*vr*(vr + 1)/2
		// both sums are folded together into the following loop
		for(int64_t n = 1, term; n <= vr; ++n){
			if(v/n <= self->y){
				term = nut_Diri_get_dense(self, v/n)*n;
			}else{
				term = nut_Diri_get_sparse(f_tbl, i*n)*n;
			}
			if(m){
				h = nut_i64_mod(h + term, m);
			}else{
				h += term;
			}
			term = nut_Diri_get_dense(f_tbl, n);
			int64_t Gvn = ((v/n)&1) ? (v/n)*((v/n + 1) >> 1) : ((v/n) >> 1)*(v/n + 1);
			if(m){
				Gvn %= m;
				h = nut_i64_mod(h + term*Gvn, m);
			}else{
				h += term*Gvn;
			}
		}
		// finally we apply the corrective term (remove double counted values)
		int64_t term = nut_Diri_get_dense(self, vr);
		int64_t Gvr = (vr&1) ? vr*((vr + 1) >> 1) : (vr >> 1)*(vr + 1);
		if(m){
			Gvr %= m;
			h = nut_i64_mod(h - term*Gvr, m);
		}else{
			h -= term*Gvr;
		}
		nut_Diri_set_sparse(self, i, h);
	}
	return nut_euler_sieve_conv_N(self->y, m, f_tbl->buf, self->buf);
}

bool nut_Diri_compute_conv(nut_Diri *restrict self, int64_t m, const nut_Diri *f_tbl, const nut_Diri *g_tbl){
	if(self->y != f_tbl->y || self->x != f_tbl->x){
		return false;
	}
	// use the dense part of self to temporarily store the sums of f for small n up to y
	self->buf[0] = 0;
	for(int64_t i = 1; i <= self->y; ++i){
		int64_t term = self->buf[i-1] + f_tbl->buf[i];
		self->buf[i] = m ? nut_i64_mod(term, m) : term;
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);
		int64_t h = 0;
		for(int64_t n = 1, term; n <= vr; ++n){
			if(v/n <= self->y){
				term = nut_Diri_get_dense(self, v/n)*nut_Diri_get_dense(g_tbl, n);
			}else{
				term = nut_Diri_get_sparse(f_tbl, i*n)*nut_Diri_get_dense(g_tbl, n);
			}
			if(m){
				h = nut_i64_mod(h + term, m);
			}else{
				h += term;
			}
		}
		nut_Diri_set_sparse(self, i, h);
	}
	// use the dense part of self to temporarily store the sums of g for small n up to y
	// simultaniously, apply the adjustments to the h values.
	// H(x/j) has an adjustment of F(i)V(i) where i = sqrt(x/j)
	// so for every i = sqrt(x/j) value, we need to adjust all H(x/j) values
	// for j in the range (x/(i + 1)**2, x/i**2] && [1, yinv]
	// So we can ignore some small i values where this intersection is empty, that is, yinv <= x/(i + 1)**2
	// (i + 1)**2 <= x/yinv <=> i + 1 <= sqrt(x/yinv)
	// that is, for i <= sqrt(x/yinv) - 1, there are no corresponding j to adjust
	// note that x/yinv = y and i <= sqrt(y) - 1 <=> i < sqrt(y)
	int64_t i_ub = nut_u64_nth_root(self->y, 2);
	for(int64_t i = 1; i < i_ub; ++i){
		int64_t term = self->buf[i-1] + g_tbl->buf[i];
		self->buf[i] = m ? nut_i64_mod(term, m) : term;
	}
	for(int64_t i = i_ub; i <= self->y; ++i){
		int64_t F = self->buf[i];
		int64_t G = self->buf[i-1] + g_tbl->buf[i];
		if(m){
			G = nut_i64_mod(G, m);
		}
		self->buf[i] = G;
		int64_t j_ub = self->x/(i*i);
		if(j_ub > self->yinv){
			j_ub = self->yinv;
		}
		for(int64_t j = self->x/((i + 1)*(i + 1)) + 1; j <= j_ub; ++j){
			int64_t h = nut_Diri_get_sparse(self, j);
			if(m){
				h = nut_i64_mod(h - F*G, m);
			}else{
				h -= F*G;
			}
			nut_Diri_set_sparse(self, j, h);
		}
	}
	for(int64_t i = 1; i <= self->yinv; ++i){
		int64_t v = self->x/i;
		int64_t vr = nut_u64_nth_root(v, 2);
		int64_t h = nut_Diri_get_sparse(self, i);
		for(int64_t n = 1, term; n <= vr; ++n){
			if(v/n <= self->y){
				term = nut_Diri_get_dense(f_tbl, n)*nut_Diri_get_dense(self, v/n);
			}else{
				term = nut_Diri_get_dense(f_tbl, n)*nut_Diri_get_sparse(g_tbl, i*n);
			}
			if(m){
				h = nut_i64_mod(h + term, m);
			}else{
				h += term;
			}
		}
		nut_Diri_set_sparse(self, i, h);
	}
	return nut_euler_sieve_conv(self->y, m, f_tbl->buf, g_tbl->buf, self->buf);
}

