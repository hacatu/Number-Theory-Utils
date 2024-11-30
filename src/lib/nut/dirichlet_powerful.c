#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <nut/debug.h>
#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/dirichlet.h>
#include <nut/dirichlet_powerful.h>
#include <nut/sieves.h>

static bool PfIt_init(nut_PfIt *self, uint64_t max, uint64_t modulus, uint64_t small_primes){
	if(!max){
		*self = (nut_PfIt){};
		return true;
	}
	self->small_primes = small_primes;
	self->cap = 63 - __builtin_clzll(max);// floor(log_2(max))
	self->rt_max = nut_u64_nth_root(max, 2);
	if(!(self->entries = malloc(self->cap*sizeof(*self->entries)))){
		return false;
	}else if(!(self->primes = nut_sieve_primes(self->rt_max, &self->num_primes))){
		free(self->entries);
		self->entries = NULL;
		return false;
	}
	self->entries[0] = (nut_PfStackEnt){.n=1, .hn=1, .i=0};
	self->len = 1;
	self->max = max;
	self->modulus = modulus;
	return true;
}

bool nut_PfIt_init_fn(nut_PfIt *self, uint64_t max, uint64_t modulus, uint64_t small_primes, int64_t (*h_fn)(uint64_t p, uint64_t pp, uint64_t e, uint64_t m)){
	if(!PfIt_init(self, max, modulus, small_primes)){
		return false;
	}
	self->h_kind = NUT_DIRI_H_FN;
	self->h_fn = h_fn;
	return true;
}

bool nut_PfIt_init_hvals(nut_PfIt *restrict self, uint64_t max, uint64_t modulus, uint64_t small_primes, const int64_t *restrict h_vals){
	if(!PfIt_init(self, max, modulus, small_primes)){
		return false;
	}
	self->h_kind = NUT_DIRI_H_VALS;
	self->h_vals = h_vals;
	return true;
}

bool nut_PfIt_init_hseqs(nut_PfIt *restrict self, uint64_t max, uint64_t modulus, uint64_t small_primes,
	int64_t (*f_fn)(uint64_t p, uint64_t pp, uint64_t e, uint64_t m), int64_t (*g_fn)(uint64_t p, uint64_t pp, uint64_t e, uint64_t m)
){
	if(!PfIt_init(self, max, modulus, small_primes)){
		return false;
	}
	uint64_t values_cap = 3*self->num_primes;
	self->h_kind = NUT_DIRI_H_SEQS;
	self->h_seqs.offsets = malloc(self->num_primes*sizeof(uint64_t));
	self->h_seqs.values = malloc(values_cap*sizeof(int64_t));
	if(!self->h_seqs.offsets || !self->h_seqs.values){
		nut_PfIt_destroy(self);
		return false;
	}
	uint64_t base_offset = 0;
	uint64_t curr_max_pow = 63 - __builtin_clzll(max);
	int64_t *f_series [[gnu::cleanup(cleanup_free)]] = malloc((curr_max_pow + 1)*sizeof(int64_t));
	int64_t *g_series [[gnu::cleanup(cleanup_free)]] = malloc((curr_max_pow + 1)*sizeof(int64_t));
	int64_t *h_series [[gnu::cleanup(cleanup_free)]] = malloc((curr_max_pow + 1)*sizeof(int64_t));
	if(!f_series || !g_series || !h_series){
		nut_PfIt_destroy(self);
		return false;
	}
	f_series[0] = g_series[0] = 1;
	uint64_t curr_max_prime = 2;
	for(uint64_t i = 0; i < self->num_primes; ++i){
		uint64_t p = self->primes[i];
		// ensure that p^curr_max_pow is still <= max
		while(curr_max_pow > 2 && p > curr_max_prime){
			// this must be a while loop because for very small primes like 2 it's possible for the max power to drop more than 1
			// also when the max power is 2 we don't need to check anymore since we've only gathered the primes up to sqrt(max)
			--curr_max_pow;
			curr_max_prime = nut_u64_nth_root(max, curr_max_pow);
		}
		// ensure that the values table is large enough
		uint64_t min_pow = i >= self->small_primes ? 2 : 1;
		if(base_offset + curr_max_pow - min_pow >= values_cap){
			uint64_t new_cap = base_offset + 2*(self->num_primes - i);
			int64_t *tmp = realloc(self->h_seqs.values, new_cap*sizeof(int64_t));
			if(!tmp){
				nut_PfIt_destroy(self);
				return false;
			}
			self->h_seqs.values = tmp;
			values_cap = new_cap;
		}
		for(uint64_t pp = 1, e = 1; e <= curr_max_pow; ++e){
			pp *= p;
			f_series[e] = f_fn(p, pp, e, modulus);
			g_series[e] = g_fn(p, pp, e, modulus);
		}
		nut_series_div(curr_max_pow + 1, modulus, h_series, f_series, g_series);
		self->h_seqs.offsets[i] = base_offset;
		memcpy(self->h_seqs.values + base_offset, h_series + min_pow, (curr_max_pow + 1 - min_pow)*sizeof(int64_t));
		base_offset += curr_max_pow + 1 - min_pow;
	}
	return true;
}

void nut_PfIt_destroy(nut_PfIt *self){
	free(self->entries);
	free(self->primes);
	if(self->h_kind == NUT_DIRI_H_SEQS){
		free(self->h_seqs.offsets);
		free(self->h_seqs.values);
	}
	*self = (nut_PfIt){};
}

bool nut_PfStack_push(nut_PfIt *restrict self, const nut_PfStackEnt *restrict ent){
	if(self->len == self->cap){
		uint64_t new_cap = (self->cap << 1) ?: 8;// it shouldn't be possible to get here with cap = 0, but work around it anyway
		void *tmp = realloc(self->entries, new_cap*sizeof(*self->entries));
		if(!tmp){
			return false;
		}
		self->cap = new_cap;
		self->entries = tmp;
	}
	memcpy(self->entries + self->len++, ent, sizeof(*self->entries));
	return true;
}

bool nut_PfStack_pop(nut_PfIt *restrict self, nut_PfStackEnt *restrict out){
	if(!self->len){
		return false;
	}
	memcpy(out, self->entries + --self->len, sizeof(*self->entries));
	return true;
}

bool nut_PfIt_next(nut_PfIt *restrict self, nut_PfStackEnt *restrict out){
	// Yields (n, h(n) mod m) where n are the O(sqrt x) powerful numbers
	// up to x, and h is any multiplicative function.
	while(self->len){
		nut_PfStack_pop(self, out);
		if(out->i >= self->num_primes){
			return true;
		}
		uint64_t p = self->primes[out->i];
		uint64_t min_pow = out->i >= self->small_primes ? 2 : 1;
		if(min_pow == 2){
			if(p*p > self->max/out->n){
				return true;
			}
		}else{
			if(p > self->max/out->n){
				return true;
			}
		}
		if(!nut_PfStack_push(self, &(nut_PfStackEnt){.n=out->n, .hn=out->hn, .i= out->i + 1})){
			return false;
		}
		for(uint64_t pp = min_pow == 2 ? p : 1, e = min_pow; !__builtin_mul_overflow(pp, p, &pp) && pp <= self->max/out->n; ++e){
			int64_t v = out->hn;
			if(self->h_kind == NUT_DIRI_H_VALS){
				v *= self->h_vals[e];
			}else if(self->h_kind == NUT_DIRI_H_FN){
				v *= self->h_fn(p, pp, e, self->modulus);
			}else if(self->h_kind == NUT_DIRI_H_SEQS){
				v *= self->h_seqs.values[self->h_seqs.offsets[out->i] + e - min_pow];
			}else{
				return false;
			}
			if(self->modulus){
				v = nut_i64_mod(v, self->modulus);
			}
			if(!nut_PfStack_push(self, &(nut_PfStackEnt){.n = out->n*pp, .hn = v, .i = out->i + 1})){
				return false;
			}
		}
	}
	return false;
}

bool nut_Diri_sum_adjusted(int64_t *restrict out, const nut_Diri *restrict g_tbl, nut_PfIt *pf_it){
	int64_t *g_dense [[gnu::cleanup(cleanup_free)]] = malloc((g_tbl->y + 1)*sizeof(int64_t));
	if(!g_dense || pf_it->modulus > INT64_MAX){
		return false;
	}
	g_dense[0] = 0;
	int64_t m = pf_it->modulus;
	for(int64_t i = 1; i <= g_tbl->y; ++i){
		int64_t term = g_dense[i - 1] + nut_Diri_get_dense(g_tbl, i);
		g_dense[i] = m ? nut_i64_mod(term, m) : term;
	}
	int64_t res = 0;
	nut_PfStackEnt ent;
	while(nut_PfIt_next(pf_it, &ent)){
		int64_t Gn;
		if(ent.n >= (uint64_t)g_tbl->yinv){
			Gn = g_dense[g_tbl->x/ent.n];
		}else{
			Gn = nut_Diri_get_sparse(g_tbl, ent.n);
		}
		res += ent.hn*Gn;
		if(m){
			res = nut_i64_mod(res, m);
		}
	}
	*out = res;
	return true;
}

bool nut_Diri_sum_u_adjusted(int64_t *restrict out, nut_PfIt *pf_it){
	int64_t m = pf_it->modulus;
	int64_t res = 0;
	uint64_t max = pf_it->max;
	nut_PfStackEnt ent;
	while(nut_PfIt_next(pf_it, &ent)){
		res += ent.hn*(max/ent.n);
		if(m){
			res = nut_i64_mod(res, m);
		}
	}
	*out = res;
	return true;
}

void nut_series_div(uint64_t n, int64_t m, int64_t h[restrict static n], int64_t f[restrict static n], int64_t g[restrict static n]){
	for(uint64_t e = 0; e < n; ++e){
		int64_t term = f[e];
		for(uint64_t k = 1; k <= e; ++k){
			term -= g[k]*h[e - k];
			if(m){
				term %= m;
			}
		}
		h[e] = (m && term < 0) ? m + term : term;
	}
}

