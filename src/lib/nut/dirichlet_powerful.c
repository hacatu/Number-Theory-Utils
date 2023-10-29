#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include <nut/debug.h>
#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/dirichlet.h>
#include <nut/dirichlet_powerful.h>
#include <nut/sieves.h>

static bool PfIt_init(nut_PfIt *self, uint64_t max, uint64_t modulus){
	if(!max){
		*self = (nut_PfIt){};
		return true;
	}
	self->cap = 64 - __builtin_clzll(max);// floor(log_2(max))
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

bool nut_PfIt_init_fn(nut_PfIt *self, uint64_t max, uint64_t modulus, int64_t (*h_fn)(uint64_t p, uint64_t pp, uint64_t e, uint64_t m)){
	if(!PfIt_init(self, max, modulus)){
		return false;
	}
	self->use_table = false;
	self->h_fn = h_fn;
	return true;
}

bool nut_PfIt_init_hvals(nut_PfIt *restrict self, uint64_t max, uint64_t modulus, const int64_t *restrict h_vals){
	if(!PfIt_init(self, max, modulus)){
		return false;
	}
	self->use_table = true;
	self->h_vals = h_vals;
	return true;
}

void nut_PfIt_destroy(nut_PfIt *self){
	free(self->entries);
	free(self->primes);
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
		if(p*p > self->max/out->n){
			return true;
		}
		if(!nut_PfStack_push(self, &(nut_PfStackEnt){.n=out->n, .hn=out->hn, .i= out->i + 1})){
			return false;
		}
		for(uint64_t pp = p, e = 2; !__builtin_mul_overflow(pp, p, &pp) && pp <= self->max/out->n; ++e){
			int64_t v = out->hn;
			if(self->use_table){
				v *= self->h_vals[e];
			}else{
				v *= self->h_fn(p, pp, e, self->modulus);
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
		if(ent.n > (uint64_t)g_tbl->yinv){
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

