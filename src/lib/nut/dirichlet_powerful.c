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

typedef struct{
	uint64_t n;
	int64_t hn;
	uint64_t i;
} PowerfulStackEntry;

typedef struct{
	uint64_t len, cap;
	PowerfulStackEntry *entries;
} PowerfulStack;

static bool PowerfulStack_init(PowerfulStack *self, uint64_t max){
	if(!max){
		*self = (PowerfulStack){};
		return true;
	}
	self->cap = 64 - __builtin_clzll(max);// floor(log_2(max))
	if(!(self->entries = malloc(self->cap*sizeof(*self->entries)))){
		return false;
	}
	self->entries[0] = (PowerfulStackEntry){.n=1, .hn=1, .i=0};
	self->len = 1;
	return true;
}

static void PowerfulStack_destroy(PowerfulStack *self){
	free(self->entries);
	*self = (PowerfulStack){};
}

static bool PowerfulStack_push(PowerfulStack *self, const PowerfulStackEntry *ent){
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

static bool PowerfulStack_pop(PowerfulStack *self, PowerfulStackEntry *out){
	if(!self->len){
		return false;
	}
	memcpy(out, self->entries + --self->len, sizeof(*self->entries));
	return true;
}

// F(max) = sum(n powerful <= max : h(n)G(max/n)), since h(n) is 0 at non-powerful values of n
bool nut_Diri_sum_adjusted_he(int64_t *restrict out, int64_t m, const nut_Diri *restrict g_tbl, const int64_t *restrict h_vals){
	// Allocate and fill a temporary array to store the sums of g up to g_tbl->y
	int64_t *g_dense [[gnu::cleanup(cleanup_free)]] = malloc((g_tbl->y + 1)*sizeof(int64_t));
	uint64_t xr = nut_u64_nth_root(g_tbl->x, 2);
	uint64_t num_primes;
	uint64_t *primes [[gnu::cleanup(cleanup_free)]] = nut_sieve_primes(xr, &num_primes);
	PowerfulStack gens [[gnu::cleanup(PowerfulStack_destroy)]] = {};
	if(!g_dense || !primes || !PowerfulStack_init(&gens, g_tbl->x)){
		return false;
	}
	g_dense[0] = 0;
	for(int64_t i = 1; i <= g_tbl->y; ++i){
		int64_t term = g_dense[i - 1] + nut_Diri_get_dense(g_tbl, i);
		g_dense[i] = m ? nut_i64_mod(term, m) : term;
	}
	int64_t res = 0;
	while(gens.len){
		PowerfulStackEntry ent;
		PowerfulStack_pop(&gens, &ent);
		uint64_t p = primes[ent.i], term;
		if(ent.i == num_primes || p*p > g_tbl->x/ent.n){
			if(ent.n > (uint64_t)g_tbl->yinv){
				term = ent.hn*g_dense[g_tbl->x/ent.n];
			}else{
				term = ent.hn*nut_Diri_get_sparse(g_tbl, ent.n);
			}
			if(m){
				res = nut_i64_mod(res + term, m);
			}else{
				res += term;
			}
			continue;
		}
		++ent.i;
		if(!PowerfulStack_push(&gens, &ent)){
			return false;
		}
		PowerfulStackEntry m_ent = ent;
		for(uint64_t pp = p, e = 2; !__builtin_mul_overflow(pp, p, &pp) && pp <= g_tbl->x/ent.n; ++e){
			m_ent.n = ent.n*pp;
			m_ent.hn = ent.hn*h_vals[e];
			if(!PowerfulStack_push(&gens, &m_ent)){
				return false;
			}
		}
	}
	*out = res;
	return true;
}

bool nut_Diri_sum_adjusted_hpe(int64_t *restrict out, int64_t m, const nut_Diri *restrict g_tbl, int64_t (*h_fn)(uint64_t p, uint64_t pp, uint64_t e, int64_t m)){
	int64_t *g_dense [[gnu::cleanup(cleanup_free)]] = malloc((g_tbl->y + 1)*sizeof(int64_t));
	uint64_t xr = nut_u64_nth_root(g_tbl->x, 2);
	uint64_t num_primes;
	uint64_t *primes [[gnu::cleanup(cleanup_free)]] = nut_sieve_primes(xr, &num_primes);
	PowerfulStack gens [[gnu::cleanup(PowerfulStack_destroy)]] = {};
	if(!g_dense || !primes || !PowerfulStack_init(&gens, g_tbl->x)){
		return false;
	}
	g_dense[0] = 0;
	for(int64_t i = 1; i <= g_tbl->y; ++i){
		int64_t term = g_dense[i - 1] + nut_Diri_get_dense(g_tbl, i);
		g_dense[i] = m ? nut_i64_mod(term, m) : term;
	}
	int64_t res = 0;
	while(gens.len){
		PowerfulStackEntry ent;
		PowerfulStack_pop(&gens, &ent);
		uint64_t p = primes[ent.i], term;
		if(ent.i == num_primes || p*p > g_tbl->x/ent.n){
			if(ent.n > (uint64_t)g_tbl->yinv){
				term = ent.hn*g_dense[g_tbl->x/ent.n];
			}else{
				term = ent.hn*nut_Diri_get_sparse(g_tbl, ent.n);
			}
			if(m){
				res = nut_i64_mod(res + term, m);
			}else{
				res += term;
			}
			continue;
		}
		++ent.i;
		if(!PowerfulStack_push(&gens, &ent)){
			return false;
		}
		PowerfulStackEntry m_ent = ent;
		for(uint64_t pp = p, e = 2; !__builtin_mul_overflow(pp, p, &pp) && pp <= g_tbl->x/ent.n; ++e){
			m_ent.n = ent.n*pp;
			m_ent.hn = ent.hn*h_fn(p, pp, e, m);
			if(!PowerfulStack_push(&gens, &m_ent)){
				return false;
			}
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

