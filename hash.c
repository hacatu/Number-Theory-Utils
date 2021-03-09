#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "hash.h"

//these are the prime infimums of powers of two in the uint64_t range so which one to use can be calculated by a bit scan (__builtin_clzll() aka bsfq)
uint64_t exp_primes[64] = {0UL, 3UL, 7UL, 13UL, 31UL, 61UL, 127UL, 251UL, 509UL, 1021UL, 2039UL, 4093UL, 8191UL, 16381UL, 32749UL, 65537UL, 131071UL, 262139UL, 524287UL, 1048573UL,
	2097143UL, 4194301UL, 8388593UL, 16777213UL, 33554393UL, 67108859UL, 134217689UL, 268435399UL, 536870909UL, 1073741789UL, 2147483647UL, 4294967291UL, 8589934583UL, 17179869143UL,
	34359738337UL, 68719476731UL, 137438953447UL, 274877906899UL, 549755813881UL, 1099511627689UL, 2199023255531UL, 4398046511093UL, 8796093022151UL, 17592186044399UL,
	35184372088777UL, 70368744177643UL, 140737488355213UL, 281474976710597UL, 562949953421231UL, 1125899906842597UL, 2251799813685119UL, 4503599627370449UL, 9007199254740881UL,
	18014398509481951UL, 36028797018963913UL, 72057594037927931UL, 144115188075855859UL, 288230376151711717UL, 576460752303423433UL, 1152921504606846883UL, 2305843009213693951UL,
	4611686018427387847UL, 9223372036854775783UL, 18446744073709551557UL};

int hash_init(hashtbl_t *self, const hashtbl_ft *ft, size_t reserve){
	*self = (hashtbl_t){};
	if(!reserve){
		return 1;
	}
	size_t i = 64 - __builtin_clzll(reserve - 1);
	self->len_a = exp_primes[i];
	self->table_a = calloc(ft->size, self->len_a);
	if(!self->table_a){
		return 0;
	}
	self->flags_a = calloc(sizeof(uint64_t), ((self->len_a - 1) >> 5) + 1);
	if(!self->flags_a){
		free(self->table_a);
		self->table_a = NULL;
		return 0;
	}
	self->cap = self->len_a*(long double)ft->load_factor;
	return 1;
}

inline static uint64_t hash_find_slot_re(hashtbl_t *self, const hashtbl_ft *ft, uint64_t a){
	for(uint64_t j = 0, i; j < self->len_a; ++j){
		i = (a + (j&1 ? j*j : self->len_a - j*j%self->len_a))%self->len_a;
		uint64_t f = (self->flags_a[i >> 5] >> (i&0x1F))&0x100000001ULL;
		if(!(f&0x100000000ULL)){
			return i;
		}
	}
	return ~0ULL;
}

inline static uint64_t hash_find_slot_ap(hashtbl_t *self, const hashtbl_ft *ft, const void *key){
	uint64_t a = ft->hash(key) % self->len_a;
	for(uint64_t j = 0, i; j < self->len_a; ++j){
		i = (a + (j&1 ? j*j : self->len_a - j*j%self->len_a))%self->len_a;
		uint64_t f = (self->flags_a[i >> 5] >> (i&0x1F))&0x100000001ULL;
		if(!(f&0x100000000ULL && ft->cmp(key, self->table_a + i*ft->size))){
			return i;
		}
	}
	return ~0ULL;
}

inline static uint64_t hash_find_slot_ex(hashtbl_t *self, const hashtbl_ft *ft, const void *key){
	uint64_t a = ft->hash(key) % self->len_a;
	for(uint64_t j = 0, i; j < self->len_a; ++j){
		i = (a + (j&1 ? j*j : self->len_a - j*j%self->len_a))%self->len_a;
		uint64_t f = (self->flags_a[i >> 5] >> (i&0x1F))&0x100000001ULL;
		if(!(f&0x100000000ULL)){
			return i;
		}else if(!ft->cmp(key, self->table_a + i*ft->size)){
			return ~0ULL;
		}
	}
	return ~0ULL;
}

inline static int hash_ix_start(hashtbl_t *self, const hashtbl_ft *ft){
	size_t i = self->len_a ? 64 - __builtin_clzll(self->len_a - 1) : 1;
	uint64_t new_len = exp_primes[i];
	void *new_table = calloc(ft->size, new_len);
	if(!new_table){
		return 0;
	}
	uint64_t *new_flags = calloc(sizeof(uint64_t), ((new_len - 1) >> 5) + 1);
	if(!new_flags){
		free(new_table);
		return 0;
	}
	size_t new_cap = new_len*(long double)ft->load_factor;
	//we currently have full entries and can fit new_cap before resizing.  This means we must move full/(new_cap - full)
	//entries on average per insertion.  Because full == cap + 1, new_cap == new_len*load_factor, cap == len_a*load_factor,
	//and new_len == len_a*growth_factor, full/(new_cap - full) ~~ cap/(new_cap - cap) ~~ len_a/(new_len - len_a)
	//~~ len_a/(len_a*growth_factor - len_a) ~~ 1/(growth_factor - 1).  I picked the largest prime number before each power
	//of 2 so full/(new_cap - full) is very close to 1.  It is always less than 2 so I might just use 2 instead of r
	self->table_b = self->table_a;
	self->table_a = new_table;
	self->flags_b = self->flags_a;
	self->flags_a = new_flags;
	self->len_b = self->len_a;
	self->len_a = new_len;
	self->cap = new_cap;
	self->r = ceill(self->full/((long double)new_cap - self->full));
	self->i = 0;
	return 1;
}

inline static void hash_ix_move(hashtbl_t *self, const hashtbl_ft *ft, uint64_t n){
	uint64_t b = self->i;
	for(uint64_t a; b < self->len_b; ++b){
		uint64_t f = (self->flags_b[b >> 5] >> (b&0x1F))&0x100000001ULL;
		if(!(f&0x100000000ULL)){
			continue;
		}
		void *ent = self->table_b + b*ft->size;
		a = ft->hash(ent) % self->len_a;
		a = hash_find_slot_re(self, ft, a);
		memcpy(self->table_a + a*ft->size, ent, ft->size);
		self->flags_b[b >> 5] |= 0x1ULL << (b&0x1F);
		self->flags_b[b >> 5] &= ~(0x100000000ULL << (b&0x1F));
		self->flags_a[a >> 5] |= 0x100000000ULL << (a&0x1F);
		self->flags_a[a >> 5] &= ~(0x1ULL << (a&0x1F));
		if(!--n){
			break;
		}
	}
	self->i = b;
	if(b == self->len_b){
		free(self->table_b);
		free(self->flags_b);
		self->table_b = NULL;
		self->flags_b = NULL;
		self->len_b = 0;
		self->i = 0;
	}
}

inline static uint64_t hash_get_index(void *table, uint64_t *flags, size_t len, uint64_t a, const hashtbl_ft *ft, const void *key){
	a %= len;
	for(uint64_t j = 0, i; j < len; ++j){
		i = (a + (j&1 ? j*j : len - j*j%len))%len;
		uint64_t f = (flags[i >> 5] >> (i&0x1F))&0x100000001ULL;
		if(!f){
			return ~0ULL;
		}
		if(f&0x100000000ULL && !ft->cmp(key, table + i*ft->size)){
			return i;
		}
	}
	return ~0ULL;
}

inline static void *hash_get_single_a(hashtbl_t *self, uint64_t i, const hashtbl_ft *ft, const void *key){
	uint64_t a = hash_get_index(self->table_a, self->flags_a, self->len_a, i, ft, key);
	return ~a ? self->table_a + a*ft->size : NULL;
}

inline static void *hash_get_single_b(hashtbl_t *self, uint64_t i, const hashtbl_ft *ft, const void *key){
	uint64_t b = hash_get_index(self->table_b, self->flags_b, self->len_b, i, ft, key);
	return ~b ? self->table_b + b*ft->size : NULL;
}

inline static void *hash_get_split(hashtbl_t *self, const hashtbl_ft *ft, const void *key){
	uint64_t i = ft->hash(key);
	if(self->i << 1 < self->len_b){
		return hash_get_single_a(self, i, ft, key) ?: hash_get_single_b(self, i, ft, key);
	}
	return hash_get_single_b(self, i, ft, key) ?: hash_get_single_a(self, i, ft, key);
}

void *hash_get(hashtbl_t *self, const hashtbl_ft *ft, const void *key){
	if(!self->table_a){
		return NULL;
	}
	return self->table_b ? hash_get_split(self, ft, key) : hash_get_single_a(self, ft->hash(key), ft, key);
}

void *hash_insert(hashtbl_t *self, const hashtbl_ft *ft, const void *key, int *_status){
	void *ret = NULL;
	int status = 0;
	if(self->full == self->cap){
		if(!hash_ix_start(self, ft)){
			if(_status){
				*_status = 0;
			}
			return NULL;
		}
	}
	if(self->table_b){
		hash_ix_move(self, ft, self->r);
		if(self->table_b){
			if(hash_get_single_b(self, ft->hash(key), ft, key)){
				if(_status){
					*_status = 2;
				}
				return NULL;
			}
		}
	}
	uint64_t i = hash_find_slot_ap(self, ft, key);
	if(!~i){
		if(_status){
			*_status = 0;
		}
		return NULL;
	}
	if((self->flags_a[i >> 5] >> (i&0x1F))&0x100000000ULL){
		status = 2;
	}else{
		ret = self->table_a + i*ft->size;
		memcpy(ret, key, ft->size);
		self->flags_a[i >> 5] |= 0x100000000ULL << (i&0x1F);
		self->flags_a[i >> 5] &= ~(0x1ULL << (i&0x1F));
		++self->full;
		status = 1;
	}
	if(_status){
		*_status = status;
	}
	return ret;
}

void *hash_append(hashtbl_t *self, const hashtbl_ft *ft, void *key, int *_status){
	void *ret = NULL;
	int status = 0;
	if(self->full == self->cap){
		if(!hash_ix_start(self, ft)){
			if(_status){
				*_status = 0;
			}
			return NULL;
		}
	}
	if(self->table_b){
		hash_ix_move(self, ft, self->r);
		if(self->table_b){
			if((ret = hash_get_single_b(self, ft->hash(key), ft, key))){
				status = ft->add(ret, key) ? 2 : 0;
				if(_status){
					*_status = status;
				}
				return ret;
			}
		}
	}
	uint64_t i = hash_find_slot_ap(self, ft, key);
	if(!~i){
		if(_status){
			*_status = 0;
		}
		return NULL;
	}
	ret = self->table_a + i*ft->size;
	if((self->flags_a[i >> 5] >> (i&0x1F))&0x100000000ULL){
		status = ft->add(ret, key) ? 2 : 0;
	}else{
		memcpy(ret, key, ft->size);
		self->flags_a[i >> 5] |= 0x100000000ULL << (i&0x1F);
		self->flags_a[i >> 5] &= ~(0x1ULL << (i&0x1F));
		++self->full;
		status = 1;
	}
	if(_status){
		*_status = status;
	}
	return ret;
}

inline static int hash_remove_single_a(hashtbl_t *self, uint64_t i, const hashtbl_ft *ft, const void *key){
	uint64_t a = hash_get_index(self->table_a, self->flags_a, self->len_a, i, ft, key);
	if(!~a){
		return 0;
	}
	self->flags_a[a >> 5] ^= 0x100000001ULL << (a&0x1F);
	--self->full;
	if(ft->del){
		ft->del(self->table_a + a*ft->size);
	}
	return 1;
}

inline static int hash_remove_single_b(hashtbl_t *self, uint64_t i, const hashtbl_ft *ft, const void *key){
	uint64_t b = hash_get_index(self->table_b, self->flags_b, self->len_b, i, ft, key);
	if(!~b){
		return 0;
	}
	self->flags_b[b >> 5] ^= 0x100000001ULL << (b&0x1F);
	--self->full;
	if(ft->del){
		ft->del(self->table_b + b*ft->size);
	}
	return 1;
}

inline static int hash_remove_split(hashtbl_t *self, const hashtbl_ft *ft, const void *key){
	uint64_t i = ft->hash(key);
	if(self->i << 1 < self->len_b){
		return hash_remove_single_a(self, i, ft, key) ?: hash_remove_single_b(self, i, ft, key);
	}
	return hash_remove_single_b(self, i, ft, key) ?: hash_remove_single_a(self, i, ft, key);
}

int hash_remove(hashtbl_t *self, const hashtbl_ft *ft, const void *key){
	if(!self->table_a){
		return 0;
	}
	return self->table_b ? hash_remove_split(self, ft, key) : hash_remove_single_a(self, ft->hash(key), ft, key);
}

void hash_delete(hashtbl_t *self, const hashtbl_ft *ft, void *ent){
	if(self->table_a <= ent && ent < self->table_a + self->len_a*ft->size){
		uint64_t i = (ent - self->table_a)/ft->size;
		self->flags_a[i >> 5] ^= 0x100000001ULL << (i&0x1F);
	}else{
		uint64_t i = (ent - self->table_b)/ft->size;
		self->flags_b[i >> 5] ^= 0x100000001ULL << (i&0x1F);
	}
	--self->full;
	if(ft->del){
		ft->del(ent);
	}
}

void hash_clear(hashtbl_t *self, const hashtbl_ft *ft){
	if(ft->del){
		for(uint64_t a = 0; a < self->len_a; ++a){
			uint64_t f = (self->flags_a[a >> 5] >> (a&0x1F))&0x100000001ULL;
			if(f&0x100000000ULL){
				ft->del(self->table_a + a*ft->size);
			}
		}
		for(uint64_t b = 0; b < self->len_b; ++b){
			uint64_t f = (self->flags_b[b >> 5] >> (b&0x1F))&0x100000001ULL;
			if(f&0x100000000ULL){
				ft->del(self->table_b + b*ft->size);
			}
		}
	}
	free(self->table_b);
	self->table_b = NULL;
	free(self->flags_b);
	self->flags_b = NULL;
	self->len_b = 0;
	memset(self->flags_a, 0, (((self->len_a - 1) >> 5) + 1)*sizeof(uint64_t));
	self->full = 0;
}

void hash_destroy(hashtbl_t *self, const hashtbl_ft *ft){
	hash_clear(self, ft);
	free(self->table_a);
	self->table_a = NULL;
	free(self->flags_a);
	self->flags_a = NULL;
	self->len_a = 0;
	self->cap = 0;
}

