#include <string.h>
#include <ctype.h>

#include <nut/modular_math.h>
#include <nut/factorization.h>
#include <nut/polynomial.h>

bool nut_Poly_init(nut_Poly *f, uint64_t reserve){
	reserve = reserve ?: 4;
	f->coeffs = malloc(reserve*sizeof(int64_t));
	if(!f->coeffs){
		return false;
	}
	f->len = 1;
	f->coeffs[0] = 0;
	f->cap = reserve;
	return true;
}

void nut_Poly_destroy(nut_Poly *f){
	free(f->coeffs);
	*f = (nut_Poly){};
}

int nut_Poly_cmp(const nut_Poly *a, const nut_Poly *b){
	uint64_t min_len;
	if(a->len < b->len){
		for(uint64_t i = b->len - 1; i >= a->len; --i){
			if(b->coeffs[i] > 0){
				return -1;
			}else if(b->coeffs[i] < 0){
				return 1;
			}else if(!i){
				return 0;
			}
		}
		min_len = a->len;
	}else if(a->len > b->len){
		for(uint64_t i = a->len - 1; i >= b->len; --i){
			if(a->coeffs[i] > 0){
				return 1;
			}else if(a->coeffs[i] < 0){
				return -1;
			}else if(!i){
				return 0;
			}
		}
		min_len = b->len;
	}else{
		min_len = a->len;
	}
	for(uint64_t i = min_len; i-- > 0;){
		if(a->coeffs[i] < b->coeffs[i]){
			return -1;
		}else if(a->coeffs[i] < b->coeffs[i]){
			return 1;
		}
	}
	return 0;
}

bool nut_Roots_init(nut_Roots *roots, uint64_t reserve){
	reserve = reserve ?: 4;
	roots->roots = malloc(reserve*sizeof(int64_t));
	if(!roots->roots){
		return false;
	}
	roots->len = 0;
	roots->cap = reserve;
	return true;
}

void nut_Roots_destroy(nut_Roots *roots){
	free(roots->roots);
	*roots = (nut_Roots){};
}

int64_t nut_Poly_eval_modn(const nut_Poly *f, int64_t x, int64_t n){
	if(!x){
		return f->coeffs[0];
	}
	int64_t r = 0;
	for(uint64_t i = f->len; i > 0;){
		--i;
		r = nut_i64_mod(r*x + f->coeffs[i], n);
	}
	return r;
}

bool nut_Poly_ensure_cap(nut_Poly *f, uint64_t cap){
	if(f->cap < cap){
		int64_t *tmp = realloc(f->coeffs, cap*sizeof(int64_t));
		if(!tmp){
			return false;
		}
		f->coeffs = tmp;
		f->cap = cap;
	}
	return true;
}

bool nut_Poly_zero_extend(nut_Poly *f, uint64_t len){
	if(!nut_Poly_ensure_cap(f, len)){
		return false;
	}
	for(uint64_t i = f->len; i < len; ++i){
		f->coeffs[i] = false;
	}
	if(len > f->len){
		f->len = len;
	}
	return true;
}

bool nut_Roots_ensure_cap(nut_Roots *roots, uint64_t cap){
	if(roots->cap < cap){
		int64_t *tmp = realloc(roots->roots, cap*sizeof(int64_t));
		if(!tmp){
			return false;
		}
		roots->roots = tmp;
		roots->cap = cap;
	}
	return true;
}

bool nut_Poly_add_modn(nut_Poly *h, const nut_Poly *f, const nut_Poly *g, int64_t n){
	if(!nut_Poly_ensure_cap(h, f->len > g->len ? f->len : g->len)){
		return false;
	}
	uint64_t i;
	for(i = 0; i < f->len && i < g->len; ++i){
		h->coeffs[i] = nut_i64_mod(f->coeffs[i] + g->coeffs[i], n);
	}
	for(; i < f->len; ++i){
		h->coeffs[i] = f->coeffs[i];
	}
	for(; i < g->len; ++i){
		h->coeffs[i] = g->coeffs[i];
	}
	int same_len = f->len == g->len;
	h->len = f->len > g->len ? f->len : g->len;
	if(same_len){
		nut_Poly_normalize(h);
	}
	return true;
}

bool nut_Poly_sub_modn(nut_Poly *h, const nut_Poly *f, const nut_Poly *g, int64_t n){
	if(!nut_Poly_ensure_cap(h, f->len > g->len ? f->len : g->len)){
		return false;
	}
	uint64_t i;
	for(i = 0; i < f->len && i < g->len; ++i){
		h->coeffs[i] = nut_i64_mod(f->coeffs[i] - g->coeffs[i], n);
	}
	for(; i < f->len; ++i){
		h->coeffs[i] = f->coeffs[i];
	}
	for(; i < g->len; ++i){
		h->coeffs[i] = g->coeffs[i] ? n - g->coeffs[i] : 0;
	}
	int same_len = f->len == g->len;
	h->len = f->len > g->len ? f->len : g->len;
	if(same_len){
		nut_Poly_normalize(h);
	}
	return true;
}

bool nut_Poly_dot_modn(nut_Poly *h, const nut_Poly *f, const nut_Poly *g, int64_t n){
	if(!nut_Poly_ensure_cap(h, f->len < g->len ? f->len : g->len)){
		return false;
	}
	uint64_t i;
	for(i = 0; i < f->len && i < g->len; ++i){
		h->coeffs[i] = nut_i64_mod(f->coeffs[i]*g->coeffs[i], n);
	}
	h->len = f->len < g->len ? f->len : g->len;
	nut_Poly_normalize(h);
	return true;
}

bool nut_Poly_copy(nut_Poly *restrict g, const nut_Poly *restrict f){
	if(!nut_Poly_ensure_cap(g, f->len)){
		return false;
	}
	memcpy(g->coeffs, f->coeffs, f->len*sizeof(int64_t));
	g->len = f->len;
	return true;
}

bool nut_Poly_setconst(nut_Poly *f, int64_t c){
	if(!nut_Poly_ensure_cap(f, 1)){
		return false;
	}
	f->coeffs[0] = c;
	f->len = 1;
	return true;
}

bool nut_Poly_scale_modn(nut_Poly *g, const nut_Poly *f, int64_t a, int64_t n){
	if(!a){
		return nut_Poly_setconst(g, 0);
	}
	if(!nut_Poly_ensure_cap(g, f->len)){
		return false;
	}
	for(uint64_t i = 0; i < f->len; ++i){
		g->coeffs[i] = nut_i64_mod(a*f->coeffs[i], n);
	}
	g->len = f->len;
	nut_Poly_normalize(g);
	return true;
}

bool nut_Poly_mul_modn(nut_Poly *restrict h, const nut_Poly *f, const nut_Poly *g, int64_t n){
	if(f->len == 1){
		return nut_Poly_scale_modn(h, g, f->coeffs[0], n);
	}else if(g->len == 1){
		return nut_Poly_scale_modn(h, f, g->coeffs[0], n);
	}
	if(!nut_Poly_ensure_cap(h, f->len + g->len - 1)){
		return false;
	}
	for(uint64_t k = 0; k < f->len + g->len - 1; ++k){
		h->coeffs[k] = 0;
		for(uint64_t i = k >= g->len ? k - g->len + 1 : 0; i < f->len && i <= k; ++i){
			h->coeffs[k] = nut_i64_mod(h->coeffs[k] + f->coeffs[i]*g->coeffs[k - i], n);
		}
	}
	h->len = f->len + g->len - 1;
	nut_Poly_normalize(h);
	return true;
}

bool nut_Poly_pow_modn(nut_Poly *restrict g, const nut_Poly *f, uint64_t e, int64_t n, uint64_t cn, nut_Poly tmps[restrict static 2]){
	//tmps: st, rt
	if(f->len > cn + 1){
		return false;
	}else if(!e){
		return nut_Poly_setconst(g, 1); // TODO: deal with 0**0 case
	}else if(e == 1){
		return nut_Poly_copy(g, f);
	}
	e = 1 + (e - 1)%cn;
	if(f->len == 1){
		return nut_Poly_setconst(g, nut_u64_pow(f->coeffs[0], e));
	}else if(!nut_Poly_copy(tmps + 0, f) || !nut_Poly_ensure_cap(tmps + 0, 2*cn + 1) || !nut_Poly_ensure_cap(tmps + 1, 2*cn + 1)){
		return false;
	}
	nut_Poly *t = g, *s = tmps + 0, *r = tmps + 1;
	while(e%2 == 0){
		nut_Poly_mul_modn(t, s, s, n);
		nut_Poly_normalize_exps_modn(t, cn);
		{
			void *tmp = t;
			t = s;
			s = tmp;
		}//s = s*s
		e >>= 1;
	}
	if(!nut_Poly_copy(r, s)){
		return false;
	}
	while((e >>= 1)){
		nut_Poly_mul_modn(t, s, s, n);
		nut_Poly_normalize_exps_modn(t, cn);
		{
			void *tmp = t;
			t = s;
			s = tmp;
		}//s = s*s
		if(e%2){
			nut_Poly_mul_modn(t, r, s, n);
			nut_Poly_normalize_exps_modn(t, cn);
			{
				void *tmp = t;
				t = r;
				r = tmp;
			}//r = r*s
		}
	}
	if(r != g){
		if(!nut_Poly_copy(g, r)){
			return 0;
		}
	}
	return true;
}

bool nut_Poly_pow_modn_tmptmp(nut_Poly *restrict g, const nut_Poly *f, uint64_t e, int64_t n, uint64_t cn){
	nut_Poly tmps[2] = {};
	bool status = true;
	for(uint64_t i = 0; status && i < 2; ++i){
		status = nut_Poly_init(tmps + i, 2*cn + 1);
	}
	if(status){
		status = nut_Poly_pow_modn(g, f, e, n, cn, tmps);
	}
	for(uint64_t i = 0; i < 2; ++i){
		if(tmps[i].cap){
			nut_Poly_destroy(tmps + i);
		}
	}
	return status;
}

bool nut_Poly_compose_modn(nut_Poly *restrict h, const nut_Poly *f, const nut_Poly *g, int64_t n, uint64_t cn, nut_Poly tmps[restrict static 2]){
	if(f->len == 1){
		return nut_Poly_setconst(h, f->coeffs[0]);
	}else if(g->len == 1){
		return nut_Poly_setconst(h, nut_Poly_eval_modn(f, g->coeffs[0], n));
	}
	uint64_t h_len = (f->len - 1)*(g->len - 1) + 1;
	if(h_len > cn + 1){
		h_len = cn + 1;
	}
	if(!nut_Poly_ensure_cap(h, h_len) || !nut_Poly_setconst(h, f->coeffs[0])){
		return false;
	}
	nut_Poly *t = tmps + 0, *p = tmps + 1;
	if(!nut_Poly_copy(p, g) || !nut_Poly_scale_modn(t, g, f->coeffs[1], n) || !nut_Poly_add_modn(h, h, t, n)){
		return false;
	}
	for(uint64_t e = 2; e <= cn && e < f->len; ++e){
		if(!nut_Poly_mul_modn(t, p, g, n)){
			return false;
		}
		nut_Poly_normalize_exps_modn(t, cn);
		void *tmp = t;
		t = p;
		p = tmp;
		if(!nut_Poly_scale_modn(t, p, f->coeffs[e], n) || !nut_Poly_add_modn(h, h, t, n)){
			return false;
		}
	}
	return true;
}

bool nut_Poly_compose_modn_tmptmp(nut_Poly *restrict h, const nut_Poly *f, const nut_Poly *g, int64_t n, uint64_t cn){
	nut_Poly tmps[2] = {};
	bool status = true;
	for(uint64_t i = 0; status && i < 2; ++i){
		status = nut_Poly_init(tmps + i, 2*cn + 1);
	}
	if(status){
		status = nut_Poly_compose_modn(h, f, g, n, cn, tmps);
	}
	for(uint64_t i = 0; i < 2; ++i){
		if(tmps[i].cap){
			nut_Poly_destroy(tmps + i);
		}
	}
	return status;
}

bool nut_Poly_quotrem_modn(nut_Poly *restrict q, nut_Poly *restrict r, const nut_Poly *restrict f, const nut_Poly *restrict g, int64_t n){
	int64_t a, d = nut_i64_egcd(g->coeffs[g->len - 1], n, &a, NULL);
	if(d != 1){
		return false;//TODO: set divide by zero error, or set remainder (?)
	}
	if(g->len == 1){//dividing by a scalar TODO: set divide by zero error if scalar is 0
		return g->coeffs[0] && nut_Poly_setconst(r, 0) && nut_Poly_scale_modn(q, f, a, n);
	}
	if(f->len < g->len){//dividing by a polynomial with higher degree
		return nut_Poly_setconst(q, 0) && nut_Poly_copy(r, f);
	}
	
	//begin extended synthetic division
	//compute max length of quotient and remainder and extend their buffers if need be
	q->len = f->len - g->len + 1;
	r->len = g->len - 1;
	if(!nut_Poly_ensure_cap(q, q->len) || !nut_Poly_ensure_cap(r, r->len)){
		return false;
	}
	
	//initialize column sums/coeffs of results
	memcpy(r->coeffs, f->coeffs, r->len*sizeof(int64_t));
	memcpy(q->coeffs, f->coeffs + r->len, q->len*sizeof(int64_t));
	
	//loop over quotient columns (coefficients in reverse order) which were initialized
	//to the first q->len dividend coefficients
	for(uint64_t k = q->len; k > 0;){
		--k;
		//finish the column by dividing the sum by the leading coefficient of the divisor
		//for monic divisors (aka most of the time) a will simply be 1 so we may optimize this
		q->coeffs[k] = nut_i64_mod(q->coeffs[k]*a, n);
		//subtract the adjusted column sum times the some of the divisor coefficients from the
		//remainder coefficients.  the remainder coefficients were initialized to the last r->len
		//dividend coefficients.  we start q->len - k columns after the current column we just finished.
		//if k == 0 we should go from coefficient 0 in the remainder to coefficient r->len - 1
		//so in general we should go from k to r->len - 1 (this can easily be an empty interval)
		for(uint64_t i = 0, j = k; j < r->len; ++i, ++j){
			r->coeffs[j] = nut_i64_mod(r->coeffs[j] - q->coeffs[k]*g->coeffs[i], n);
		}
		//j goes from max(k - g->len + 1, 0) to k - 1 (both inclusive)
		//i ends at g->len - 2 (inclusive) and should go through the same number of values
		//so i starts at g->len - 2 - (k - 1 - j) = g->len + j - k - 1
		for(uint64_t j = k > g->len - 1 ? k - g->len + 1 : 0, i = g->len + j - k - 1; i < g->len - 1; ++i, ++j){
			q->coeffs[j] = nut_i64_mod(q->coeffs[j] - q->coeffs[k]*g->coeffs[i], n);
		}
	}
	nut_Poly_normalize(r);
	return true;
}

void nut_Poly_normalize(nut_Poly *f){
	while(f->len > 1 && !f->coeffs[f->len - 1]){
		--f->len;
	}
}

void nut_Poly_normalize_modn(nut_Poly *f, int64_t n, bool use_negatives){
	int64_t offset = use_negatives ? (1-n)/2 : 0;
	for(uint64_t i = 0; i < f->len; ++i){
		f->coeffs[i] = offset + nut_i64_mod(f->coeffs[i] - offset, n);
	}
	nut_Poly_normalize(f);
}

void nut_Poly_normalize_exps_modn(nut_Poly *f, uint64_t cn){
	for(uint64_t i = f->len - 1; i > cn; --i){
		uint64_t j = 1 + (i - 1)%cn;
		f->coeffs[j] += f->coeffs[i];
	}
	if(f->len > cn + 1){
		f->len = cn + 1;
	}
	nut_Poly_normalize(f);
}

int nut_Poly_fprint(FILE *file, const nut_Poly *f, const char *vname, const char *add, const char *sub, const char *pow, bool descending){
	int res = 0;
	for(uint64_t i = 0; i < f->len; ++i){
		uint64_t j = descending ? f->len - 1 - i : i;
		int64_t coeff = f->coeffs[j];
		if(!coeff){
			continue;
		}
		if(res){
			res += fprintf(file, "%s", coeff > 0 ? add : sub);
			coeff = coeff < 0 ? -coeff : coeff;
			if(coeff != 1 || !j){
				res += fprintf(file, "%"PRId64, coeff);
			}
		}else{
			if(coeff != 1 || !j){
				if(coeff == -1 && j){
					res += fprintf(file, "-");
				}else{
					res += fprintf(file, "%"PRId64, coeff);
				}
			}
		}
		if(j){
			res += fprintf(file, "%s", vname);
		}
		if(j > 1){
			res += fprintf(file, "%s%"PRIu64, pow, j);
		}
	}
	return res;
}

static inline void skip_whitespace(const char *restrict *restrict _str){
	while(isspace(**_str)){
		++*_str;
	}
}

static inline int parse_monomial(nut_Poly *restrict f, const char *restrict *restrict _str){
	const char *str = *_str;
	int64_t sign = 1, coeff = 0;
	uint64_t x = 0;
	while(*str == '+' || *str == '-' || isspace(*str)){
		if(*str == '-'){
			sign = -sign;
		}
		++str;
	}
	bool need_vpow = true;
	while(isdigit(*str)){
		coeff = 10*coeff + *str - '0';
		need_vpow = false;
		++str;
	}
	skip_whitespace(&str);
	if(need_vpow){
		coeff = 1;
	}
	if(!need_vpow && *str == '*'){
		++str;
		need_vpow = true;
		skip_whitespace(&str);
	}
	bool have_vpow = false;
	if(strncmp(str, "mod", 3) || !isspace(str[3])){
		while(isalpha(*str)){
			have_vpow = true;
			++str;
		}
	}
	if(have_vpow){
		skip_whitespace(&str);
		if(*str == '^' || !strncmp(str, "**", 2)){
			str += *str == '^' ? 1 : 2;
			skip_whitespace(&str);
			bool have_pow = false;
			while(isdigit(*str)){
				x = 10*x + *str - '0';
				have_pow = true;
				++str;
			}
			if(!have_pow){
				return 0;
			}
		}else{
			x = 1;
		}
	}
	skip_whitespace(&str);
	if(!coeff){
		*_str = str;
		return 1;
	}
	if(!nut_Poly_zero_extend(f, x + 1)){
		return 0;
	}
	f->coeffs[x] += sign*coeff;
	*_str = str;
	return 1;
}

int nut_Poly_parse(nut_Poly *restrict f, int64_t *restrict n, const char *restrict str, const char *restrict *restrict end){
	if(!nut_Poly_setconst(f, 0) || !parse_monomial(f, &str)){
		return 0;
	}
	while(*str == '+' || *str == '-'){
		if(!parse_monomial(f, &str)){
			if(end){
				skip_whitespace(&str);
				*end = str;
			}
			nut_Poly_normalize(f);
			return 1;
		}
	}
	if(!strncmp(str, "mod", 3) && isspace(str[3])){
		str += 4;
		skip_whitespace(&str);
		int64_t _n = 0;
		bool have_mod = false;
		while(isdigit(*str)){
			_n = 10*_n + *str - '0';
			have_mod = true;
			++str;
		}
		if(have_mod){
			*n = _n;
			if(end){
				skip_whitespace(&str);
				*end = str;
			}
			nut_Poly_normalize_modn(f, *n, 0);
			return 2;
		}else{
			str -= 4;
		}
	}
	if(end){
		skip_whitespace(&str);
		*end = str;
	}
	nut_Poly_normalize(f);
	return 1;
}

bool nut_Poly_rand_modn(nut_Poly *f, uint64_t max_len, int64_t n){
	if(!max_len){
		return nut_Poly_setconst(f, 0);
	}
	if(!nut_Poly_ensure_cap(f, max_len)){
		return false;
	}
	for(uint64_t i = 0; i < max_len; ++i){
		f->coeffs[i] = nut_u64_rand(0, n);
	}
	f->len = max_len;
	nut_Poly_normalize(f);
	return true;
}

bool nut_Poly_gcd_modn(nut_Poly *restrict d, const nut_Poly *restrict f, const nut_Poly *restrict g, int64_t n, nut_Poly tmps[restrict static 3]){
	//tmps: qt, r0t, r1t
	nut_Poly *tmp, *r0, *r1, *r2;
	bool status = true;
	if(g->len > f->len){//Ensure the degree of f is >= the degree of g
		const nut_Poly *tmp = f;//Shadow the other tmp with a const version
		f = g;
		g = tmp;
	}
	if(g->len == 1){//If one of the inputs is a constant, we either have the gcd of something and a unit or something and zero
		status = g->coeffs[0] ? nut_Poly_setconst(d, 1) : nut_Poly_copy(d, f);
		goto CLEANUP;
	}
	
	r0 = d, r1 = tmps + 1, r2 = tmps + 2;
	//Unroll first two remainder calculations to prevent copying
	if(!nut_Poly_quotrem_modn(tmps + 0, r0, f, g, n)){
		status = false;
		goto CLEANUP;
	}
	if(r0->len == 1 && !r0->coeffs[0]){
		status = nut_Poly_copy(d, g);
		goto CLEANUP;
	}
	
	if(!nut_Poly_quotrem_modn(tmps + 0, r1, g, r0, n)){
		status = false;
		goto CLEANUP;
	}
	if(r1->len == 1 && !r1->coeffs[0]){
		goto CLEANUP;// r0 == d in this case, no need to copy
	}
	
	//Euclidean algorithm: take remainders until we reach 0, then the last nonzero remainder is the gcd
	while(1){
		if(!nut_Poly_quotrem_modn(tmps + 0, r2, r0, r1, n)){
			status = false;
			goto CLEANUP;
		}
		if(r2->len == 1 && !r2->coeffs[0]){
			if(r1 != d){//Make sure the result is in the output and not a temporary
				status = nut_Poly_copy(d, r1);
			}
			goto CLEANUP;
		}
		tmp = r0;
		r0 = r1;
		r1 = r2;
		r2 = tmp;
	}
	
	CLEANUP:;
	if(status){//If the gcd was found, make it monic
		int64_t a = d->coeffs[d->len - 1];
		if(a == n - 1){//Inverse of -1 is always -1
			nut_Poly_scale_modn(d, d, a, n);
		}else if(a > 1){//If the leading coefficient is 0 or 1 the gcd is already monic
			int64_t c, g = nut_i64_egcd(a, n, &c, NULL);//This g shadows the input
			if(g != 1){//The leading coefficient cannot be inverted mod n
				return false;
			}
			nut_Poly_scale_modn(d, d, c, n);
		}
	}
	return status;
}

bool nut_Poly_gcd_modn_tmptmp(nut_Poly *restrict d, const nut_Poly *restrict f, const nut_Poly *restrict g, int64_t n){
	nut_Poly tmps[3] = {};
	uint64_t min_len, quot_len;
	if(f->len <= g->len){
		min_len = f->len;
		quot_len = g->len - f->len + 1;
	}else{
		min_len = g->len;
		quot_len = f->len - g->len + 1;
	}
	bool status = nut_Poly_init(tmps + 0, quot_len) &&
	         nut_Poly_init(tmps + 1, min_len) &&
	         nut_Poly_init(tmps + 2, min_len) &&
	         nut_Poly_gcd_modn(d, f, g, n, tmps);
	for(uint64_t i = 0; i < 3; ++i){
		if(tmps[i].cap){
			nut_Poly_destroy(tmps + i);
		}
	}
	return status;
}

bool nut_Poly_powmod_modn(nut_Poly *restrict h, const nut_Poly *restrict f, uint64_t e, const nut_Poly *restrict g, int64_t n, nut_Poly tmps[restrict static 3]){
	//tmps: qt, st, rt
	if(g->len <= f->len){
		if(!nut_Poly_quotrem_modn(tmps + 0, tmps + 1, f, g, n)){
			return false;
		}
	}else if(!nut_Poly_copy(tmps + 1, f)){
		return false;
	}
	if(tmps[1].len == 1){
		if(!tmps[1].coeffs[0]){
			return nut_Poly_setconst(h, 0);
		}else if(!e){
			return nut_Poly_setconst(h, 1);
		}
	}
	if(!nut_Poly_ensure_cap(tmps + 1, 2*g->len - 3) || !nut_Poly_ensure_cap(tmps + 2, 2*g->len - 3) || !nut_Poly_ensure_cap(h, 2*g->len - 3)){
		return false;
	}
	
	nut_Poly *t = h, *s = tmps + 1, *r = tmps + 2;
	while(e%2 == 0){
		nut_Poly_mul_modn(t, s, s, n);
		if(!nut_Poly_quotrem_modn(tmps + 0, s, t, g, n)){
			return false;
		}//s = s*s%g
		e >>= 1;
	}
	if(!nut_Poly_copy(r, s)){
		return false;
	}
	while((e >>= 1)){
		nut_Poly_mul_modn(t, s, s, n);
		if(!nut_Poly_quotrem_modn(tmps + 0, s, t, g, n)){
			return false;
		}//s = s*s%g
		if(e%2){
			nut_Poly_mul_modn(t, r, s, n);
			if(!nut_Poly_quotrem_modn(tmps + 0, r, t, g, n)){
				return false;
			}//r = r*s%g
		}
	}
	
	if(r != h){// TODO: this check is always true, fix buffer juggling so we never need to copy
		if(!nut_Poly_copy(h, r)){
			return false;
		}
	}
	return true;
}

bool nut_Poly_powmod_modn_tmptmp(nut_Poly *restrict h, const nut_Poly *restrict f, uint64_t e, const nut_Poly *restrict g, int64_t n){
	nut_Poly tmps[3] = {};
	bool status = true;
	for(uint64_t i = 0; status && i < 3; ++i){
		status = nut_Poly_init(tmps + i, 2*g->len - 1);
	}
	if(status){
		status = nut_Poly_powmod_modn(h, f, e, g, n, tmps);
	}
	for(uint64_t i = 0; i < 3; ++i){
		if(tmps[i].cap){
			nut_Poly_destroy(tmps + i);
		}
	}
	return status;
}

bool nut_Poly_factors_d_modn(nut_Poly *restrict f_d, const nut_Poly *restrict f, uint64_t d, int64_t n, nut_Poly tmps[restrict static 4]){
	//tmps: xt, qt, st, rt
	if(!nut_Poly_ensure_cap(f_d, f->len)){
		return false;
	}
	f_d->coeffs[0] = 0;
	f_d->coeffs[1] = 1;
	f_d->len = 2;
	return nut_Poly_powmod_modn(tmps + 0, f_d, nut_u64_pow(n, d), f, n, tmps + 1) &&
		nut_Poly_sub_modn(tmps + 0, tmps + 0, f_d, n) &&
		nut_Poly_gcd_modn(f_d, tmps + 0, f, n, tmps + 1);
}

bool nut_Poly_factor1_modn(nut_Poly *restrict g, const nut_Poly *restrict f, uint64_t d, int64_t n, nut_Poly tmps[restrict static 4]){
	//tmps: xt, qt, st, rt
	while(1){
		if(!nut_Poly_rand_modn(tmps + 0, f->len - 1, n) || !nut_Poly_gcd_modn(g, tmps + 0, f, n, tmps + 1)){
			return false;
		}
		if(g->len > 1 && g->len < f->len){
			return true;
		}
		if(
			!nut_Poly_powmod_modn(g, tmps + 0, (nut_u64_pow(n, d)-1)/2, f, n, tmps + 1) ||
			!nut_Poly_setconst(tmps + 0, 1) ||
			!nut_Poly_add_modn(tmps + 0, g, tmps + 0, n) ||
			!nut_Poly_gcd_modn(g, tmps + 0, f, n, tmps + 1)
		){
			return false;
		}
		if(g->len > 1 && g->len < f->len){
			return true;
		}
	}
}

[[gnu::nonnull(1, 3)]]
NUT_ATTR_ACCESS(read_write, 1)
static inline bool roots_polyn_modn_rec(nut_Roots *restrict roots, int64_t n, nut_Poly tmps[restrict static 6]){
	//tmps: gt, ft, xt, qt, st, rt
	while(1){
		//fprintf(stderr, "\e[1;33mFactoring (");
		//nut_Poly_fprint(stderr, tmps + 1, "x", " + ", " - ", "**", 0);
		//fprintf(stderr, ") mod %"PRId64"\e[0m\n", n);
		if(tmps[1].len == 2){//linear factor
			//fprintf(stderr, "\e[1;33m polynomial is linear\e[0m\n");
			if(tmps[1].coeffs[0]){//nonzero root
				if(tmps[1].coeffs[1] == 1){//monic
					roots->roots[roots->len++] = n - tmps[1].coeffs[0];
				}else{//not monic
					int64_t a;
					nut_i64_egcd(tmps[1].coeffs[0], n, &a, NULL);
					roots->roots[roots->len++] = nut_i64_mod(-tmps[1].coeffs[0]*a, n);
				}
			}else{//zero root
				roots->roots[roots->len++] = 0;
			}
			return true;
		}else if(tmps[1].len == 3){//quadratic factor
			//fprintf(stderr, "\e[1;33m polynomial is quadratic\e[0m\n");
			if(tmps[1].coeffs[2] != 1){//make monic
				int64_t a;
				nut_i64_egcd(tmps[1].coeffs[2], n, &a, NULL);
				nut_Poly_scale_modn(tmps + 1, tmps + 1, a, n);
			}
			int64_t c = tmps[1].coeffs[0];
			int64_t b = tmps[1].coeffs[1];
			int64_t r = nut_i64_sqrt_mod(nut_i64_mod(b*b - 4*c, n), n);
			roots->roots[roots->len++] = nut_i64_mod((n+1)/2*(-b + r), n);
			roots->roots[roots->len++] = nut_i64_mod((n+1)/2*(-b - r), n);
			return true;
		}
		if(!nut_Poly_factor1_modn(tmps + 0, tmps + 1, 1, n, tmps + 2) || !nut_Poly_quotrem_modn(tmps + 3, tmps + 5, tmps + 1, tmps + 0, n)){
			return false;
		}
		//tmps[0] and tmps[3] hold the nontrivial factors we found
		//we have to do a recursive call on the factor with smaller degree
		//if it is linear or quadratic, we do not need to copy the larger factor
		//otherwise we do
		nut_Poly tmp;
		tmp = tmps[0];
		tmps[0] = tmps[1];
		if(tmp.len <= tmps[3].len){
			tmps[1] = tmp;
		}else{
			tmps[1] = tmps[3];
			tmps[3] = tmp;
		}
		//now, tmps[1] and tmps[3] hold the nontrivial factors, with tmps[1] having smallest degree
		//fprintf(stderr, "\e[1;33m got factors (");
		//nut_Poly_fprint(stderr, tmps + 1, "x", " + ", " - ", "**", 0);
		//fprintf(stderr, ")(");
		//nut_Poly_fprint(stderr, tmps + 3, "x", " + ", " - ", "**", 0);
		//fprintf(stderr, ")\e[0m\n");
		bool have_small_factor = tmps[1].len <= 3;
		if(!have_small_factor){
			if(!nut_Poly_init(&tmp, tmps[3].len) || !nut_Poly_copy(&tmp, tmps + 3)){
				return false;
			}
		}
		if(!roots_polyn_modn_rec(roots, n, tmps)){
			return false;
		}
		if(!have_small_factor){
			if(!nut_Poly_copy(tmps + 3, &tmp)){
				return false;
			}
			nut_Poly_destroy(&tmp);
		}
		tmp = tmps[3];
		tmps[3] = tmps[1];
		tmps[1] = tmp;
	}
}

bool nut_Poly_roots_modn(const nut_Poly *restrict f, int64_t n, nut_Roots *restrict roots, nut_Poly tmps[restrict static 6]){
	//tmps: gt, ft, xt, qt, st, rt
	if(!nut_Poly_factors_d_modn(tmps + 1, f, 1, n, tmps + 2)){
		return false;
	}
	roots->len = 0;
	if(tmps[1].len == 1){
		return true;
	}
	return nut_Roots_ensure_cap(roots, tmps[1].len - 1) && roots_polyn_modn_rec(roots, n, tmps);
}

bool nut_Poly_roots_modn_tmptmp(const nut_Poly *restrict f, int64_t n, nut_Roots *restrict roots){
	nut_Poly tmps[6] = {};
	bool status = true;
	for(uint64_t i = 0; status && i < 6; ++i){
		status = nut_Poly_init(tmps + i, f->len);
	}
	if(status){
		status = nut_Poly_roots_modn(f, n, roots, tmps);
	}
	for(uint64_t i = 0; i < 6; ++i){
		if(tmps[i].cap){
			nut_Poly_destroy(tmps + i);
		}
	}
	return status;
}

