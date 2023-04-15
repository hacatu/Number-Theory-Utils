#include <string.h>
#include <ctype.h>

#include <nut/factorization.h>
#include <nut/polynomial.h>

int init_poly(poly_t *f, uint64_t reserve){
	reserve = reserve ?: 4;
	f->coeffs = malloc(reserve*sizeof(int64_t));
	if(!f->coeffs){
		return 0;
	}
	f->len = 1;
	f->coeffs[0] = 0;
	f->cap = reserve;
	return 1;
}

void destroy_poly(poly_t *f){
	free(f->coeffs);
	*f = (poly_t){};
}

int cmp_polys(const poly_t *a, const poly_t *b){
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

int init_poly_roots(poly_roots_t *roots, uint64_t reserve){
	reserve = reserve ?: 4;
	roots->roots = malloc(reserve*sizeof(int64_t));
	if(!roots->roots){
		return 0;
	}
	roots->len = 0;
	roots->cap = reserve;
	return 1;
}

void destroy_poly_roots(poly_roots_t *roots){
	free(roots->roots);
	*roots = (poly_roots_t){};
}

int64_t eval_poly_modn(const poly_t *f, int64_t x, int64_t n){
	if(!x){
		return f->coeffs[0];
	}
	int64_t r = 0;
	for(uint64_t i = f->len; i-- > 0;){
		r = mod(r*x + f->coeffs[i], n);
	}
	return r;
}

int ensure_poly_cap(poly_t *f, uint64_t cap){
	if(f->cap < cap){
		int64_t *tmp = realloc(f->coeffs, cap*sizeof(int64_t));
		if(!tmp){
			return 0;
		}
		f->coeffs = tmp;
		f->cap = cap;
	}
	return 1;
}

int zero_extend_poly(poly_t *f, uint64_t len){
	if(!ensure_poly_cap(f, len)){
		return 0;
	}
	for(uint64_t i = f->len; i < len; ++i){
		f->coeffs[i] = 0;
	}
	if(len > f->len){
		f->len = len;
	}
	return 1;
}

int ensure_poly_roots_cap(poly_roots_t *roots, uint64_t cap){
	if(roots->cap < cap){
		int64_t *tmp = realloc(roots->roots, cap*sizeof(int64_t));
		if(!tmp){
			return 0;
		}
		roots->roots = tmp;
		roots->cap = cap;
	}
	return 1;
}

int add_poly_modn(poly_t *h, const poly_t *f, const poly_t *g, int64_t n){
	if(!ensure_poly_cap(h, f->len > g->len ? f->len : g->len)){
		return 0;
	}
	uint64_t i;
	for(i = 0; i < f->len && i < g->len; ++i){
		h->coeffs[i] = mod(f->coeffs[i] + g->coeffs[i], n);
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
		normalize_poly(h);
	}
	return 1;
}

int sub_poly_modn(poly_t *h, const poly_t *f, const poly_t *g, int64_t n){
	if(!ensure_poly_cap(h, f->len > g->len ? f->len : g->len)){
		return 0;
	}
	uint64_t i;
	for(i = 0; i < f->len && i < g->len; ++i){
		h->coeffs[i] = mod(f->coeffs[i] - g->coeffs[i], n);
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
		normalize_poly(h);
	}
	return 1;
}

int dot_poly_modn(poly_t *h, const poly_t *f, const poly_t *g, int64_t n){
	if(!ensure_poly_cap(h, f->len < g->len ? f->len : g->len)){
		return 0;
	}
	uint64_t i;
	for(i = 0; i < f->len && i < g->len; ++i){
		h->coeffs[i] = mod(f->coeffs[i]*g->coeffs[i], n);
	}
	h->len = f->len < g->len ? f->len : g->len;
	normalize_poly(h);
	return 1;
}

int copy_poly(poly_t *g, const poly_t *f){
	if(!ensure_poly_cap(g, f->len)){
		return 0;
	}
	memcpy(g->coeffs, f->coeffs, f->len*sizeof(int64_t));
	g->len = f->len;
	return 1;
}

int const_poly(poly_t *f, int64_t c){
	if(!ensure_poly_cap(f, 1)){
		return 0;
	}
	f->coeffs[0] = c;
	f->len = 1;
	return 1;
}

int scale_poly_modn(poly_t *g, const poly_t *f, int64_t a, int64_t n){
	if(!a){
		return const_poly(g, 0);
	}
	if(!ensure_poly_cap(g, f->len)){
		return 0;
	}
	for(uint64_t i = 0; i < f->len; ++i){
		g->coeffs[i] = mod(a*f->coeffs[i], n);
	}
	g->len = f->len;
	normalize_poly(g);
	return 1;
}

int mul_poly_modn(poly_t *restrict h, const poly_t *f, const poly_t *g, int64_t n){
	if(f->len == 1){
		return scale_poly_modn(h, g, f->coeffs[0], n);
	}else if(g->len == 1){
		return scale_poly_modn(h, f, g->coeffs[0], n);
	}
	if(!ensure_poly_cap(h, f->len + g->len - 1)){
		return 0;
	}
	for(uint64_t k = 0; k < f->len + g->len - 1; ++k){
		h->coeffs[k] = 0;
		for(uint64_t i = k >= g->len ? k - g->len + 1 : 0; i < f->len && i <= k; ++i){
			h->coeffs[k] = mod(h->coeffs[k] + f->coeffs[i]*g->coeffs[k - i], n);
		}
	}
	h->len = f->len + g->len - 1;
	normalize_poly(h);
	return 1;
}

int pow_poly_modn(poly_t *restrict g, const poly_t *f, uint64_t e, int64_t n, uint64_t cn, poly_t tmps[static 2]){
	//tmps: st, rt
	if(f->len > cn + 1){
		return 0;
	}else if(!e){
		return const_poly(g, 1); // TODO: deal with 0**0 case
	}else if(e == 1){
		return copy_poly(g, f);
	}
	e = 1 + (e - 1)%cn;
	if(f->len == 1){
		return const_poly(g, pow_u64(f->coeffs[0], e));
	}else if(!copy_poly(tmps + 0, f) || !ensure_poly_cap(tmps + 0, 2*cn + 1) || !ensure_poly_cap(tmps + 1, 2*cn + 1)){
		return 0;
	}
	poly_t *t = g, *s = tmps + 0, *r = tmps + 1;
	while(e%2 == 0){
		mul_poly_modn(t, s, s, n);
		normalize_exponents_modn(t, cn);
		{
			void *tmp = t;
			t = s;
			s = tmp;
		}//s = s*s
		e >>= 1;
	}
	copy_poly(r, s);
	while((e >>= 1)){
		mul_poly_modn(t, s, s, n);
		normalize_exponents_modn(t, cn);
		{
			void *tmp = t;
			t = s;
			s = tmp;
		}//s = s*s
		if(e%2){
			mul_poly_modn(t, r, s, n);
			normalize_exponents_modn(t, cn);
			{
				void *tmp = t;
				t = r;
				r = tmp;
			}//r = r*s
		}
	}
	if(r != g){
		if(!copy_poly(g, r)){
			return 0;
		}
	}
	return 1;
}

int pow_poly_modn_tmptmp(poly_t *restrict g, const poly_t *f, uint64_t e, int64_t n, uint64_t cn){
	poly_t tmps[2] = {};
	int status = 1;
	for(uint64_t i = 0; status && i < 2; ++i){
		status = init_poly(tmps + i, 2*cn + 1);
	}
	if(status){
		status = pow_poly_modn(g, f, e, n, cn, tmps);
	}
	for(uint64_t i = 0; i < 2; ++i){
		if(tmps[i].cap){
			destroy_poly(tmps + i);
		}
	}
	return status;
}

int compose_poly_modn(poly_t *restrict h, const poly_t *f, const poly_t *g, int64_t n, uint64_t cn, poly_t tmps[static 2]){
	if(f->len == 1){
		return const_poly(h, f->coeffs[0]);
	}else if(g->len == 1){
		return const_poly(h, eval_poly_modn(f, g->coeffs[0], n));
	}
	uint64_t h_len = (f->len - 1)*(g->len - 1) + 1;
	if(h_len > cn + 1){
		h_len = cn + 1;
	}
	if(!ensure_poly_cap(h, h_len) || !const_poly(h, f->coeffs[0])){
		return 0;
	}
	poly_t *t = tmps + 0, *p = tmps + 1;
	if(!copy_poly(p, g) || !scale_poly_modn(t, g, f->coeffs[1], n) || !add_poly_modn(h, h, t, n)){
		return 0;
	}
	for(uint64_t e = 2; e <= cn && e < f->len; ++e){
		if(!mul_poly_modn(t, p, g, n)){
			return 0;
		}
		normalize_exponents_modn(t, cn);
		void *tmp = t;
		t = p;
		p = tmp;
		if(!scale_poly_modn(t, p, f->coeffs[e], n) || !add_poly_modn(h, h, t, n)){
			return 0;
		}
	}
	return 1;
}

int compose_poly_modn_tmptmp(poly_t *restrict h, const poly_t *f, const poly_t *g, int64_t n, uint64_t cn){
	poly_t tmps[2] = {};
	int status = 1;
	for(uint64_t i = 0; status && i < 2; ++i){
		status = init_poly(tmps + i, 2*cn + 1);
	}
	if(status){
		status = compose_poly_modn(h, f, g, n, cn, tmps);
	}
	for(uint64_t i = 0; i < 2; ++i){
		if(tmps[i].cap){
			destroy_poly(tmps + i);
		}
	}
	return status;
}

int quotrem_poly_modn(poly_t *restrict q, poly_t *restrict r, const poly_t *restrict f, const poly_t *restrict g, int64_t n){
	int64_t a, d = egcd(g->coeffs[g->len - 1], n, &a, NULL);
	if(d != 1){
		return 0;//TODO: set divide by zero error, or set remainder (?)
	}
	if(g->len == 1){//dividing by a scalar
		if(!g->coeffs[0]){
			return 0;//TODO: set divide by zero error
		}
		if(!const_poly(r, 0)){
			return 0;
		}
		return scale_poly_modn(q, f, a, n);
	}
	if(f->len < g->len){//dividing by a polynomial with higher degree
		if(!const_poly(q, 0)){
			return 0;
		}
		return copy_poly(r, f);
	}
	
	//begin extended synthetic division
	//compute max length of quotient and remainder and extend their buffers if need be
	q->len = f->len - g->len + 1;
	r->len = g->len - 1;
	if(!ensure_poly_cap(q, q->len)){
		return 0;
	}else if(!ensure_poly_cap(r, r->len)){
		return 0;
	}
	
	//initialize column sums/coeffs of results
	memcpy(r->coeffs, f->coeffs, r->len*sizeof(int64_t));
	memcpy(q->coeffs, f->coeffs + r->len, q->len*sizeof(int64_t));
	
	//loop over quotient columns (coefficients in reverse order) which were initialized
	//to the first q->len dividend coefficients
	for(uint64_t k = q->len; k-- > 0;){
		//finish the column by dividing the sum by the leading coefficient of the divisor
		//for monic divisors (aka most of the time) a will simply be 1 so we may optimize this
		q->coeffs[k] = mod(q->coeffs[k]*a, n);
		//subtract the adjusted column sum times the some of the divisor coefficients from the
		//remainder coefficients.  the remainder coefficients were initialized to the last r->len
		//dividend coefficients.  we start q->len - k columns after the current column we just finished.
		//if k == 0 we should go from coefficient 0 in the remainder to coefficient r->len - 1
		//so in general we should go from k to r->len - 1 (this can easily be an empty interval)
		for(uint64_t i = 0, j = k; j < r->len; ++i, ++j){
			r->coeffs[j] = mod(r->coeffs[j] - q->coeffs[k]*g->coeffs[i], n);
		}
		//j goes from max(k - g->len + 1, 0) to k - 1 (both inclusive)
		//i ends at g->len - 2 (inclusive) and should go through the same number of values
		//so i starts at g->len - 2 - (k - 1 - j) = g->len + j - k - 1
		for(uint64_t j = k > g->len - 1 ? k - g->len + 1 : 0, i = g->len + j - k - 1; i < g->len - 1; ++i, ++j){
			q->coeffs[j] = mod(q->coeffs[j] - q->coeffs[k]*g->coeffs[i], n);
		}
	}
	normalize_poly(r);
	return 1;
}

void normalize_poly(poly_t *f){
	while(f->len > 1 && !f->coeffs[f->len - 1]){
		--f->len;
	}
}

void normalize_poly_modn(poly_t *f, int64_t n, int use_negatives){
	int64_t offset = use_negatives ? (1-n)/2 : 0;
	for(uint64_t i = 0; i < f->len; ++i){
		f->coeffs[i] = offset + mod(f->coeffs[i] - offset, n);
	}
	normalize_poly(f);
}

void normalize_exponents_modn(poly_t *f, uint64_t cn){
	for(uint64_t i = f->len - 1; i > cn; --i){
		uint64_t j = 1 + (i - 1)%cn;
		f->coeffs[j] += f->coeffs[i];
	}
	if(f->len > cn + 1){
		f->len = cn + 1;
	}
	normalize_poly(f);
}

int fprint_poly(FILE *file, const poly_t *f, const char *vname, const char *add, const char *sub, const char *pow, int descending){
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

static void skip_whitespace(const char **_str){
	while(isspace(**_str)){
		++*_str;
	}
}

static int parse_monomial(poly_t *f, const char **_str){
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
	if(!zero_extend_poly(f, x + 1)){
		return 0;
	}
	f->coeffs[x] += sign*coeff;
	*_str = str;
	return 1;
}

int str_to_poly(poly_t *f, int64_t *n, const char *str, const char **end){
	if(!f || !n){
		return 0;
	}
	if(!const_poly(f, 0) || !parse_monomial(f, &str)){
		return 0;
	}
	while(*str == '+' || *str == '-'){
		if(!parse_monomial(f, &str)){
			if(end){
				skip_whitespace(&str);
				*end = str;
			}
			normalize_poly(f);
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
			normalize_poly_modn(f, *n, 0);
			return 2;
		}else{
			str -= 4;
		}
	}
	if(end){
		skip_whitespace(&str);
		*end = str;
	}
	normalize_poly(f);
	return 1;
}

int rand_poly_modn(poly_t *f, uint64_t max_len, int64_t n){
	if(!max_len){
		return const_poly(f, 0);
	}
	if(!ensure_poly_cap(f, max_len)){
		return 0;
	}
	for(uint64_t i = 0; i < max_len; ++i){
		f->coeffs[i] = rand_u64(0, n);
	}
	f->len = max_len;
	normalize_poly(f);
	return 1;
}

int gcd_poly_modn(poly_t *d, const poly_t *restrict f, const poly_t *restrict g, int64_t n, poly_t tmps[static 3]){
	//tmps: qt, r0t, r1t
	poly_t *tmp, *r0, *r1, *r2;
	int status = 1;
	if(g->len > f->len){//Ensure the degree of f is >= the degree of g
		const poly_t *tmp = f;//Shadow the other tmp with a const version
		f = g;
		g = tmp;
	}
	if(g->len == 1){//If one of the inputs is a constant, we either have the gcd of something and a unit or something and zero
		status = g->coeffs[0] ? const_poly(d, 1) : copy_poly(d, f);
		goto CLEANUP;
	}
	
	r0 = d, r1 = tmps + 1, r2 = tmps + 2;
	//Unroll first two remainder calculations to prevent copying
	if(!quotrem_poly_modn(tmps + 0, r0, f, g, n)){
		status = 0;
		goto CLEANUP;
	}
	if(r0->len == 1 && !r0->coeffs[0]){
		status = copy_poly(d, g);
		goto CLEANUP;
	}
	
	if(!quotrem_poly_modn(tmps + 0, r1, g, r0, n)){
		status = 0;
		goto CLEANUP;
	}
	if(r1->len == 1 && !r1->coeffs[0]){
		goto CLEANUP;// r0 == d in this case, no need to copy
	}
	
	//Euclidean algorithm: take remainders until we reach 0, then the last nonzero remainder is the gcd
	while(1){
		if(!quotrem_poly_modn(tmps + 0, r2, r0, r1, n)){
			status = 0;
			goto CLEANUP;
		}
		if(r2->len == 1 && !r2->coeffs[0]){
			if(r1 != d){//Make sure the result is in the output and not a temporary
				status = copy_poly(d, r1);
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
			scale_poly_modn(d, d, a, n);
		}else if(a > 1){//If the leading coefficient is 0 or 1 the gcd is already monic
			int64_t c, g = egcd(a, n, &c, NULL);//This g shadows the input
			if(g != 1){//The leading coefficient cannot be inverted mod n
				return 0;
			}
			scale_poly_modn(d, d, c, n);
		}
	}
	return status;
}

int gcd_poly_modn_tmptmp(poly_t *d, const poly_t *restrict f, const poly_t *restrict g, int64_t n){
	poly_t tmps[3] = {};
	uint64_t min_len, quot_len;
	if(f->len <= g->len){
		min_len = f->len;
		quot_len = g->len - f->len + 1;
	}else{
		min_len = g->len;
		quot_len = f->len - g->len + 1;
	}
	int status = init_poly(tmps + 0, quot_len) &&
	         init_poly(tmps + 1, min_len) &&
	         init_poly(tmps + 2, min_len) &&
	         gcd_poly_modn(d, f, g, n, tmps);
	for(uint64_t i = 0; i < 3; ++i){
		if(tmps[i].cap){
			destroy_poly(tmps + i);
		}
	}
	return status;
}

int powmod_poly_modn(poly_t *h, const poly_t *restrict f, uint64_t e, const poly_t *restrict g, int64_t n, poly_t tmps[static 3]){
	//tmps: qt, st, rt
	if(g->len <= f->len){
		if(!quotrem_poly_modn(tmps + 0, tmps + 1, f, g, n)){
			return 0;
		}
	}else if(!copy_poly(tmps + 1, f)){
		return 0;
	}
	if(tmps[1].len == 1){
		if(!tmps[1].coeffs[0]){
			return const_poly(h, 0);
		}else if(!e){
			return const_poly(h, 1);
		}
	}
	if(!ensure_poly_cap(tmps + 1, 2*g->len - 3) || !ensure_poly_cap(tmps + 2, 2*g->len - 3) || !ensure_poly_cap(h, 2*g->len - 3)){
		return 0;
	}
	
	poly_t *t = h, *s = tmps + 1, *r = tmps + 2;
	while(e%2 == 0){
		mul_poly_modn(t, s, s, n);
		if(!quotrem_poly_modn(tmps + 0, s, t, g, n)){
			return 0;
		}//s = s*s%g
		e >>= 1;
	}
	copy_poly(r, s);
	while((e >>= 1)){
		mul_poly_modn(t, s, s, n);
		if(!quotrem_poly_modn(tmps + 0, s, t, g, n)){
			return 0;
		}//s = s*s%g
		if(e%2){
			mul_poly_modn(t, r, s, n);
			if(!quotrem_poly_modn(tmps + 0, r, t, g, n)){
				return 0;
			}//r = r*s%g
		}
	}
	
	if(r != h){// TODO: this check is always true, fix buffer juggling so we never need to copy
		if(!copy_poly(h, r)){
			return 0;
		}
	}
	return 1;
}

int powmod_poly_modn_tmptmp(poly_t *h, const poly_t *restrict f, uint64_t e, const poly_t *restrict g, int64_t n){
	poly_t tmps[3] = {};
	int status = 1;
	for(uint64_t i = 0; status && i < 3; ++i){
		status = init_poly(tmps + i, 2*g->len - 1);
	}
	if(status){
		status = powmod_poly_modn(h, f, e, g, n, tmps);
	}
	for(uint64_t i = 0; i < 3; ++i){
		if(tmps[i].cap){
			destroy_poly(tmps + i);
		}
	}
	return status;
}

int factors_d_poly_modn(poly_t *f_d, const poly_t *f, uint64_t d, int64_t n, poly_t tmps[static 4]){
	//tmps: xt, qt, st, rt
	if(!ensure_poly_cap(f_d, f->len)){
		return 0;
	}
	f_d->coeffs[0] = 0;
	f_d->coeffs[1] = 1;
	f_d->len = 2;
	if(!powmod_poly_modn(tmps + 0, f_d, pow_u64(n, d), f, n, tmps + 1)){
		return 0;
	}
	if(!sub_poly_modn(tmps + 0, tmps + 0, f_d, n)){
		return 0;
	}
	return gcd_poly_modn(f_d, tmps + 0, f, n, tmps + 1);
}

int factor1_poly_modn(poly_t *g, const poly_t *f, uint64_t d, int64_t n, poly_t tmps[static 4]){
	//tmps: xt, qt, st, rt
	while(1){
		if(!rand_poly_modn(tmps + 0, f->len - 1, n)){
			return 0;
		}
		if(!gcd_poly_modn(g, tmps + 0, f, n, tmps + 1)){
			return 0;
		}
		if(g->len > 1 && g->len < f->len){
			return 1;
		}
		if(!powmod_poly_modn(g, tmps + 0, (pow_u64(n, d)-1)/2, f, n, tmps + 1)){
			return 0;
		}
		const_poly(tmps + 0, 1);
		add_poly_modn(tmps + 0, g, tmps + 0, n);
		if(!gcd_poly_modn(g, tmps + 0, f, n, tmps + 1)){
			return 0;
		}
		if(g->len > 1 && g->len < f->len){
			return 1;
		}
	}
}

static inline int roots_polyn_modn_rec(poly_roots_t *roots, int64_t n, poly_t tmps[static 6]){
	//tmps: gt, ft, xt, qt, st, rt
	while(1){
		//fprintf(stderr, "\e[1;33mFactoring (");
		//fprint_poly(stderr, tmps + 1, "x", " + ", " - ", "**", 0);
		//fprintf(stderr, ") mod %"PRId64"\e[0m\n", n);
		if(tmps[1].len == 2){//linear factor
			//fprintf(stderr, "\e[1;33m polynomial is linear\e[0m\n");
			if(tmps[1].coeffs[0]){//nonzero root
				if(tmps[1].coeffs[1] == 1){//monic
					roots->roots[roots->len++] = n - tmps[1].coeffs[0];
				}else{//not monic
					int64_t a;
					egcd(tmps[1].coeffs[0], n, &a, NULL);
					roots->roots[roots->len++] = mod(-tmps[1].coeffs[0]*a, n);
				}
			}else{//zero root
				roots->roots[roots->len++] = 0;
			}
			return 1;
		}else if(tmps[1].len == 3){//quadratic factor
			//fprintf(stderr, "\e[1;33m polynomial is quadratic\e[0m\n");
			if(tmps[1].coeffs[2] != 1){//make monic
				int64_t a;
				egcd(tmps[1].coeffs[2], n, &a, NULL);
				scale_poly_modn(tmps + 1, tmps + 1, a, n);
			}
			int64_t c = tmps[1].coeffs[0];
			int64_t b = tmps[1].coeffs[1];
			int64_t r = sqrt_mod(mod(b*b - 4*c, n), n);
			roots->roots[roots->len++] = mod((n+1)/2*(-b + r), n);
			roots->roots[roots->len++] = mod((n+1)/2*(-b - r), n);
			return 1;
		}
		if(!factor1_poly_modn(tmps + 0, tmps + 1, 1, n, tmps + 2)){
			return 0;
		}
		if(!quotrem_poly_modn(tmps + 3, tmps + 5, tmps + 1, tmps + 0, n)){
			return 0;
		}
		//tmps[0] and tmps[3] hold the nontrivial factors we found
		//we have to do a recursive call on the factor with smaller degree
		//if it is linear or quadratic, we do not need to copy the larger factor
		//otherwise we do
		poly_t tmp;
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
		//fprint_poly(stderr, tmps + 1, "x", " + ", " - ", "**", 0);
		//fprintf(stderr, ")(");
		//fprint_poly(stderr, tmps + 3, "x", " + ", " - ", "**", 0);
		//fprintf(stderr, ")\e[0m\n");
		int have_small_factor = tmps[1].len <= 3;
		if(!have_small_factor){
			if(!init_poly(&tmp, tmps[3].len)){
				return 0;
			}
			copy_poly(&tmp, tmps + 3);
		}
		if(!roots_polyn_modn_rec(roots, n, tmps)){
			return 0;
		}
		if(!have_small_factor){
			copy_poly(tmps + 3, &tmp);
			destroy_poly(&tmp);
		}
		tmp = tmps[3];
		tmps[3] = tmps[1];
		tmps[1] = tmp;
	}
}

int roots_poly_modn(const poly_t *f, int64_t n, poly_roots_t *roots, poly_t tmps[static 6]){
	//tmps: gt, ft, xt, qt, st, rt
	if(!factors_d_poly_modn(tmps + 1, f, 1, n, tmps + 2)){
		return 0;
	}
	roots->len = 0;
	if(tmps[1].len == 1){
		return 1;
	}
	if(!ensure_poly_roots_cap(roots, tmps[1].len - 1)){
		return 0;
	}
	return roots_polyn_modn_rec(roots, n, tmps);
}

int roots_poly_modn_tmptmp(const poly_t *f, int64_t n, poly_roots_t *roots){
	poly_t tmps[6] = {};
	int status = 1;
	for(uint64_t i = 0; status && i < 6; ++i){
		status = init_poly(tmps + i, f->len);
	}
	if(status){
		status = roots_poly_modn(f, n, roots, tmps);
	}
	for(uint64_t i = 0; i < 6; ++i){
		if(tmps[i].cap){
			destroy_poly(tmps + i);
		}
	}
	return status;
}

