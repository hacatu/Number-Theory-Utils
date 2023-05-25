#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nut/polynomial.h>
#include <nut/factorization.h>
#include <nut/debug.h>

int main(){
	uint64_t l = 10;
	int64_t p = 1048573;
	uint64_t trials = 1000;
	uint64_t approx_eq_pts = 5;
	nut_Poly f[1], g[1], h[1], q[1], r[1];
	if(!nut_Poly_init(f, l) || !nut_Poly_init(g, l) || !nut_Poly_init(h, l) || !nut_Poly_init(q, l) || !nut_Poly_init(r, l)){
		fprintf(stderr, "\e[1;31mERROR: Could not allocate polynomials.\e[0m\n");
		exit(EXIT_FAILURE);
	}
	uint64_t passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		"addition",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, l, p);
		nut_Poly_rand_modn(g, l, p);
		nut_Poly_add_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if((nut_Poly_eval_modn(f, x, p) + nut_Poly_eval_modn(g, x, p))%p != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		if(approx_eq){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		"subtraction",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, l, p);
		nut_Poly_rand_modn(g, l, p);
		nut_Poly_sub_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if(nut_i64_mod(nut_Poly_eval_modn(f, x, p) - nut_Poly_eval_modn(g, x, p), p) != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		if(approx_eq){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" random polynomials...\e[0m\n",
		"scalar product",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, l, p);
		int64_t a = nut_u64_prand(0, p);
		nut_Poly_scale_modn(g, f, a, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if(nut_Poly_eval_modn(f, x, p)*a%p != nut_Poly_eval_modn(g, x, p)){
				approx_eq = 0;
				break;
			}
		}
		if(approx_eq){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		"multiplication",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, l, p);
		nut_Poly_rand_modn(g, l, p);
		nut_Poly_mul_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if(nut_Poly_eval_modn(f, x, p)*nut_Poly_eval_modn(g, x, p)%p != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		if(approx_eq){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting division of degree <=8 polynomials by degree <=4 polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, 8, p);
		nut_Poly_rand_modn(g, 4, p);
		nut_Poly_quotrem_modn(q, r, f, g, p);
		//should have q*g + r = f mod p
		nut_Poly_mul_modn(h, q, g, p);
		nut_Poly_add_modn(h, h, r, p);
		if(f->len == h->len && !memcmp(h->coeffs, f->coeffs, f->len*sizeof(int64_t))){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	trials = 3;
	passed = 0;
	fprintf(stderr, "\e[1;34mTesting \"nut_Poly_setconst\"/\"nut_Poly_eval_modn\" for 0, 1, and 2 mod 3...\e[0m\n");
	for(int64_t i = 0; i < (int64_t)trials; ++i){
		if(nut_Poly_setconst(f, i) && nut_Poly_eval_modn(f, i, 3) == i){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	trials = 1000;
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=4+8 polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		"commutativity of addition",
		p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, 4, p);
		nut_Poly_rand_modn(g, 8, p);
		nut_Poly_add_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if((nut_Poly_eval_modn(f, x, p) + nut_Poly_eval_modn(g, x, p))%p != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		nut_Poly_add_modn(g, g, f, p);
		if(approx_eq && g->len == h->len && !memcmp(g->coeffs, h->coeffs, g->len*sizeof(int64_t))){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=4-8 polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		"reversing subtraction",
		p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, 4, p);
		nut_Poly_rand_modn(g, 8, p);
		nut_Poly_sub_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if(nut_i64_mod(nut_Poly_eval_modn(f, x, p) - nut_Poly_eval_modn(g, x, p), p) != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		nut_Poly_sub_modn(g, g, f, p);
		nut_Poly_scale_modn(g, g, p-1, p);
		if(approx_eq && g->len == h->len && !memcmp(g->coeffs, h->coeffs, g->len*sizeof(int64_t))){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" random polynomials...\e[0m\n",
		"simple multiplication",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, 10, p);
		nut_Poly_rand_modn(g, 1, p);
		nut_Poly_mul_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if(nut_Poly_eval_modn(f, x, p)*nut_Poly_eval_modn(g, x, p)%p != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		nut_Poly_mul_modn(h, g, f, p);
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = nut_u64_prand(0, p);
			if(nut_Poly_eval_modn(f, x, p)*nut_Poly_eval_modn(g, x, p)%p != nut_Poly_eval_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		if(approx_eq){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" random polynomials...\e[0m\n",
		"scalar division",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, 10, p);
		do{
			nut_Poly_rand_modn(g, 1, p);
		}while(!g->coeffs[0]);
		nut_Poly_quotrem_modn(q, r, f, g, p);
		nut_Poly_scale_modn(q, q, g->coeffs[0], p);
		if(r->len == 1 && !r->coeffs[0] && q->len == f->len && !memcmp(q->coeffs, f->coeffs, f->len*sizeof(int64_t))){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" random polynomials...\e[0m\n",
		"simple division",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		nut_Poly_rand_modn(f, 10, p);
		nut_Poly_rand_modn(g, 1, p);
		nut_Poly_quotrem_modn(q, r, g, f, p);
		if(q->len == 1 && !q->coeffs[0] && r->len == g->len && !memcmp(r->coeffs, g->coeffs, g->len*sizeof(int64_t))){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	g->coeffs[0] = 1;
	g->coeffs[1] = 7;
	g->coeffs[2] = 5;
	g->coeffs[3] = 21;
	g->coeffs[4] = 1;
	g->len = 5;
	
	check_alloc("poly f", (void*)nut_Poly_ensure_cap(f, 24));
	memset(f->coeffs, 0, 23*sizeof(int64_t));
	f->coeffs[23] = 1;
	f->coeffs[1] = 22;
	f->len = 24;
	
	fprintf(stderr, "\e[1;34mComputing gcd of (");
	nut_Poly_fprint(stderr, f, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") and (");
	nut_Poly_fprint(stderr, g, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") mod 23...\e[0m\n");
	
	if(!nut_Poly_gcd_modn_tmptmp(h, f, g, 23)){
		fprintf(stderr, "\e[1;31mERROR: \"nut_Poly_gcd_modn_tmptmp\" failed.\e[0m\n");
	}else{
		fprintf(stderr, "\e[1;33m got (");
		nut_Poly_fprint(stderr, h, "x", " + ", " - ", "**", 1);
		fprintf(stderr, ")\e[0m\n");
	}
	
	
	
	f->coeffs[1] = 1;
	f->len = 2;
	
	fprintf(stderr, "\e[1;34mComputing gcd of ((");
	nut_Poly_fprint(stderr, f, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ")**23 - x) and (");
	nut_Poly_fprint(stderr, g, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") mod 23, with an initial modular exponentation step...\e[0m\n");
	
	if(!nut_Poly_powmod_modn_tmptmp(h, f, 23, g, 23)){
		fprintf(stderr, "\e[1;31mERROR: \"nut_Poly_powmod_modn_tmptmp\" failed.\e[0m\n");
	}else{
		nut_Poly_sub_modn(f, h, f, 23);
		if(!nut_Poly_gcd_modn_tmptmp(h, f, g, 23)){
			fprintf(stderr, "\e[1;31mERROR: \"nut_Poly_gcd_modn_tmptmp\" failed.\e[0m\n");
		}else{
			fprintf(stderr, "\e[1;33m got (");
			nut_Poly_fprint(stderr, h, "x", " + ", " - ", "**", 1);
			fprintf(stderr, ")\e[0m\n");
		}
	}
	
	
	
	nut_Poly_destroy(f);
	nut_Poly_destroy(g);
	nut_Poly_destroy(h);
	nut_Poly_destroy(q);
	nut_Poly_destroy(r);
}

