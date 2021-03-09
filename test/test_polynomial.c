#include "../polynomial.h"
#include "../factorization.h"

#include <string.h>

int main(){
	uint64_t l = 10;
	uint64_t p = 1048573;
	uint64_t trials = 1000;
	uint64_t approx_eq_pts = 5;
	poly_t f[1], g[1], h[1], q[1], r[1];
	if(!init_poly(f, l) || !init_poly(g, l) || !init_poly(h, l) || !init_poly(q, l) || !init_poly(r, l)){
		fprintf(stderr, "\e[1;31mERROR: Could not allocate polynomials.\e[0m\n");
		exit(EXIT_FAILURE);
	}
	uint64_t
	
	
	
	passed = 0;
	fprintf(stderr,
		"\e[1;34mTesting %s of degree <=%"PRIu64" polynomials mod %"PRIu64" for %"PRIu64" pairs of random polynomials...\e[0m\n",
		"addition",
		l, p, trials
	);
	for(uint64_t i = 0; i < trials; ++i){
		rand_poly_modn(f, l, p);
		rand_poly_modn(g, l, p);
		add_poly_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if((eval_poly_modn(f, x, p) + eval_poly_modn(g, x, p))%p != eval_poly_modn(h, x, p)){
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
		rand_poly_modn(f, l, p);
		rand_poly_modn(g, l, p);
		sub_poly_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if(mod(eval_poly_modn(f, x, p) - eval_poly_modn(g, x, p), p) != eval_poly_modn(h, x, p)){
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
		rand_poly_modn(f, l, p);
		int64_t a = prand_u64(0, p);
		scale_poly_modn(g, f, a, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if(eval_poly_modn(f, x, p)*a%p != eval_poly_modn(g, x, p)){
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
		rand_poly_modn(f, l, p);
		rand_poly_modn(g, l, p);
		mul_poly_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if(eval_poly_modn(f, x, p)*eval_poly_modn(g, x, p)%p != eval_poly_modn(h, x, p)){
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
		rand_poly_modn(f, 8, p);
		rand_poly_modn(g, 4, p);
		quotrem_poly_modn(q, r, f, g, p);
		//should have q*g + r = f mod p
		mul_poly_modn(h, q, g, p);
		add_poly_modn(h, h, r, p);
		if(f->len == h->len && !memcmp(h->coeffs, f->coeffs, f->len*sizeof(int64_t))){
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
	
	
	
	trials = 3;
	passed = 0;
	fprintf(stderr, "\e[1;34mTesting \"const_poly\"/\"eval_poly_modn\" for 0, 1, and 2 mod 3...\e[0m\n");
	for(uint64_t i = 0; i < trials; ++i){
		const_poly(f, i);
		if(eval_poly_modn(f, i, 3) == i){
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
		rand_poly_modn(f, 4, p);
		rand_poly_modn(g, 8, p);
		add_poly_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if((eval_poly_modn(f, x, p) + eval_poly_modn(g, x, p))%p != eval_poly_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		add_poly_modn(g, g, f, p);
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
		rand_poly_modn(f, 4, p);
		rand_poly_modn(g, 8, p);
		sub_poly_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if(mod(eval_poly_modn(f, x, p) - eval_poly_modn(g, x, p), p) != eval_poly_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		sub_poly_modn(g, g, f, p);
		scale_poly_modn(g, g, p-1, p);
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
		rand_poly_modn(f, 10, p);
		rand_poly_modn(g, 1, p);
		mul_poly_modn(h, f, g, p);
		int approx_eq = 1;
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if(eval_poly_modn(f, x, p)*eval_poly_modn(g, x, p)%p != eval_poly_modn(h, x, p)){
				approx_eq = 0;
				break;
			}
		}
		mul_poly_modn(h, g, f, p);
		for(uint64_t j = 0; j < approx_eq_pts; ++j){
			uint64_t x = prand_u64(0, p);
			if(eval_poly_modn(f, x, p)*eval_poly_modn(g, x, p)%p != eval_poly_modn(h, x, p)){
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
		rand_poly_modn(f, 10, p);
		do{
			rand_poly_modn(g, 1, p);
		}while(!g->coeffs[0]);
		quotrem_poly_modn(q, r, f, g, p);
		scale_poly_modn(q, q, g->coeffs[0], p);
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
		rand_poly_modn(f, 10, p);
		rand_poly_modn(g, 1, p);
		quotrem_poly_modn(q, r, g, f, p);
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
	
	ensure_poly_cap(f, 24);
	memset(f->coeffs, 0, 23*sizeof(int64_t));
	f->coeffs[23] = 1;
	f->coeffs[1] = 22;
	f->len = 24;
	
	fprintf(stderr, "\e[1;34mComputing gcd of (");
	fprint_poly(stderr, f, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") and (");
	fprint_poly(stderr, g, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") mod 23...\e[0m\n");
	
	if(!gcd_poly_modn_tmptmp(h, f, g, 23)){
		fprintf(stderr, "\e[1;31mERROR: \"gcd_poly_modn_tmptmp\" failed.\e[0m\n");
	}else{
		fprintf(stderr, "\e[1;33m got (");
		fprint_poly(stderr, h, "x", " + ", " - ", "**", 1);
		fprintf(stderr, ")\e[0m\n");
	}
	
	
	
	f->coeffs[1] = 1;
	f->len = 2;
	
	fprintf(stderr, "\e[1;34mComputing gcd of ((");
	fprint_poly(stderr, f, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ")**23 - x) and (");
	fprint_poly(stderr, g, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") mod 23, with an initial modular exponentation step...\e[0m\n");
	
	if(!powmod_poly_modn_tmptmp(h, f, 23, g, 23)){
		fprintf(stderr, "\e[1;31mERROR: \"powmod_poly_modn_tmptmp\" failed.\e[0m\n");
	}else{
		sub_poly_modn(f, h, f, 23);
		if(!gcd_poly_modn_tmptmp(h, f, g, 23)){
			fprintf(stderr, "\e[1;31mERROR: \"gcd_poly_modn_tmptmp\" failed.\e[0m\n");
		}else{
			fprintf(stderr, "\e[1;33m got (");
			fprint_poly(stderr, h, "x", " + ", " - ", "**", 1);
			fprintf(stderr, ")\e[0m\n");
		}
	}
	
	
	
	destroy_poly(f);
	destroy_poly(g);
	destroy_poly(h);
	destroy_poly(q);
	destroy_poly(r);
}

