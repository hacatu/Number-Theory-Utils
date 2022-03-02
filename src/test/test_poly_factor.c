#include <string.h>

#include <nut/polynomial.h>
#include <nut/factorization.h>

int cmp_i64(const void *_a, const void *_b){
	int64_t a = *(const int64_t*)_a, b = *(const int64_t*)_b;
	if(a < b){
		return -1;
	}else if(a > b){
		return 1;
	}
	return 0;
}

int main(){
	poly_t f[1];
	if(!init_poly(f, 9)){
		fprintf(stderr, "\e[1;31mERROR: Could not allocate polynomials.\e[0m\n");
		exit(EXIT_FAILURE);
	}
	
	
	
	/*
	f->coeffs[0] = 1;
	f->coeffs[1] = -1;
	f->coeffs[2] = 1;
	f->coeffs[3] = -1;
	f->coeffs[4] = 1;
	f->len = 5;
	*/
	
	f->coeffs[0] = 1;
	f->coeffs[1] = 1;
	f->coeffs[2] = 0;
	f->coeffs[3] = -1;
	f->coeffs[4] = -1;
	f->coeffs[5] = -1;
	f->coeffs[6] = 0;
	f->coeffs[7] = 1;
	f->coeffs[8] = 1;
	f->len = 9;
	
	normalize_poly_modn(f, 5, 0);
	normalize_poly_modn(f, 5, 1);
	uint64_t trials = 5000;
	fprintf(stderr, "\e[1;34mComputing roots of (");
	fprint_poly(stderr, f, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") mod p for %"PRIu64" random primes...\e[0m\n", trials);
	
	poly_roots_t roots[1];
	if(!init_poly_roots(roots, 8)){
		fprintf(stderr, "\e[1;31mERROR: Could not allocate polynomial root buffer.\e[0m\n");
		exit(EXIT_FAILURE);
	}
	// Remember if O(12*n**2) >= O(p/lb(p)) linear search is faster than cantor zassenhaus so
	// for a degree 4 polynomial linear search is fast for 192 >= p/lb(p) --> p <= 2113
	// and for a degree 8 polynomial it is fast for 768 >= p/lb(p) --> p <= 10223
	uint64_t passed = 0;
	for(uint64_t i = 0; i < trials; ++i){
		int64_t p = rand_u64(2, 10224);
		while(!is_prime_dmr(p)){
			++p;
		}
		//fprintf(stderr, "\e[1;34mp=%"PRId64"\e[0m\n", p);
		int64_t brute_roots[8];
		uint64_t brute_roots_len = 0;
		for(int64_t r = 0; r < p; ++r){
			int64_t y = eval_poly_modn(f, r, p);
			if(!y){
				if(brute_roots_len >= 8){
					abort();
				}
				brute_roots[brute_roots_len++] = r;
			}
		}
		if(!roots_poly_modn_tmptmp(f, p, roots)){
			fprintf(stderr, "\e[1;31mERROR: \"roots_poly_modn\" failed for p=%"PRId64".\e[0m\n", p);
			continue;
		}
		if(roots->len != brute_roots_len){
			fprintf(stderr, "\e[1;31mERROR: \"roots_poly_modn\" gave wrong number of roots for p=%"PRId64".\e[0m\n", p);
			continue;
		}
		qsort(roots->roots, roots->len, sizeof(int64_t), cmp_i64);
		if(memcmp(brute_roots, roots->roots, brute_roots_len*sizeof(int64_t))){
			fprintf(stderr, "\e[1;31mERROR: \"roots_poly_modn\" gave wrong roots for p=%"PRId64".\e[0m\n", p);
			continue;
		}
		++passed;
		//fprintf(stderr, "\e[1;32mSUCCESS\e[0m\n");
	}
	destroy_poly_roots(roots);
	destroy_poly(f);
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
}

