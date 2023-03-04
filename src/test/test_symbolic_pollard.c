#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nut/polynomial.h>
#include <nut/factorization.h>

int main(){
	int64_t n;
	poly_t tmps[6] = {{}, {}, {}, {}, {}, {}};
	poly_t *x = tmps + 2, *y = tmps + 3, *g = tmps + 4, *h = tmps + 5;
	const char *end = "";
	if(str_to_poly(g, &n, "x^2 + 1 mod 9", &end) != 2){
		fprintf(stderr, "\e[1;31mFailed to parse polynomial\e[0m\n");
		exit(0);
	}
	if(n < 2){
		fprintf(stderr, "\e[1;31mn too small!\e[0m\n");
		exit(0);
	}
	factors_t *factors = init_factors_t_ub(n, 25, small_primes);
	if(!factors){
		fprintf(stderr, "\e[1;31mCouldn't allocate factorization struct\e[0m\n");
		exit(0);
	}
	{
		uint64_t res = factor_heuristic(n, 25, small_primes, &default_factor_conf, factors);
		if(res != 1){
			fprintf(stderr, "\e[1;31mFactoring %"PRIu64" failed\e[0m\n", n);
			exit(0);
		}
	}
	uint64_t cn = carmichael_lambda(factors);
	free(factors);
	if(!ensure_poly_cap(x, cn + 1) || !ensure_poly_cap(y, cn + 1)){
		fprintf(stderr, "\e[1;31mCouldn't allocate polynomials\e[0m");
		exit(0);
	}
	x->coeffs[0] = y->coeffs[0] = 0;
	x->coeffs[1] = y->coeffs[1] = 1;
	x->len = y->len = 2;
	for(uint64_t i = 0; i < (uint64_t)n; ++i){
		printf("%"PRIu64") x: ", i);
		fprint_poly(stdout, x, "x", " + ", " - ", "^", 1);
		printf("; y: ");
		fprint_poly(stdout, y, "x", " + ", " - ", "^", 1);
		putc('\n', stdout);
		compose_poly_modn(h, g, x, n, cn, tmps);
		normalize_poly_modn(h, n, 0);
		void *tmp = h;
		h = x;
		x = tmp;
		compose_poly_modn(h, g, y, n, cn, tmps);
		normalize_poly_modn(h, n, 0);
		compose_poly_modn(y, g, h, n, cn, tmps);
		normalize_poly_modn(y, n, 0);
	}
	for(uint64_t i = 0; i < 6; ++i){
		destroy_poly(tmps + i);
	}
}

