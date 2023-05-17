#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nut/polynomial.h>
#include <nut/factorization.h>

int main(){
	int64_t n;
	nut_Poly tmps[6] = {{}, {}, {}, {}, {}, {}};
	nut_Poly *x = tmps + 2, *y = tmps + 3, *g = tmps + 4, *h = tmps + 5;
	const char *end = "";
	if(nut_Poly_parse(g, &n, "x^2 + 1 mod 9", &end) != 2){
		fprintf(stderr, "\e[1;31mFailed to parse polynomial\e[0m\n");
		exit(0);
	}
	if(n < 2){
		fprintf(stderr, "\e[1;31mn too small!\e[0m\n");
		exit(0);
	}
	nut_Factors *factors = nut_make_Factors_ub(n, 25, nut_small_primes);
	if(!factors){
		fprintf(stderr, "\e[1;31mCouldn't allocate factorization struct\e[0m\n");
		exit(0);
	}
	{
		uint64_t res = nut_u64_factor_heuristic(n, 25, nut_small_primes, &nut_default_factor_conf, factors);
		if(res != 1){
			fprintf(stderr, "\e[1;31mFactoring %"PRIu64" failed\e[0m\n", n);
			exit(0);
		}
	}
	uint64_t cn = nut_Factor_carmichael(factors);
	free(factors);
	if(!nut_Poly_ensure_cap(x, cn + 1) || !nut_Poly_ensure_cap(y, cn + 1)){
		fprintf(stderr, "\e[1;31mCouldn't allocate polynomials\e[0m");
		exit(0);
	}
	x->coeffs[0] = y->coeffs[0] = 0;
	x->coeffs[1] = y->coeffs[1] = 1;
	x->len = y->len = 2;
	for(uint64_t i = 0; i < (uint64_t)n; ++i){
		printf("%"PRIu64") x: ", i);
		nut_Poly_fprint(stdout, x, "x", " + ", " - ", "^", 1);
		printf("; y: ");
		nut_Poly_fprint(stdout, y, "x", " + ", " - ", "^", 1);
		putc('\n', stdout);
		nut_Poly_compose_modn(h, g, x, n, cn, tmps);
		nut_Poly_normalize_modn(h, n, 0);
		void *tmp = h;
		h = x;
		x = tmp;
		nut_Poly_compose_modn(h, g, y, n, cn, tmps);
		nut_Poly_normalize_modn(h, n, 0);
		nut_Poly_compose_modn(y, g, h, n, cn, tmps);
		nut_Poly_normalize_modn(y, n, 0);
	}
	for(uint64_t i = 0; i < 6; ++i){
		nut_Poly_destroy(tmps + i);
	}
}

