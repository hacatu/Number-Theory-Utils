#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include <nut/modular_math.h>
#include <nut/polynomial.h>

const int64_t M = 20092010ull;
const int64_t K = 1000000000000000000ull;
const int64_t N = 2000ull;

int main(){
	nut_Poly f[1], a[1], x[1];
	if(!nut_Poly_init(f, N + 1) || !nut_Poly_init(a, N) || !nut_Poly_init(x, 2)){
		fprintf(stderr, "\e[1;31mERROR: Could not allocate polynomials.\e[0m\n");
		exit(EXIT_FAILURE);
	}
	memset(f->coeffs + 2, 0, (N - 1 - 2 + 1)*sizeof(int64_t));
	f->coeffs[0] = M - 1;
	f->coeffs[1] = M - 1;
	f->coeffs[N] = 1;
	f->len = N + 1;
	
	x->coeffs[0] = 0;
	x->coeffs[1] = 1;
	x->len = 2;
	
	nut_Poly_powmod_modn_tmptmp(a, x, K - N + 1, f, M);
	
	uint64_t s = a->coeffs[0];
	for(uint64_t j = 1; j < (uint64_t)N; ++j){
		s = nut_i64_mod(s + a->coeffs[j]*2, M);
	}
	
	printf("%"PRIu64"\n", s);
	
	nut_Poly_destroy(f);
	nut_Poly_destroy(x);
	nut_Poly_destroy(a);
}

