#include <stdio.h>
#include <nut/dirichlet.h>
#include <nut/modular_math.h>
#include <nut/debug.h>
#include <nut/sieves.h>

static const uint64_t N = 1'000'000;

int main(){
	nut_Diri pi_tbl [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri_init(&pi_tbl, N, 0);
	nut_Diri_compute_pi(&pi_tbl);
	uint8_t *is_composite [[gnu::cleanup(cleanup_free)]] = nut_sieve_is_composite(N);
	uint64_t *pi_packed [[gnu::cleanup(cleanup_free)]] = nut_compute_pi_range(N, is_composite);

	for(int64_t i = 1; i <= pi_tbl.y; ++i){
		if((int64_t)nut_compute_pi_from_tables(i, pi_packed, is_composite) != nut_Diri_get_dense(&pi_tbl, i)){
			fprintf(stderr, "\e[1;31j163 table was wrong at dense %"PRIi64"\e[0m\n", i);
		}
	}
	for(int64_t i = 1; i < pi_tbl.yinv; ++i){
		if((int64_t)nut_compute_pi_from_tables(pi_tbl.x/i, pi_packed, is_composite) != nut_Diri_get_sparse(&pi_tbl, i)){
			fprintf(stderr, "\e[1;31j163 table was wrong at sparse %"PRIi64"\e[0m\n", i);
		}
	}
}

