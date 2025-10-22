#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

#include <nut/factorization.h>
#include <nut/sieves.h>
#include <nut/modular_math.h>
#include <nut/debug.h>

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

static int dvisitor_find_order(const nut_Factors *dfactors, uint64_t d, void *_data){
	uint64_t *data = _data;
	if(d < data[1] && nut_u64_powmod(10, d, data[0]) == 1){
		data[1] = d;
	}
	return 0;
}

int main(){
	uint64_t max = 10000, passed = 0, failed = 0;
	nut_Factors *factors [[gnu::cleanup(cleanup_free)]] = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	uint64_t *largest_factors [[gnu::cleanup(cleanup_free)]] = nut_sieve_largest_factors(max);
	uint64_t *carmichael_vals [[gnu::cleanup(cleanup_free)]] = nut_sieve_carmichael(max);
	fprintf(stderr, "\e[1;34mChecking order of 10 mod all coprime numbers from 2 to %"PRIu64"...\e[0m\n", max);
	for(uint64_t n = 2; n <= max; ++n){
		if(nut_i64_egcd(n, 10, NULL, NULL) != 1){
			continue;
		}
		uint64_t cn = carmichael_vals[n];
		nut_fill_factors_from_largest(factors, cn, largest_factors);
		uint64_t data[2] = {n, cn};
		nut_Factor_forall_divs_tmptmp(factors, dvisitor_find_order, data);
		if(nut_u64_order_mod(10, n, cn, factors) != data[1]){
			fprintf(stderr, "\e[1;31morder_mod does not match the value from forall_divs reduction\e[0m\n");
			++failed;
		}else{
			++passed;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", !failed ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, passed + failed);
}

