#include <stdio.h>
#include <inttypes.h>

#include <nut/debug.h>
#include <nut/modular_math.h>
#include <nut/factorization.h>

NUT_ATTR_NO_SAN("unsigned-integer-overflow")
static void test_modinv(uint64_t trials){
	uint64_t passed = 0;
	fprintf(stderr, "\e[1;34mInverting %"PRIu64" random odd numbers mod powers of 2...\e[0m\n", trials);
	for(uint64_t i = 0, t; i < trials; ++i){
		uint64_t x = (nut_u64_rand(0, 1ull << 63) << 1) + 1;
		for(t = 1; t <= 64; ++t){
			uint64_t xinv = nut_u64_modinv_2t(x, t);
			bool failed;
			if(t < 64){
				failed = (xinv & ~((1ull << t) - 1)) || // xinv is not reduced mod 2**t
				(((xinv*x) & ((1ull << t) - 1)) != 1); // xinv is not actually an inverse of x mod 2**t
			}else{
				failed = xinv*x != 1;
			}
			if(failed){
				fprintf(stderr, "\e[1;31mmodinv_2t failed to invert %"PRIu64" mod 2**%"PRIu64"\e[0m\n", x, t);
				break;
			}
		}
		if(t == 65){
			++passed;
		}
	}
	print_summary("modular inverse", passed, trials);
}

static void test_binom_simple(uint64_t trials, uint64_t binom_tbl[static 63]){
	uint64_t passed = 0;
	fprintf(stderr, "\e[1;34mTesting binomial coeffs until (61 choose k), the last n when the simple algorithm doesn't overflow...\e[0m\n");
	for(uint64_t n = 1; n <= 61; ++n){
		for(uint64_t k = n - 1; k; --k){
			binom_tbl[k] += binom_tbl[k - 1];
		}
		binom_tbl[n] = 1;
		bool success = true;
		uint64_t prev = 1;
		for(uint64_t k = 0; k <= n && success; ++k){
			if(nut_u64_binom(n, k) != binom_tbl[k]){
				fprintf(stderr, "\e[1;31mbinom(%"PRIu64", %"PRIu64") failed\e[0m\n", n, k);
				success = false;
			}
			if(k && (prev = nut_u64_binom_next(n, k, prev)) != binom_tbl[k]){
				fprintf(stderr, "\e[1;31mbinom_next(%"PRIu64", %"PRIu64") failed\e[0m\n", n, k);
				success = false;
			}
		}
		if(success){
			++passed;
		}
	}
	print_summary("binomial coefficients", passed, trials);
}

NUT_ATTR_NO_SAN("unsigned-integer-overflow")
static void test_binom_mod2t(uint64_t trials, uint64_t binom_tbl[static 63]){
	uint64_t passed = 0;
	fprintf(stderr, "\e[1;34mTesting (n choose k) mod 2**t for n = 62, ..., 561; k = 0, ..., 61; and t = 2, ..., 61...\e[0m\n");
	for(uint64_t n = 62; n <= 561; ++n){
		for(uint64_t k = 61; k; --k){
			binom_tbl[k] += binom_tbl[k - 1]; // this will overflow, but this ok because we build with -fwrapv and are computing mod powers of 2
		}
		bool success = true;
		for(uint64_t t = 2; t <= 64; ++t){
			uint64_t v2 = 0, p2 = 1;
			for(uint64_t k = 1; k <= 61 && success; ++k){
				uint64_t c = nut_u64_binom_next_mod_2t(n, k, t, &v2, &p2);
				if(c != (binom_tbl[k] & (t < 64 ? (1ull << t) - 1 : ~0ull))){
					fprintf(stderr, "\e[1;31mbinom_next_mod_2t(%"PRIu64", %"PRIu64") mod 2**%"PRIu64" failed\e[0m\n", n, k, t);
					success = false;
				}
			}
		}
		if(success){
			++passed;
		}else{
			break;
		}
	}
	print_summary("binomial coefficients", passed, trials);
}

int main(){
	test_modinv(1000);
	uint64_t binom_tbl[63] = {[0] = 1};
	test_binom_simple(61, binom_tbl);
	test_binom_mod2t(500, binom_tbl);
}

