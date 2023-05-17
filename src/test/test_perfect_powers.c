#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

#include <nut/factorization.h>

static const char *get_ord_suffix(uint64_t x){
	x %= 100;
	if(4 <= x && x <= 20){
		return "th";
	}
	switch(x%10){
		case 1: return "st";
		case 2: return "nd";
		case 3: return "rd";
	}
	return "th";
}

static bool check_power(uint64_t b, uint64_t x){
	uint64_t a = nut_u64_pow(b, x);
	uint64_t base_a, exp_a;
	uint64_t base_b, exp_b;
	if(!nut_u64_is_perfect_power(a, 63, &base_a, &exp_a)){
		fprintf(stderr, "\e[1;31mis_perfect_power failed to recognize %"PRIu64"**%"PRIu64" = %"PRIu64".\e[0m\n", b, x, a);
		return false;
	}
	if(base_a != b){
		if(!nut_u64_is_perfect_power(b, 31, &base_b, &exp_b)){
			fprintf(stderr, "\e[1;31mis_perfect_power got wrong base for %"PRIu64"**%"PRIu64" = %"PRIu64" or misidentified the original base as not a perfect power.\e[0m\n", b, x, a);
			return false;
		}
		if(base_a != base_b || exp_b*x != exp_a){
			fprintf(stderr, "\e[1;31mis_perfect_power incompatible results for %"PRIu64"**%"PRIu64" = %"PRIu64" and the original base.\e[0m\n", b, x, a);
			return false;
		}
	}else if(exp_a != x){
		fprintf(stderr, "\e[1;31mis_perfect_power got wrong exponent (but correct base) for %"PRIu64"**%"PRIu64" = %"PRIu64".\e[0m\n", b, x, a);
		return false;
	}else if(nut_u64_is_perfect_power(b, 31, &base_b, &exp_b)){
		fprintf(stderr, "\e[1;31mis_perfect_power inconsistent behavior for %"PRIu64"**%"PRIu64" = %"PRIu64" and the original base.\e[0m\n", b, x, a);
		return false;
	}
	return true;
}

int main(){
	uint64_t correct, tested;
	fprintf(stderr, "\e[1;34mTesting nth root and perfect power for cubes and higher powers...\e[0m\n");
	for(uint64_t x = 3, prev_r = 0; x <= 64; ++x){
		uint64_t r = nut_u64_nth_root(UINT64_MAX, x);
		if(powl(r, x) > (long double)UINT64_MAX || powl(r + 1, x) <= (long double)UINT64_MAX){
			fprintf(stderr, "\e[1;31m%"PRIu64" is not the largest integer <= %"PRIu64"-%s root of UINT64_MAX.\e[0m\n", r, x, get_ord_suffix(x));
			continue;
		}
		if(r == prev_r || r < 2){
			continue;
		}
		correct = 0, tested = 0;
		for(uint64_t b = 2; b <= r; ++b){
			++tested;
			if(check_power(b, x)){
				++correct;
			}
		}
		fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" roots correct for %"PRIu64"-%s powers)\e[0m\n", correct == tested ? "\e[1;32mPASSED" : "\e[1;31mFAILED", correct, tested, x, get_ord_suffix(x));
	}
	correct = 0, tested = 0;
	for(uint64_t b = 2; b <= 1000000; ++b){
		++tested;
		if(check_power(b, 2)){
			++correct;
		}
	}
	for(uint64_t b = (1ull << 32) - 1000000; b < (1ull << 32); ++b){
		++tested;
		if(check_power(b, 2)){
			++correct;
		}
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" roots correct for squares)\e[0m\n", correct == tested ? "\e[1;32mPASSED" : "\e[1;31mFAILED", correct, tested);
}

