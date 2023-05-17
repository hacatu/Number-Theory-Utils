#include <stdio.h>
#include <inttypes.h>

#include <nut/factorization.h>

int main(){
	uint64_t trials = 1000000, passed = 0;
	fprintf(stderr, "\e[1;34mFinding roots of n**2 - n + 1 mod p for %"PRIu64" random primes...\e[0m\n", trials);
	for(uint64_t i = 0; i < trials; ++i){
		int64_t p = nut_u64_rand(2, 1ull << 30);
		while(!nut_u64_is_prime_dmr(p)){
			++p;
		}
		//Want to solve n**2 - n + 1 == 0 mod p
		//<=> 4*n**2 - 4*n + 4 == 0 mod p
		//<=> (2*n - 1)**2 == -3 mod p
		//<=> n == 2**-1*(1 +- sqrt(-3)) mod p
		switch(nut_i64_jacobi(p - 3, p)){
			case -1:
				//printf("n**2 - n + 1 has no roots mod %"PRId64"\n", p);
				break;
			case 0://-3 is a multiple of p, ie p == 3 so we only have 1 solution
				//printf("n**2 - n + 1 has a root at %"PRId64" mod %"PRId64"\n", (p + 1)/2, p);
				break;
			default: {
				//WARNING: here we know p - 3 is a quadratic residue, but nut_i64_sqrt_mod does not check this and thus if a nonresidue is given and
				//Tonelli-Shanks is selected as the optimal algorithm for this p, this would be an infinite loop
				int64_t r = nut_i64_sqrt_mod(p - 3, p);
				//printf("n**2 - n + 1 has roots at %"PRId64" and %"PRId64" mod %"PRId64"\n", mod((p + 1)/2*(1 + r), p), mod((p + 1)/2*(1 - r), p), p);
				int64_t n = nut_i64_mod((p + 1)/2*(1 + r), p);
				if(nut_i64_mod(n*n - n + 1, p)){
					fprintf(stderr, "\e[1;31mFound invalid root of -3 mod %"PRId64"\e[0m\n", p);
					continue;
				}
				n = nut_i64_mod((p + 1)/2*(1 - r), p);
				if(nut_i64_mod(n*n - n + 1, p)){
					fprintf(stderr, "\e[1;31mFound invalid root of -3 mod %"PRId64"\e[0m\n", p);
					continue;
				}
			}
		}
		++passed;
	}
	fprintf(stderr, "%s (%"PRIu64"/%"PRIu64" trials passed)\e[0m\n", passed == trials ? "\e[1;32mPASSED" : "\e[1;31mFAILED", passed, trials);
}