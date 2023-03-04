#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>

#include <nut/sieves.h>

int main(){
	fprintf(stderr, "\e[1;34mSieving Mobius function up to 1000000...\e[0m\n");
	uint64_t *mobius = sieve_mobius(1000000);
	fprintf(stderr, "\e[1;34mAccumulating Mertens function up to 1000000...\e[0m\n");
	int64_t *mertens = compute_mertens_range(1000000, mobius);
	free(mobius);
	if(mertens[1000000] == 212){
		fprintf(stderr, "\e[1;32mPASSED\e[0m\n");
	}else{
		fprintf(stderr, "\e[1;31mFAILED\e[0m\n");
	}
	free(mertens);
}

