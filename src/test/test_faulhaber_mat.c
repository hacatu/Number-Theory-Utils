#include <nut/matrix.h>

int main(){
	nut_i64_Matrix a1, a2;
	nut_i64_Matrix_init(&a1, 8, 8);
	nut_i64_Matrix_init(&a2, 8, 8);
	nut_i64_Matrix_fill_short_pascal(&a1);
	int64_t denom = nut_i64_Matrix_invert_ltr(&a1, &a2);
	nut_i64_Matrix_destroy(&a1);
	int64_t sums[8] = {};
	for(int64_t k = 0; k <= 100; ++k){
		for(int64_t e = 0, ke = 1; e < 8; ++e, ke *= k){
			sums[e] += ke;
		}
	}
	int64_t vand_vec[8] = {};
	for(int64_t e = 1, ke = 101; e <= 8; ++e, ke *= 101){
		vand_vec[e - 1] = ke;
	}
	int64_t faulhaber_vec[8] = {};
	nut_i64_Matrix_mul_vec(&a2, vand_vec, faulhaber_vec);
	bool status = true;
	for(int64_t i = 0; i < 8; ++i){
		if(faulhaber_vec[i]%denom){
			fprintf(stderr, "\e[1;31mERROR: got non-integer power sum for %"PRIi64"\e[0m\n", i);
			status = false;
		}
		int64_t s = faulhaber_vec[i]/denom;
		if(s != sums[i]){
			fprintf(stderr, "\e[1;31mERROR: got %"PRIi64" instead of %"PRIi64" for sum of %"PRIi64" powers\e[0m\n", s, sums[i], i);
			status = false;
		}
	}
	if(status){
		fprintf(stderr, "\e[1;32mSUCCESS: Got sums of 0-7th powers of integers 0-100 correctly\e[0m\n");
	}
	nut_i64_Matrix_destroy(&a2);
}

