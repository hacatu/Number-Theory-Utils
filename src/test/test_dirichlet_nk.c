#include <nut/modular_math.h>
#include <nut/debug.h>
#include <nut/sieves.h>
#include <nut/dirichlet.h>
#include <nut/dirichlet_powerful.h>
#include <stdint.h>

static void test_dirichlet_Nk(uint64_t k){
	nut_Diri nk_tbl [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri_init(&nk_tbl, 1000, 0);
	nut_Diri_compute_Nk(&nk_tbl, k, 0);
	int64_t Nk = 0;
	bool passed = true;
	for(int64_t n = 1; n <= nk_tbl.x; ++n){
		int64_t nk = nut_u64_pow(n, k);
		Nk += nk;
		if(n <= nk_tbl.y){
			if(nk != nut_Diri_get_dense(&nk_tbl, n)){
				fprintf(stderr, "\e[1;31mNk_tbl.dense[%"PRIi64"] is wrong\e[0m\n", n);
				passed = false;
			}
		}else if((nk_tbl.x/(n + 1) != nk_tbl.x/n) || n == nk_tbl.x){
			if(Nk != nut_Diri_get_sparse(&nk_tbl, nk_tbl.x/n)){
				fprintf(stderr, "\e[1;31mNk_tbl.sparse[%"PRIi64"] is wrong\e[0m\n", nk_tbl.x/n);
				passed = false;
			}
		}
	}
	if(passed){
		fprintf(stderr, "\e[1;32mPASSED test_dirichlet_Nk(%"PRIu64")\e[0m\n", k);
	}
}

static void test_convdiv(){
	uint64_t n = 100;
	int64_t modulus = 0;
	nut_Diri diri_1 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri diri_2 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri diri_3 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri diri_4 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri_init(&diri_1, n, 0);
	nut_Diri_init(&diri_2, n, 0);
	nut_Diri_init(&diri_3, n, 0);
	nut_Diri_init(&diri_4, n, 0);
	nut_Diri_compute_dk(&diri_1, 5, modulus, &diri_2, &diri_3);
	nut_Diri_compute_dk(&diri_2, 3, modulus, &diri_3, &diri_4);
	nut_Diri_convdiv(&diri_3, modulus, &diri_1, &diri_2);
	nut_Diri_compute_u(&diri_1, modulus);
	nut_Diri_compute_conv_u(&diri_2, modulus, &diri_1); // now diri_2 is d2, which should be the same as diri_3
	bool passed = true;
	for(int64_t n = 1; n <= diri_3.y; ++n){
		int64_t correct_val = nut_Diri_get_dense(&diri_2, n);
		int64_t val = nut_Diri_get_dense(&diri_3, n);
		if(val != correct_val){
			fprintf(stderr, "\e[1;31m(d5 </> d3).dense[%"PRIi64"] is wrong (got %"PRIi64", expected %"PRIi64")\e[0m\n", n, val, correct_val);
			passed = false;
		}
	}
	for(int64_t i = diri_3.yinv - 1; i >= 1; --i){
		int64_t correct_val = nut_Diri_get_sparse(&diri_2, i);
		int64_t val = nut_Diri_get_sparse(&diri_3, i);
		if(val != correct_val){
			fprintf(stderr, "\e[1;31m(d5 </> d3).sparse[%"PRIi64"] is wrong (got %"PRIi64", expected %"PRIi64")\e[0m\n", i, val, correct_val);
			passed = false;
		}
	}
	if(passed){
		fprintf(stderr, "\e[1;32mPASSED test_convdiv\e[0m\n");
	}
}

static void test_convdiv_2(){uint64_t n = 100;
	int64_t modulus = 1001961001;
	nut_Diri diri_1 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri diri_2 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri diri_3 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri diri_4 [[gnu::cleanup(nut_Diri_destroy)]];
	nut_Diri_init(&diri_1, n, 0);
	nut_Diri_init(&diri_2, n, 0);
	nut_Diri_init(&diri_3, n, 0);
	nut_Diri_init(&diri_4, n, 0);
	nut_Diri_compute_Nk(&diri_1, 7, modulus);
	nut_Diri_compute_Nk(&diri_2, 4, modulus);
	nut_Diri_compute_conv(&diri_3, modulus, &diri_1, &diri_2);
	nut_Diri_compute_Nk(&diri_1, 3, modulus);
	nut_Diri_convdiv(&diri_2, modulus, &diri_3, &diri_1);
	nut_Diri_compute_conv(&diri_4, modulus, &diri_1, &diri_2);
	bool passed = true;
	for(int64_t n = 1; n <= diri_3.y; ++n){
		int64_t correct_val = nut_Diri_get_dense(&diri_3, n);
		int64_t val = nut_Diri_get_dense(&diri_4, n);
		if(val != correct_val){
			fprintf(stderr, "\e[1;31m(N^7 <*> N^4 </> N^3 <*> N^3).dense[%"PRIi64"] is wrong (got %"PRIi64", expected %"PRIi64")\e[0m\n", n, val, correct_val);
			passed = false;
		}
	}
	for(int64_t i = diri_3.yinv - 1; i >= 1; --i){
		int64_t correct_val = nut_Diri_get_sparse(&diri_3, i);
		int64_t val = nut_Diri_get_sparse(&diri_4, i);
		if(val != correct_val){
			fprintf(stderr, "\e[1;31m(N^7 <*> N^4 </> N^3 <*> N^3).sparse[%"PRIi64"] is wrong (got %"PRIi64", expected %"PRIi64")\e[0m\n", i, val, correct_val);
			passed = false;
		}
	}
	if(passed){
		fprintf(stderr, "\e[1;32mPASSED test_convdiv_2\e[0m\n");
	}
}

static int64_t pp_h_fn(uint64_t p, uint64_t pp, uint64_t e, uint64_t m){
	uint64_t res;
	if(e < 2){
		return 1 - e;
	}else if(e == 2){
		res = pp - 1;
	}else{
		res = pp - pp/p;
	}
	return m ? res%m : res;
}

static int64_t pp_f_fn(uint64_t p, uint64_t pp, uint64_t e, uint64_t m){
	return e > 1 ? m ? pp%m : pp : 1;
}

static int64_t pp_g_fn(uint64_t p, uint64_t pp, uint64_t e, uint64_t m){
	return 1;
}

static void test_pp_pf_it(){
	nut_PfIt pf_it_h_fn [[gnu::cleanup(nut_PfIt_destroy)]];
	nut_PfIt pf_it_h_seqs [[gnu::cleanup(nut_PfIt_destroy)]];
	nut_PfIt_init_fn(&pf_it_h_fn, 10, 0, 0, pp_h_fn);
	nut_PfIt_init_hseqs(&pf_it_h_seqs, 10, 0, 0, pp_f_fn, pp_g_fn);
	while(1){
		nut_PfStackEnt a, b;
		bool status1 = nut_PfIt_next(&pf_it_h_fn, &a), status2 = nut_PfIt_next(&pf_it_h_seqs, &b);
		if(status1 != status2){
			fprintf(stderr, "\e[1;31mIterators did not end at the same time!\e[0m\n");
			break;
		}
		if(!status1){
			break;
		}
		if(b.n != a.n || b.hn != a.hn){
			fprintf(stderr, "\e[1;31mIterators returned different values!\e[0m\n");
			break;
		}
	}
}

int main(){
	test_dirichlet_Nk(4);
	test_convdiv();
	test_convdiv_2();
	test_pp_pf_it();
}

