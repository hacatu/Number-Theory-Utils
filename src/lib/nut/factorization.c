#include <stddef.h>
#include <sys/random.h>
#include <string.h>
#include <stdbool.h>

#include <nut/modular_math.h>
#include <nut/factorization.h>

static inline void cleanup_free(void *_p){
	free(*(void**)_p);
	*(void**)_p = NULL;
}

nut_Factors *nut_make_Factors_w(uint64_t max_primes){
	static const nut_Factors dummy;
	return calloc(1, sizeof(nut_Factors) + sizeof(dummy.factors[0])*max_primes);
}

nut_Factors *nut_make_Factors_ub(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes]){
	uint64_t p = 1, w = 0;
	for(uint64_t i = 0; i < num_primes; ++i){
		p *= primes[i];
		if(p > n){
			return nut_make_Factors_w(w);
		}
		++w;
	}
	return NULL;
}

nut_Factors *nut_Factors_copy(const nut_Factors *factors){
	nut_Factors *ret = nut_make_Factors_w(factors->num_primes);
	if(ret){
		memcpy(ret, factors, offsetof(nut_Factors, factors) + factors->num_primes*sizeof(factors->factors[0]));
	}
	return ret;
}

uint64_t nut_Factors_prod(const nut_Factors *factors){
	uint64_t r = 1;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		r *= nut_u64_pow(factors->factors[i].prime, factors->factors[i].power);
	}
	return r;
}

uint64_t nut_Factor_divcount(const nut_Factors *factors){
	uint64_t s = 1;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		s *= factors->factors[i].power + 1;
	}
	return s;
}

uint64_t nut_Factor_divsum(const nut_Factors *factors){
	uint64_t s = 1;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		uint64_t p = factors->factors[i].prime;
		uint64_t a = factors->factors[i].power;
		s *= (nut_u64_pow(p, a + 1) - 1)/(p - 1);
	}
	return s;
}

uint64_t nut_Factor_divpowsum(const nut_Factors *factors, uint64_t power){
	if(power == 0){
		return nut_Factor_divcount(factors);
	}else if(power == 1){
		return nut_Factor_divsum(factors);
	}
	uint64_t s = 1;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		uint64_t p = factors->factors[i].prime;
		uint64_t a = factors->factors[i].power;
		s *= (nut_u64_pow(p, (a + 1)*power) - 1)/(nut_u64_pow(p, power) - 1);
	}
	return s;
}

uint64_t nut_Factor_divtupcount(const nut_Factors *factors, uint64_t k){
	if(k == 0){
		return !factors->num_primes;
	}else if(k == 1){
		return 1;
	}else if(k == 2){
		return nut_Factor_divcount(factors);
	}
	uint64_t s = 1;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		s *= nut_u64_binom(factors->factors[i].power + k - 1, k - 1);
	}
	return s;
}

void nut_Factor_ipow(nut_Factors *factors, uint64_t power){
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		factors->factors[i].power *= power;
	}
}

uint64_t nut_Factor_phi(const nut_Factors *factors){
	uint64_t s = 1;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		uint64_t p = factors->factors[i].prime;
		uint64_t a = factors->factors[i].power;
		s *= (nut_u64_pow(p, a - 1))*(p - 1);
	}
	return s;
}

uint64_t nut_Factor_carmichael(const nut_Factors *factors){
	uint64_t s = 1;
	if(!factors->num_primes){
		return s;
	}
	uint64_t p = factors->factors[0].prime;
	uint64_t a = factors->factors[0].power;
	if(p == 2){
		if(a >= 3){
			s = 1ull << (a - 2);
		}else{
			s = 1ull << (a - 1);
		}
	}else{
		s = nut_u64_pow(p, a - 1)*(p - 1);
	}
	for(uint64_t i = 1; i < factors->num_primes; ++i){
		p = factors->factors[i].prime;
		a = factors->factors[i].power;
		uint64_t phi_pk = (nut_u64_pow(p, a - 1))*(p - 1);
		s = nut_i64_lcm(s, phi_pk);
	}
	return s;
}

int nut_Factor_forall_divs_tmptmp(const nut_Factors *factors, int (*f)(const nut_Factors*, uint64_t, void*), void *data){
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
	nut_Factors *dfactors __attribute__((cleanup(cleanup_free))) = nut_Factors_copy(factors);
	nut_Factors *pfactors __attribute__((cleanup(cleanup_free))) = nut_Factors_copy(factors);
//#pragma GCC pop
	return nut_Factor_forall_divs(factors, f, data, dfactors, pfactors);
}

int nut_Factor_forall_divs(const nut_Factors *factors, int (*f)(const nut_Factors*, uint64_t, void*), void *data, nut_Factors *dfactors, nut_Factors *pfactors){
	dfactors->num_primes = factors->num_primes;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		dfactors->factors[i].prime = factors->factors[i].prime;
		dfactors->factors[i].power = 0;
		pfactors->factors[i].prime = 1;
	}
	uint64_t d = 1;
	while(1){
		if(f(dfactors, d, data)){
			return 1;
		}
		uint64_t i;
		for(i = 0; i < factors->num_primes; ++i){
			if(dfactors->factors[i].power < factors->factors[i].power){
				++dfactors->factors[i].power;
				pfactors->factors[i].prime *= dfactors->factors[i].prime;
				d *= dfactors->factors[i].prime;
				break;
			}
			dfactors->factors[i].power = 0;
			d /= pfactors->factors[i].prime;
			pfactors->factors[i].prime = 1;
		}
		if(i == factors->num_primes){
			return 0;
		}
	}
}

int nut_Factor_forall_divs_le_tmptmp(const nut_Factors *factors, uint64_t d_max, int (*f)(const nut_Factors*, uint64_t, void*), void *data){
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
	nut_Factors *dfactors __attribute__((cleanup(cleanup_free))) = nut_Factors_copy(factors);
	nut_Factors *pfactors __attribute__((cleanup(cleanup_free))) = nut_Factors_copy(factors);
//#pragma GCC pop
	return nut_Factor_forall_divs_le(factors, d_max, f, data, dfactors, pfactors);
}

int nut_Factor_forall_divs_le(const nut_Factors *factors, uint64_t d_max, int (*f)(const nut_Factors*, uint64_t, void*), void *data, nut_Factors *dfactors, nut_Factors *pfactors){
	dfactors->num_primes = factors->num_primes;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		dfactors->factors[i].prime = factors->factors[i].prime;
		dfactors->factors[i].power = 0;
		pfactors->factors[i].prime = 1;
	}
	uint64_t d = 1;
	if(d_max <= 1){
		return d_max && f(dfactors, d, data);
	}
	while(1){
		if(f(dfactors, d, data)){
			return 1;
		}
		uint64_t i;
		for(i = 0; i < factors->num_primes; ++i){
			if(dfactors->factors[i].power < factors->factors[i].power){
				++dfactors->factors[i].power;
				pfactors->factors[i].prime *= dfactors->factors[i].prime;
				d *= dfactors->factors[i].prime;
				if(d <= d_max){
					break;
				}
			}
			dfactors->factors[i].power = 0;
			d /= pfactors->factors[i].prime;
			pfactors->factors[i].prime = 1;
		}
		if(i == factors->num_primes){
			return 0;
		}
	}
}

void nut_Factor_append(nut_Factors *factors, uint64_t m, uint64_t k){
	for(uint64_t i = 0;; ++i){
		if(i == factors->num_primes){
			factors->factors[i].prime = m;
			factors->factors[i].power = k;
			++factors->num_primes;
			break;
		}else if(factors->factors[i].prime == m){
			factors->factors[i].power += k;
			break;
		}else if(factors->factors[i].prime > m){
			memmove(factors->factors + i + 1, factors->factors + i, (factors->num_primes - i)*sizeof(*factors->factors));
			++factors->num_primes;
			factors->factors[i].prime = m;
			factors->factors[i].power = k;
			break;
		}
	}
}

void nut_Factor_combine(nut_Factors *factors, const nut_Factors *factors2, uint64_t k){
	nut_Factors *factors3 = nut_make_Factors_w(factors->num_primes + factors2->num_primes);
	uint64_t i = 0, i2 = 0, i3 = 0;
	while(i < factors->num_primes && i2 < factors2->num_primes){
		if(factors->factors[i].prime < factors2->factors[i2].prime){
			factors3->factors[i3++] = factors->factors[i++];
		}else if(factors->factors[i].prime > factors2->factors[i2].prime){
			factors3->factors[i3].prime = factors2->factors[i2].prime;
			factors3->factors[i3++].power = factors2->factors[i2++].power*k;
		}else{
			factors3->factors[i3].prime = factors->factors[i].prime;
			factors3->factors[i3++].power = factors->factors[i++].power + factors2->factors[i2++].power*k;
		}
	}
	while(i < factors->num_primes){
		factors3->factors[i3++] = factors->factors[i++];
	}
	while(i2 < factors2->num_primes){
		factors3->factors[i3].prime = factors2->factors[i2].prime;
		factors3->factors[i3++].power = factors2->factors[i2++].power*k;
	}
	//TODO: elide copying initial factors from factors to factors3 and back to factors
	memcpy(factors->factors, factors3->factors, i3*sizeof(*factors->factors));
	factors->num_primes = i3;
	free(factors3);
}

int nut_Factor_fprint(FILE *file, const nut_Factors *factors){
	int res = 0;
	for(uint64_t i = 0; i < factors->num_primes; ++i){
		uint64_t power = factors->factors[i].power;
		if(!power){
			continue;
		}
		if(res){
			res += fprintf(file, "*");
		}
		res += fprintf(file, "%"PRIu64, factors->factors[i].prime);
		if(power != 1){
			res += fprintf(file, "^%"PRIu64, power);
		}
	}
	if(!res){
		res += fprintf(file, "1");
	}
	return res;
}

uint64_t nut_u64_powmod(uint64_t b, uint64_t e, uint64_t n){
	uint64_t r = 1;
	b %= n;
	while(e){
		if(e&1){
			r = (uint128_t)r*b%n;
		}
		e >>= 1;
		b = (uint128_t)b*b%n;
	}
	return (uint64_t)r;
}

int nut_u64_is_prime_dmr(uint64_t n){
	//static const uint64_t DMR_PRIMES[7] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};//sufficient for 64 bit numbers
	//had to disable 7 base check because the last base is too large to be squared under 64 bit multiplication
	static const uint64_t DMR_PRIMES[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};//sufficient for 64 bit numbers
	//static const uint64_t DMR_PRIMES_C = 7;
	static const uint64_t DMR_PRIMES_C = 12;
	uint64_t s, d;//s, d | 2^s*d = n - 1
	if(n%2 == 0){
		return n == 2;
	}
	--n;
	s = __builtin_ctz(n);
	d = n>>s;
	++n;
	for(uint64_t i = 0, a, x; i < DMR_PRIMES_C; ++i){
		a = DMR_PRIMES[i];
		if(a >= n){
			break;
		}
		x = nut_u64_powmod(a, d, n);
		if(x == 1 || x == n - 1){
			goto CONTINUE_WITNESSLOOP;
		}
		for(a = 0; a < s - 1; ++a){
			x = nut_u64_powmod(x, 2, n);
			if(x == 1){
				return 0;
			}
			if(x == n - 1){
				goto CONTINUE_WITNESSLOOP;
			}
		}
		return 0;
		CONTINUE_WITNESSLOOP:;
	}
	return 1;
}

/*
int is_prime_ecam(uint64_t n){
	static const uint64_t nDs[13] = {3, 4, 7, 8, 11, 12, 16, 19, 27, 28, 43, 67, 163};
	static const uint64_t nDs_len = 13;
	//Check if n is prime using a simplified version of the Atkin-Morain elliptic curve prime proving algorithm
	//First, we need to find a quadratic residue D so that a^2 + |D|b^2 = 4n and D is in (negative) nDs.
	//We can use Cornacchia's algorithm to help solve a^2 + |D|b^2 = 4n.
	//The number of potential solutions is very limited: a solution can be a primitive solution, or a primitive
	//solution to a^2 + |D|b^2 = n, since n is supposed to be prime and hence squarefree.
	//In the first case, there are at most 4 square roots of -|D| mod 4n and hence 2 starting points for Cornacchia.
	//In the second case, there are at most 2 square roots of -|D| mod n and hence 1 starting point for Cornacchia.
	//We can obtain a solution for each of these three starting points, although some may not lead to a solution.
	//Any solution a, b will have up to 3 associated solutions by negating components, thus there are at most
	//12 curves for every D, for at most 156 curves total (although it is very unlikely all Ds will be quadratic residues).
	//
	//Now, supposing we have a and b, we need to actually find a point on the curve.  But for this,
	//trying random x until x^3 + ax + b is a quadratic residue mod N and then pick y as a corresponding root.
	//We have that m = n + 1 - a is the number of points on the curve.
	//So the final thing to check is that we can factor m quickly and into a good form.
	//We want a prime divisor q of m such that q > (n^(1/4)+1)^2.
	//If we find one, demonstrating a point P such that [m/q]P is not identity but [m]P is identity proves
	//q is prime implies n is prime.
	//
}
*/

uint64_t nut_u64_rand(uint64_t a, uint64_t b){
	uint64_t l = b - a, r = 0, bytes = (71 - __builtin_clzll(l))/8;
	uint64_t ub;
	if(bytes == 8){
		ub = 0x7FFFFFFFFFFFFFFFull%l*-2;
	}else{
		ub = 1ull << (bytes*8);
		ub -= ub%l;
	}
	do{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
		getrandom(&r, bytes, 0);
#pragma GCC diagnostic pop
	}while(r >= ub);
	return r%l + a;
}

uint64_t nut_u64_prand(uint64_t a, uint64_t b){
	return nut_u64_rand(a, b);//TODO: test if this is a bottleneck, we don't need calls to this function to be secure
}

const uint64_t nut_small_primes[25] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

const nut_FactorConf nut_default_factor_conf = {
	.pollard_max= 100000,    //maximum number to use Pollard's Rho algorithm for
	.pollard_stride= 10,     //number of gcd operations to coalesce, decreases time for a single iteration at the cost of potentially doing twice this many extra iterations
	.lenstra_max= UINT64_MAX,//maximum number to use Lenstra's Elliptic Curve algorithm for
	.lenstra_bfac= 10        //roughly speaking, the number of iterations to try before picking a new random point and curve
};

uint64_t nut_u64_factor_trial_div(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes], nut_Factors *factors){
	factors->num_primes = 0;
	for(uint64_t i = 0; i < num_primes; ++i){
		uint64_t p = primes[i];
		if(n%p == 0){
			uint64_t j = factors->num_primes++;
			factors->factors[j].prime = p;
			factors->factors[j].power = 1;
			n /= p;
			while(n%p == 0){
				++factors->factors[j].power;
				n /= p;
			}
		}
		if(p*p > n){
			if(n > 1){
				uint64_t j = factors->num_primes++;
				factors->factors[j].prime = n;
				factors->factors[j].power = 1;
			}
			return 1;
		}
	}
	return n;
}

uint64_t nut_u64_factor_heuristic(uint64_t n, uint64_t num_primes, const uint64_t primes[static num_primes], const nut_FactorConf *conf, nut_Factors *factors){
	n = nut_u64_factor_trial_div(n, num_primes, primes, factors);
	if(n == 1){
		return 1;
	}
	uint64_t exponent = 1;
	nut_u64_is_perfect_power(n, 9, &n, &exponent);
	if(nut_u64_is_prime_dmr(n)){//TODO: get a better primality test and possibly run it later (ie after some pollard-rho)
		nut_Factor_append(factors, n, exponent);
		return 1;
	}
	uint64_t smoothness = 101*101;
	nut_Factors *factors2 = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	uint64_t m;
	while(1){
		if(n <= conf->pollard_max){//TODO: allow iteration count based stopping of pollard-rho brent so we can get small factors of big numbers that way?
			do{
				uint64_t x = nut_u64_prand(0, n);
				m = nut_u64_factor1_pollard_rho_brent(n, x, conf->pollard_stride);
			}while(m == n);
			uint64_t k = 1;
			n /= m;
			while(n%m == 0){
				k += 1;
				n /= m;
			}
			if(m < smoothness || nut_u64_is_prime_dmr(m)){
				nut_Factor_append(factors, m, k*exponent);
			}else{
				factors2->num_primes = 0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
				m = nut_u64_factor_heuristic(m, 0, NULL, conf, factors2);
#pragma GCC diagnostic pop
				if(m != 1){
					free(factors2);
					return m;//TODO: abort
				}
				nut_Factor_combine(factors, factors2, k*exponent);
			}
			if(n == 1){
				free(factors2);
				return 1;
			}
			if(n < smoothness || nut_u64_is_prime_dmr(n)){
				nut_Factor_append(factors, n, exponent);
				free(factors2);
				return 1;
			}
		}else if(n <= conf->lenstra_max){
			do{
				uint64_t x = nut_u64_prand(0, n);
				uint64_t y = nut_u64_prand(0, n);
				uint64_t a = nut_u64_prand(0, n);
				m = nut_u64_factor1_lenstra(n, x, y, a, conf->lenstra_bfac);
			}while(m == n);
			uint64_t k = 1;
			n /= m;
			while(n%m == 0){
				k += 1;
				n /= m;
			}
			if(m < smoothness || nut_u64_is_prime_dmr(m)){
				nut_Factor_append(factors, m, k*exponent);
			}else{
				factors2->num_primes = 0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
				m = nut_u64_factor_heuristic(m, 0, NULL, conf, factors2);
#pragma GCC diagnostic pop
				if(m != 1){
					free(factors2);
					return m;//TODO: abort
				}
				nut_Factor_combine(factors, factors2, k*exponent);
			}
			if(n == 1){
				free(factors2);
				return 1;
			}
			if(n < smoothness || nut_u64_is_prime_dmr(n)){
				nut_Factor_append(factors, n, exponent);
				free(factors2);
				return 1;
			}
		}else{//TODO: implement quadratic sieve and number field sieve
			free(factors2);
			return n;
		}
	}
}

uint64_t nut_u64_nth_root(uint64_t a, uint64_t n){
	if(a < 2){
		return a;
	}
	uint64_t x;
	if(n < 13){
		uint64_t floor_log_2_a = 63 - __builtin_clzll(a);
		x = 1ull << (floor_log_2_a/n);
		if(x == 1){
			return 1;
		}
	}else for(x = 1;; ++x){
		uint128_t x_pow = nut_u128_pow(x + 1, n);
		if(x_pow > a){
			return x;
		}
	}
	for(uint64_t i = 0;; ++i){
		uint128_t x_pow = nut_u128_pow(x, n - 1);
		uint64_t x_next = ((n - 1)*x_pow*x + a)/(n*x_pow);
		if(i > 1 && x_next >= x){
			return x;
		}
		x = x_next;
	}
}

bool nut_u64_is_perfect_power(uint64_t a, uint64_t max, uint64_t *_base, uint64_t *_exp){
	if(a < 2){
		return false;
	}
	/// the maximum exponent ever for a uint64 is 63, but if the exponent is composite we can find it in stages,
	/// so the maximum exponent we will ever have to take is 61.
	uint64_t x = 1;
	for(uint64_t i = 0; i < 18; ++i){
		uint64_t p = nut_small_primes[i];
		if(p > max){
			break;
		}
		while(1){
			uint64_t r = nut_u64_nth_root(a, p);
			if(nut_u64_pow(r, p) != a){
				break;
			}
			x *= p;
			a = r;
			max /= p;
			if(p > max){
				goto OUTPUT_RESULT;
			}
		}
	}
	OUTPUT_RESULT:;
	if(x != 1){
		*_base = a;
		*_exp = x;
		return true;
	}
	return false;
}

uint64_t nut_u64_factor1_pollard_rho(uint64_t n, uint64_t x){
	uint64_t y = x, d = 1;
	while(d == 1){
		x = (x*x + 1)%n;
		y = (y*y + 1)%n;
		y = (y*y + 1)%n;
		d = nut_i64_egcd(x > y ? x - y : y - x, n, NULL, NULL);
	}
	return d;
}

uint64_t nut_u64_factor1_pollard_rho_brent(uint64_t n, uint64_t x, uint64_t m){
	uint64_t y = x, ys = x;
	uint64_t d = 1;//gcd of n and the difference of the terms in the sequence
	uint64_t r = 1;//power of two in brent's algorithm
	uint64_t q = 1;//product of all differences so far, allowing us to process m at a time when checking for a nontrival factor d
	while(d == 1){
		x = y;
		for(uint64_t i = 0; i < r; ++i){
			y = (y*y + 1)%n;
		}
		for(uint64_t k = 0; k < r && d == 1; k += m){
			ys = y;
			for(uint64_t i = 0; i < m && i < r - k; ++i){
				y = (y*y + 1)%n;
				q = q*(x > y ? x - y : y - x)%n;
			}
			d = nut_i64_egcd(q, n, NULL, NULL);
		}
		r *= 2;
	}
	if(d == n){
		do{
			ys = (ys*ys + 1)%n;
			d = nut_i64_egcd(x > ys ? x - ys : ys - x, n, NULL, NULL);
		}while(d == 1);
	}
	return d;
}

//convenience function to double a point on an elliptic curve.  _xr and _yr are out params
//returns 0 for an ordinary point, 1 for identity, and -1 if a nontrivial factor was found and placed in _xr
static inline int ecg_double(int64_t n, int64_t a, int64_t x, int64_t y, int is_id, int64_t *_xr, int64_t *_yr){
	if(is_id || y == 0){
		return 1;
	}
	int64_t s, dy, dx;
	dy = nut_i64_mod(3*x*x + a, n);
	dx = nut_i64_mod(2*y, n);
	int64_t d = nut_i64_egcd(dx, n, &s, NULL);
	if(d != 1){
		*_xr = d;
		return -1;
	}
	s = s*dy%n;
	*_xr = nut_i64_mod(s*s - 2*x, n);
	*_yr = nut_i64_mod(y + s*(*_xr - x), n);
	return 0;
}

//convenience function to add two points on an elliptic curve.  _xr and _yr are out params
//returns 0 for an ordinary point, 1 for identity, and -1 if a nontrivial factor was found and placed in _xr
static inline int ecg_add(int64_t n, int64_t a, int64_t xp, int64_t yp, int is_id_p, int64_t xq, int64_t yq, int is_id_q, int64_t *_xr, int64_t *_yr){
	if(is_id_p){
		if(is_id_q){
			return 1;
		}
		*_xr = xq, *_yr = yq;
		return 0;
	}else if(is_id_q){
		*_xr = xp, *_yr = yp;
		return 0;
	}
	int64_t s, dy, dx;
	if(xp != xq){
		dy = nut_i64_mod(yp - yq, n);
		dx = nut_i64_mod(xp - xq, n);
	}else if(yp != yq || yp == 0){
		return 1;
	}else{
		dy = nut_i64_mod(3*xp*xp + a, n);
		dx = nut_i64_mod(2*yp, n);
	}
	int64_t d = nut_i64_egcd(dx, n, &s, NULL);
	if(d != 1){
		*_xr = d;
		return -1;
	}
	s = s*dy%n;
	*_xr = nut_i64_mod(s*s - xp - xq, n);
	*_yr = nut_i64_mod(yp + s*(*_xr - xp), n);
	return 0;
}

//convenience function to add k copies of a point on an elliptic curve together.  _xr and _yr are out params
//returns 0 for an ordinary point, 1 for identity, and -1 if a nontrivial factor was found and placed in _xr
static inline int ecg_scale(int64_t n, int64_t a, int64_t x, int64_t y, int is_id, int64_t k, int64_t *_xr, int64_t *_yr){
	if(is_id || !k){
		return 1;
	}
	int is_id_r = 0, is_id_s = 0;
	int64_t xs = x, ys = y;
	while(!(k&1)){
		is_id_s = ecg_double(n, a, xs, ys, is_id_s, &xs, &ys);
		if(is_id_s == -1){
			*_xr = xs;
			return -1;
		}else if(is_id_s == 1){
			return 1;
		}
		k >>= 1;
	}
	*_xr = xs, *_yr = ys;
	is_id_r = is_id_s;
	while(1){
		k >>= 1;
		if(!k){
			return is_id_r;
		}
		is_id_s = ecg_double(n, a, xs, ys, is_id_s, &xs, &ys);
		if(is_id_s == -1){
			*_xr = xs;
			return -1;
		}else if(is_id_s == 1){
			return is_id_r;
		}else if(k&1){
			is_id_r = ecg_add(n, a, xs, ys, is_id_s, *_xr, *_yr, is_id_r, _xr, _yr);
			if(is_id_r == -1){
				return -1;
			}
		}
	}
}

//TODO: this whole function needs to be rewritten to use projective montgomery curves instead of affine wierstrass curves
//TODO: consider adding brent's second stage for fun
int64_t nut_u64_factor1_lenstra(int64_t n, int64_t x, int64_t y, int64_t a, int64_t B){
	int64_t d = nut_i64_egcd(n, 6, NULL, NULL);//ensure Z/nZ doesn't have a field of characteristic 2 or 3 as a subgroup
	if(d != 1){
		return d;
	}
	int64_t b = (x*x + a)%n;
	b = nut_i64_mod(y*y - x*b, n);
	d = a*a%n;
	d = 4*a*d%n;
	d = nut_i64_egcd(d + 27*(b*b%n), n, NULL, NULL);//ensure the rhs of the elliptic curve does not have a repeated zero, leading to a cusp
	if(d != 1){
		return d;
	}
	int is_id = 0;
	for(int64_t k = 2; k <= B; ++k){
		is_id = ecg_scale(n, a, x, y, is_id, k, &x, &y);
		if(is_id == -1){
			return x;
		}else if(is_id == 1){
			return n;
		}
	}
	return n;
}

/* There are many possible forms for elliptic curves, similar to how there are many possible forms for parabolas (y=ax^2+bx+c, y=a(x-l)(x-m), etc).
 * The basic form, used above, is called Wierstrass form and looks like y^2=x^3+ax+b.  Other common forms include twisted Edwards and Montgomery.
 * Montgomery form is by^2=x^3+ax^2+x.  Aside from what form we use for the curve, there are many possible ways to store points on the curve.
 * The simplest is to simply store x and y and a boolean indicating if the point is actually identity or not.  However, this creates branching
 * and complexity due to handling points at infinity and finite points on the curve separately.  So instead we can use projective coordinates,
 * where x and y are expressed as ratios x=X/Z and y=Y/Z.  If Z is zero, we have a point at infinity.  Variations of this exist where X and Y are
 * divided by different powers of Z.
 * 
 * Obviously, x and y are closely related by whatever curve we are using.  For any x, there are at most 2 solutions y placing it on the curve.
 * For any y, there are at most 3 solutions x putting it on the curve.  So it's pretty concievable that X and Z might be sufficient to compute
 * the x coordinate of any scalar multiple of a point on the curve.  In fact, it turns out this is true.  Using projective coordinates also lets
 * us do fewer gcd computations, hence the motivation to have a Z coordinate, yet we do not really need the y coordinate per se, so it is great
 * we can do away with it.
 * 
 * To compute a scalar multiple of a point P, we can use a double and add algorithm, just like binary exponentiation.  For explicit x, y coordinates
 * on a Wierstrass curve this is trivially easy.  For a Montgomery curve, we get the following rules for addition and doubling:
 * 
 * X_{m+n} = Z_{m-n}((X_m-Z_m)(X_n+Z_n)+(X_m+Z_m)(X_n-Z_n))^2
 * Z_{m+n} = X_{m-n}((X_m-Z_m)(X_n+Z_n)-(X_m+Z_m)(X_n-Z_n))^2
 * 
 * X_{2m} = (X_n+Z_n)^2(X_n-Z_n)^2
 * Z_{2m} = (4X_n Z_n)((X_n-Z_n)^2+((a+20)/4)(4X_n Z_n))
 * where
 * 4X_n Z_n = (X_n + Z_n)^2 - (X_n - Z_n)^2
 * 
 * The rule for doubling is great, but the rule for addition seems pretty much useless: we don't know a thing about X_{m-n} or Z_{m-n}!!!!
 * The simplest way around this, which is very good in our case where we need to compute kP for many arbitrary P and consecutive small k, is
 * to consider the case where m-n = 1 so the equation for P_{m+n} in terms of P_{m-n}, P_m, and P_n becomes an equation for P_{2n+1} in terms of
 * P_{n+1} and P_n.
 * 
 * Thus, we can build up P_k by scanning over the binary representation of k from most significant to least.  At each bit index, we have P_n and
 * P_{n+1} where n is the portion of k more significant than the current position, rounded down to be even, and from them we compute P_{2n+1}
 * and P_{2n}
 * 
 * Ex:      L=(0)P   H=(1)P
 * 10110    L=L+H=(1)P=(2*0+1)P
 * ^        H=2H=(2)P=(2*1)P
 * 10110    H=L+H=(3)P=(2*1+1)P
 *  ^       L=2L=(2)P=(2*1)P
 * 10110    L=L+H=(5)P=(2*2+1)P
 *   ^      H=2H=(6)P=(2*3)P
 * 10110    L=L+H=(11)P=(2*5+1)P
 *    ^     H=2H=(12)P=(2*6)P
 * 10110    H=L+H=(23)P=(2*11+1)P
 *     ^    L=2L=(22)P=(2*11)P
 * and the result is L
 */
int64_t nut_u64_factor1_lenstra_montgomery(int64_t n, int64_t x, int64_t y, int64_t a, int64_t B){
	int64_t d = nut_i64_egcd(n, 6, NULL, NULL);//check n does not contain a degenerate prime field
	if(d != 1){
		return d;
	}
	int64_t b = nut_i64_mod(y*y, n);
	d = nut_i64_egcd(b, n, &b, NULL);//check that x, y is on some montgomery curve mod n with a given
	if(d != 1){
		return d;
	}
	b = nut_i64_mod(b*x, n);
	d = nut_i64_mod(x + a, n);
	d = nut_i64_mod(x*d + 1, n);
	b = nut_i64_mod(b*d, n);
	d = nut_i64_mod(a*a - 4, n);
	d = nut_i64_egcd(b*d, n, NULL, NULL);//check that the curve does not have a sharp cusp
	if(d != 1){
		return d;
	}
	int64_t C;
	nut_i64_egcd(4, n, &C, NULL);
	C = nut_i64_mod((a + 2)*C, n);
	int64_t Zh = 1, Xh = x;
	int64_t Z1 = 1, X1 = x;
	for(int64_t k = 2; k <= B; ++k){
		int64_t Zl = 0, Xl = 1;
		for(int64_t t = 1ll << (63 - __builtin_clzll(k)); t; t >>= 1){
			if(k & t){
				//L = L + H
				//H = 2*H
				int64_t dh = nut_i64_mod(Xh - Zh, n);
				int64_t sl = nut_i64_mod(Xl + Zl, n);
				int64_t sh = nut_i64_mod(Xh + Zh, n);
				int64_t dl = nut_i64_mod(Xl - Zl, n);
				int64_t dhsl = nut_i64_mod(dh*sl, n);
				int64_t shdl = nut_i64_mod(sh*dl, n);
				Xl = nut_i64_mod(dhsl + shdl, n);
				Xl = nut_i64_mod(Xl*Xl, n);
				Xl = nut_i64_mod(Z1*Xl, n);
				Zl = nut_i64_mod(dhsl - shdl, n);
				Zl = nut_i64_mod(Zl*Zl, n);
				Zl = nut_i64_mod(X1*Zl, n);
				int64_t sh2 = nut_i64_mod(sh*sh, n);
				int64_t dh2 = nut_i64_mod(dh*dh, n);
				int64_t ch = nut_i64_mod(sh2 - dh2, n);
				Xh = nut_i64_mod(sh2*dh2, n);
				Zh = nut_i64_mod(C*ch, n);
				Zh = nut_i64_mod(ch*(dh2 + Zh), n);
			}else{
				//H = L + H
				//L = 2*L
				int64_t dh = nut_i64_mod(Xh - Zh, n);
				int64_t sl = nut_i64_mod(Xl + Zl, n);
				int64_t sh = nut_i64_mod(Xh + Zh, n);
				int64_t dl = nut_i64_mod(Xl - Zl, n);
				int64_t dhsl = nut_i64_mod(dh*sl, n);
				int64_t shdl = nut_i64_mod(sh*dl, n);
				Xh = nut_i64_mod(dhsl + shdl, n);
				Xh = nut_i64_mod(Xh*Xh, n);
				Zh = nut_i64_mod(dhsl - shdl, n);
				Zh = nut_i64_mod(Zh*Zh, n);
				Zh = nut_i64_mod(x*Zh, n);
				int64_t sl2 = nut_i64_mod(sl*sl, n);
				int64_t dl2 = nut_i64_mod(dl*dl, n);
				int64_t cl = nut_i64_mod(sl2 - dl2, n);
				Xl = nut_i64_mod(sl2*dl2, n);
				Zl = nut_i64_mod(C*cl, n);
				Zl = nut_i64_mod(cl*(dl2 + Zl), n);
			}
		}
		if(!Zl){
			return n;
		}
		d = nut_i64_egcd(Zl, n, NULL, NULL);
		if(d != 1){
			return d;
		}
		Z1 = Zl;
		X1 = Xl;
	}
	return n;
}

int64_t nut_i64_jacobi(int64_t n, int64_t k){
	int64_t j = 1;
	while(1){
		int64_t s = __builtin_ctzll(n);
		int64_t q = n >> s;
		if((s&1) && ((k&7) == 3 || ((k&7) == 5))){
			j = -j;
		}
		if(q == 1){
			return j;
		}else if(q == k - 1){
			return (k&3)==1 ? j : -j;
		}else if((q&2) && (k&2)){
			j = -j;
		}
		n = k%q;
		k = q;
	}
}

int64_t nut_i64_rand_nr(int64_t p){
	while(1){
		int64_t z = nut_u64_rand(2, p);
		if(nut_i64_jacobi(z, p) == -1){
			return z;
		}
	}
}

int64_t nut_i64_sqrt_shanks(int64_t n, int64_t p){
	int64_t s = __builtin_ctzll(p-1);
	int64_t q = p >> s;//p-1 = q*2^s
	int64_t z = nut_i64_rand_nr(p);
	//printf("trying \"nonresidue\" %"PRIu64"\n", z);
	int64_t m = s;
	int64_t c = nut_u64_powmod(z, q, p);
	int64_t t = nut_u64_powmod(n, q, p);
	int64_t r = nut_u64_powmod(n, (q + 1) >> 1, p);
	while(t != 1){
		int64_t i = 1;
		for(int64_t s = (int128_t)t*t%p; s != 1; s = (int128_t)s*s%p, ++i);
		int64_t b = c;
		for(int64_t j = 0; j < m - i - 1; ++j){
			b = (int128_t)b*b%p;
		}
		m = i;
		c = (int128_t)b*b%p;
		t = (int128_t)t*c%p;
		r = (int128_t)r*b%p;
	}
	return r;
}

int64_t nut_i64_sqrt_cipolla(int64_t n, int64_t p){
	int64_t a, w;
	do{
		a = nut_u64_rand(2, p);
		w = nut_i64_mod((int128_t)a*a%p - n, p);
	}while(nut_i64_jacobi(w, p) != -1);
	int64_t u_s = a, w_s = 1;
	int64_t u_r = 1, w_r = 0;
	for(int64_t k = (p + 1) >> 1; k; k >>= 1){
		if(k&1){
			int64_t _w_r = (int128_t)u_r*w_s%p;
			_w_r = (_w_r + (int128_t)w_r*u_s)%p;
			u_r = (int128_t)u_r*u_s%p;
			w_r = (int128_t)w_r*w_s%p;
			w_r = (int128_t)w_r*w%p;
			u_r = ((int128_t)u_r + w_r)%p;
			w_r = _w_r;
		}
		int64_t _w_s = (int128_t)2*u_s*w_s%p;
		u_s = (int128_t)u_s*u_s%p;
		w_s = (int128_t)w_s*w_s%p;
		w_s = (int128_t)w*w_s%p;
		u_s = ((int128_t)u_s + w_s)%p;
		w_s = _w_s;
	}
	return u_r;
}

int64_t nut_i64_sqrt_mod(int64_t n, int64_t p){
	int64_t r;
	if((p&3) == 3){
		r = nut_u64_powmod(n, (p + 1) >> 2, p);
	}else if((p&7) == 5){
		r = nut_u64_powmod(n, (p + 3) >> 3, p);
		if(r*r%p != n){
			r = (int128_t)r*nut_u64_powmod(2, (p - 1) >> 2, p)%p;
		}//can add 9 mod 16 case
	}else{
		int64_t m = 64 - __builtin_clzll(p);
		int64_t s = __builtin_ctzll(p - 1);
		if(8*m + 20 >= s*(s - 1)){
			r = nut_i64_sqrt_shanks(n, p);
		}else{
			r = nut_i64_sqrt_cipolla(n, p);
		}
	}
	return r;
}

