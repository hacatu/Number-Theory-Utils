# Number Theory Utils (NUT)

[![GitHub CI](https://github.com/hacatu/Number-Theory-Utils/actions/workflows/cov_and_docs.yml/badge.svg)](https://github.com/hacatu/Number-Theory-Utils/actions/workflows/cov_and_docs.yml)
[![Documentation](https://img.shields.io/badge/-documentation-gray)](https://hacatu.github.io/Number-Theory-Utils/docs)
[![Coverage](https://hacatu.github.io/Number-Theory-Utils/coverage.svg)](https://hacatu.github.io/Number-Theory-Utils/cov)

A collection of functions for working modular arithmetic, polynomials over finite fields, and related things.

Implements factorization of 64 bit numbers using trial division, Pollard's Rho algorithm with Brent or Floyd
cycle finding, and Lenstra's Elliptic Curve algorithm with Wierstrass or Montgomery curves.  Planned support
for arbitrary precision integers and Kraitcheck style methods (ie quadratic sieve and number field sieve),
although Flint and its extension library Arb may be more suited for your needs.

For polynomials, implements finding roots of a polynomial over a prime field using the Cantor-Zassenhaus algorithm.

Also includes general purpose modular arithmetic routines, including powers, Miller-Rabin prime checking,
random numbers, random quadratic nonresidues, Jacobi symbols, modular square root including Tonelli-Shanks
and Cippola, extended gcd, Chinese remainder theorem, and Euclidean remainder.

For polynomials, also includes arithmetic (addition, subtraction, scalar multiplication, coefficientwise multiplication,
multiplication, quotient and remainder, Horner evaluation, gcd, and powers), printing/normalization,
and distinct degree factorization.  Currently there is a simple driver function that will find all roots
of a polynomial over a prime field, but if a full factorization is desired the squarefree step must still be done manually.
Like the factorization part of the library, polynomials use 64 bit integers, so computations involving numbers larger than
about 2<sup>30</sup> potentially can fail due to overflow.

## Building

This library uses variant makefiles for each build type.  To build the release variant, simply run `make`
in the project root directory.  The resulting files will be in `build/release/lib` and `build/release/bin`,
or you can run `sudo make install` after `make` finishes.

To build a different variant, change the `$(BUILD_ROOT)` variable from `build/release` to `build/debug`
or some other value, for instance `make BUILD_ROOT=build/release`.  You can check what variants exist
by looking at the subdirectories of `build`, and you can create your own by copying and modifying `build/debug`.

The currently extant variants are
- `build/coverage`: debug variant, includes address sanitizer and some ub sanitizers.  builds using `gcc`.  default `make` target: `coverage`
- `build/debug`: debug variant, includes memory sanitizer and most ub sanitizers.  builds using `clang`.  default `make` target: `test`
- `build/release`: release variant, no sanitizers.  builds using `gcc` with `-O3`.  default `make` target: `all`
- `build/valgrind`: debug variant, no sanitizers.  builds using `clang` with `-O1`.  default `make` target: `test`.  passes `-g` to `test.py` to invoke the valgrind wrapper

The files in a variant build directory are arranged as follows: `Makefile` is the makefile, variants
probably won't have to modify this much aside from removing coverage information; `cflags.txt`,
`ldflags.txt`, and `makedeps_cflags.txt` contain configuration flags that should be passed to the compiler
and linker.  Flag files are used to consolidate compiler flags, reduce how often the makefile needs to be
tweaked, and clean up build logs.

Based on the source files in `include` and `src`, together with the configuration in a variant build directory,
output files are produced in the following subdirectories of `$(BUILD_ROOT)`: `bin` for executables, `bin/test`
for test executables, `lib` for static and dynamic libraries, `obj` for object files and dependency files,
`notes` for coverage information, `log` for test logs, and `cov` for human readable coverage reports.

Executables, static libraries, and test executables are automatically created based on the contents of `src`:
every C file or directory in `src/bin` is turned into its own executable in `$(BUILD_ROOT)/bin`, every C file
or directory in `src/lib` is turned into its own static library in `$(BUILD_ROOT)/lib`, and every C file or
directory in `src/test` is turned into its own test executable in `$(BUILD_ROOT)/bin/test`.

`make clean` should remove all output files in all variant build directories, as well as all generated
documentation.

Finally, `make docs` generates the documentation in the `docs` directory.  Like `make clean`, this is not
tied to a build variant and even if you specify one the same thing will happen.

`make debug_makefile` simply exists to facillitate printing make variables, don't worry about it.

Only Linux is properly supported.  To build, only `gcc`, `ar`, and `make` are strictly required, but `lcov`
and `doxygen` are required for coverage and documentation, and `python`, `clang`, and `valgrind` are required
for testing.

If running as root in an automated environment, you will need to either pass `ALLOW_ROOT=1` as an option to
`make`, or pipe `yes` into `make` to auto accept the safety prompt for running as root.

## Linking
Once the library has been built, it can be linked with C programs by adding the flags `-Lbuild/debug -lnut`
or `-L$(BUILD_ROOT) -lnut` to use a different variant build.  You can also install
it to your system libraries by doing `sudo make BUILD_ROOT=build/release install`.
This will place the libraries in `/usr/lib/nut/` and the headers in `/usr/include/nut/`.
Then the `-L` flag can be omitted and the library can be linked with simply `-lnut`.

## Examples

The tests in `src/test` offer a lot more interesting examples, with basically every function being at least used.

### Solve a quadratic equation mod some primes
```C
	for(uint64_t i = 0; i < 100; ++i){
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
				printf("n**2 - n + 1 has no roots mod %"PRId64"\n", p);
				break;
			case 0://-3 is a multiple of p, ie p == 3 so we only have 1 solution
				printf("n**2 - n + 1 has a root at %"PRId64" mod %"PRId64"\n", (p + 1)/2, p);
				break;
			default: {
				//WARNING: here we know p - 3 is a quadratic residue, but nut_i64_sqrt_mod does not check this and thus if a nonresidue is given and
				//Tonelli-Shanks is selected as the optimal algorithm for this p, this would be an infinite loop
				int64_t r = nut_i64_sqrt_mod(p - 3, p);
				printf("n**2 - n + 1 has roots at %"PRId64" and %"PRId64" mod %"PRId64"\n", nut_i64_mod((p + 1)/2*(1 + r), p), nut_i64_mod((p + 1)/2*(1 - r), p), p);
			}
		}
	}
```

### Factor some numbers
```C
	factors_t *factors = nut_make_Factors_w(NUT_MAX_PRIMES_64);
	for(uint64_t i = 0; i < 100; ++i){
		uint64_t n = nut_u64_rand(2, 1ull << 30);
		//nut_u64_factor_heuristic will check if the input and intermediate factors are prime
		//it returns the input n divided by all factors found, ie if n is factored completely it returns 1,
		//and it should not return > 1 unless n is greater than all the *_max fields in factor_conf and composite.
		//the primes and num_primes arguments here are null because primes for trial division must be generated by
		//some other function
		if(nut_u64_factor_heuristic(n, 0, NULL, &nut_default_factor_conf, &factors) != 1){
			printf("Failed to factor %"PRIu64"!\n", n);
		}else{
			printf("%"PRIu64" = ", n);
			nut_Factor_fprint(stdout, factors);
			printf("\n");
		}
	}
	free(factors);
```

### Factor a polynomial mod some primes
```C
	nut_Poly f[1];
	nut_Poly_init(f, 9);
	const char *end = "";
	int64_t modulus;
	if(!nut_Poly_parse(f, &modulus, "x^8 + x^7 - x^5 - x^4 - x^3 + x + 1", &end) || *end){
		printf("Failed to parse polynomial\n");
	}
	
	fprintf(stderr, "\e[1;34mComputing roots of (");
	nut_Poly_fprint(stderr, f, "x", " + ", " - ", "**", 1);
	fprintf(stderr, ") nut_i64_mod p for %"PRIu64" random primes...\e[0m\n", trials);
	
	nut_Roots roots[1];
	nut_Roots_init(roots, 8);
	
	for(uint64_t i = 0; i < 50; ++i){
		int64_t p = nut_u64_rand(2, 1ull << 30);
		while(!nut_u64_is_prime_dmr(p)){
			++p;
		}
		if(!nut_Poly_roots_modn_tmptmp(f, p, roots)){
			printf("Failed to factor (");
			nut_Poly_fprint(stdout, f, "x", " + ", " - ", "**", 1);
			printf(") nut_i64_mod %"PRId64"\n", p)
			continue;
		}
		printf("Found %"PRIu64" roots of (", roots->len);
		nut_Poly_fprint(stdout, f, "x", " + ", " - ", "**", 1);
		printf(") nut_i64_mod %"PRId64, p);
		if(roots->len){
			printf(": ");
		}else{
			printf("\n");
		}
		for(uint64_t i = 0; i < roots->len; ++i){
			printf();
			if(i != roots->len - 1){
				printf("%"PRId64", ", roots->roots[i]);
			}else{
				printf("and %"PRId64"\n", roots->roots[i]);
			}
		}
	}
	nut_Roots_destroy(roots);
	nut_Poly_destroy(f);
```

