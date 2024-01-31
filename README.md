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

This library uses [Waf](https://waf.io) as the build system.
To build the release variant, simply run `./waf configure build_release` in the project root directory.
The resulting files will be in `build/release` and `build/release/bin`,
or you can run `sudo ./waf install_release` after building.

The general usage of Waf is `./waf subcommand_1 subcommand_2 ...`, where each subcommand is either user defined/modified
or builtin, and many subcommands can be specified.
Additionally, environment variables can be changed as needed normally by doing `CC=clang-13 ./waf ...` or similar.

The main subcommands to be aware of are
- `configure`: set up Waf, must be run before most other commands
- `build_<variant>`: build the indicated variant.  Results are stored in `build/<variant>`.  The currently extant variants are:
  - `coverage`: debug variant, includes address sanitizer, some ub sanitizers, and coverage collection.  builds using `gcc`
  - `debug`: debug variant, includes memory sanitizer and most ub sanitizers.  builds using `clang`
  - `release`: release variant, no sanitizers.  builds using `gcc` with `-O3`
  - `valgrind`: debug variant, no sanitizers.  builds using `clang` with `-O1`
  - `windows`: release variant, no sanitizers.  builds using `x86_64-w64-mingw32-gcc`
- `clean_<variant>`: delete build files
- `install_<variant>`: copy headers, built libraries, and build binaries (but not test binaries) to the system.  Typically requires `sudo`.  Consider running without `sudo` first to ensure it's doing what you want.
- `uninstall_<variant>`: remove installed files.  Typically requires `sudo`.  Consider running without `sudo` first to ensure it's doing what you want.
- `test_<variant>`: run tests (invoke `tests/test.py`).  This is a special rule, it will always run tests, unlike `build_*` which only builds outputs if their inputs have changed.  Requires `build_<variant>` to have been run successfully first.
  - For the `coverage` variant, `geninfo` and `genhtml` are automatically invoked after tests succeed.
- `distclean`: delete the whole build directory (build files for all variants)
- `docs`: run `doxygen` to create the documentation.  Similar to `test_<variant>` this always executes, but it is not tied to a variant.

To introduce new variants, look in `wscript` and add an appropriate section in `configure` as well as the `for` loop at the end of the file.

Waf depends on two main files: `waf`, the compressed python script comprising Waf itself, and `wscript`, the main user config.
These are somewhat analagous to the systemwide `make` executable and the `Makefile`.
Waf is intended to be distributed with the project, with a separate copy in each project.

Executables, static libraries, and test executables are automatically created based on the contents of `src`:
every C file or directory in `src/bin` is turned into its own executable in `$build/<variant>/bin`, every C file
or directory in `src/lib` is turned into its own static library in `build/<variant>`, and every C file or
directory in `src/test` is turned into its own test executable in `build/<variant>`.

`./waf distclean` will not remove `cov` and `docs`, the directories generated by `test_coverage` and `docs` respectively;
these can be removed with `rm -rf cov docs`.

Only Linux is properly supported.  To build, only (`gcc` and `ar` or `clang`) and `python` are strictly required, but `lcov`
and `doxygen` are required for coverage and documentation, `valgrind` is required for testing, and it's recommended to have both `gcc` and `clang`.
Cross compiling for windows requires `x86_64-w64-mingw32-gcc` and `libwinpthread` and isn't fully supported.  WSL or similar may be better for windows users.

## Linking
Once the library has been built, it can be linked with C programs by adding the flags `-Lbuild/release -lnut` (
or `-Lbuild/<variant> -lnut` for a different variant).
You can also install it to your system libraries by doing `sudo ./waf install_release`.
This will place the libraries in `/usr/local/lib/`, the headers in `/usr/local/include/nut/`, and non-test binaries in `/usr/local/bin/`.
Then the `-L` flag can be omitted and the library can be linked with simply `-lnut`.
Having the release variant installed systemwide and linking manually against the debug variant if some need arises is probably most convenient.

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

