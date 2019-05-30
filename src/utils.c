#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "RWG2.h"
#include "RWG7.h"
#include "ptestio.h"

/*			multiplication for numbers a + b * sqrt(c) where x_c = y_c
Parameters:
	x	an array of 3 mpz_t representing the number x[0] + x[1]*sqrt(x[2])
	y	an array of 3 mpz_t with the same third element as x representing the same number format
	p	modulus for arithmetic
Returns:
	z	= x*y (mod p)
Note that this method doesn't check that x_c = y_c, unexpected results may occur if this is not the case
Runtime:	M(p)
*/
void mpzrz_mul (mpz_t z[], mpz_t x[], mpz_t y[], mpz_t p) {
	mpz_t tmp_val[3]; mpz_init(tmp_val[0]); mpz_init(tmp_val[1]); mpz_init(tmp_val[2]);
	mpz_t result; mpz_init(result);
	mpz_mul (tmp_val[0], x[0], y[0]);			// x_a * y_a
	mpz_mul (tmp_val[1], x[1], y[1]);			// x_b * y_b
	mpz_mul (tmp_val[1], tmp_val[1], x[2]);		// x_b * y_b * x/y_c
	mpz_add (tmp_val[2], tmp_val[0], tmp_val[1]);
	mpz_mul (tmp_val[0], x[0], y[1]);			// x_a * y_b
	mpz_mul (tmp_val[1], x[1], y[0]);			// x_b * y_a
	
	mpz_set (z[0], tmp_val[2]);
	mpz_add (z[1], tmp_val[0], tmp_val[1]);
	mpz_set (z[2], x[2]);						// z = (x_a * y_a + x_b * y_b * x/y_c) + (x_a * y_b + x_b * y_a) * sqrt(x/y_c)
	
	mpz_clear(tmp_val[0]);
	mpz_clear(tmp_val[1]);
	mpz_clear(tmp_val[2]);
	
	mpz_mod (z[0], z[0], p);
	mpz_mod (z[1], z[1], p);
}

/*			squaring for numbers a + b * sqrt(c)
Parameters:
	x	an array of 3 mpz_t representing the number x[0] + x[1]*sqrt(x[2])
	p	modulus for arithmetic
Returns:
	z	= x*x (mod p)
Runtime:	M(p)
*/
void mpzrz_sqr (mpz_t z[], mpz_t x[], mpz_t p) {
	mpz_t tmp_val[2]; mpz_init(tmp_val[0]); mpz_init(tmp_val[1]);
	
	mpz_mul (tmp_val[0], x[0], x[0]);				// = x_a^2
	mpz_mul (tmp_val[1], x[1], x[1]);				// = x_b^2
	mpz_mul (tmp_val[1], tmp_val[1], x[2]);			// = x_b^2 * x_c
	mpz_add (tmp_val[1], tmp_val[0], tmp_val[1]);
	
	mpz_mul (tmp_val[0], x[0], x[1]);				// = x_a * x_b
	
	mpz_set (z[0], tmp_val[1]);
	mpz_mul_ui (z[1], tmp_val[0], 2);
	mpz_set (z[2], x[2]);							// z = (x_a^2 + x_b^2 * x_c) + (2 * x_a * x_b) * sqrt(x_c)
	
	mpz_clear(tmp_val[0]);
	mpz_clear(tmp_val[1]);
	
	mpz_mod (z[0], z[0], p);
	mpz_mod (z[1], z[1], p);
}

/*			binary exponentiation for numbers a + b * sqrt(c)
Parameters:
	base	an array of 3 mpz_t representing the number x[0] + x[1]*sqrt(x[2])
	exp		mpz_t exponent to raise base to
	p		modulus for arithmetic
Returns:
	int_result		= integer part of base^exp (mod p)
	root_result		= non-integer part of base^exp (mod p) (may not be zero)
Runtime:	log(exp)M(p)
*/
void mpzrz_exp (mpz_t int_result, mpz_t root_result, mpz_t base[3], mpz_t exp, mpz_t p) {
	mpz_t cexp; mpz_init(cexp);
	mpz_t term[3]; mpz_init(term[0]); mpz_init(term[1]); mpz_init(term[2]);
	mpz_t product[3]; mpz_init(product[0]); mpz_init(product[1]); mpz_init(product[2]);
	mpz_set (cexp, exp);									// need a copy of the exponent to work with
	mpz_set (term[0], base[0]);								// initialize term
	mpz_set (term[1], base[1]);
	mpz_set (term[2], base[2]);
	while (! mpz_tstbit (cexp, 0)) {						// before first set bit of exponent
		mpzrz_sqr (term, term, p);							// square the term to multiply by
		mpz_tdiv_q_2exp (cexp, cexp, 1);					// right shift exponent by 1
	}														// now at first set bit of exponent => can initialize product
	mpz_set (product[0], term[0]);							// initialize product
	mpz_set (product[1], term[1]);
	mpz_set (product[2], term[2]);
	mpz_tdiv_q_2exp (cexp, cexp, 1);						// right shift exponent by 1
	while (mpz_cmp_ui (cexp, 0) != 0) {						// while exponent has more bits
		mpzrz_sqr (term, term, p);							// square the term to multiply by
		if (mpz_tstbit (cexp, 0)) {							// if rightmost bit is 1
			mpzrz_mul (product, product, term, p);			// do multiplication
		}
		mpz_tdiv_q_2exp (cexp, cexp, 1);					// right shift exponent by 1
	}
	mpz_set (int_result, product[0]);
	mpz_set (root_result, product[1]);
	mpz_clear(cexp);
	mpz_clear(term[0]); mpz_clear(term[1]); mpz_clear(term[2]);
	mpz_clear(product[0]); mpz_clear(product[1]); mpz_clear(product[2]);
}

/*			Trial division by the first million prime numbers
Parameters:
	N	the number to check primality of
	ndiv		the number of prime numbers to try dividing N by (0 => do all 1,000,000)

Returns:
	integer		-1 if failed to open "primessmall.txt" (this file should exist and contain the first million primes on separate lines)
				0 if no prime divisor of N was found (doesn't guarentee primality!)
				>1 => the returned value divides N
Runtime:	M(N) for division
*/
int trial_div (mpz_t N, int ndiv) {
	if (ndiv > 1000000 || ndiv < 1) {
		ndiv = 1000000;
	}
	FILE *f = fopen ("primessmall.txt", "r");
	if (f == NULL) {printf("Error opening file.\n"); return 0;}
	int prime;
	fscanf (f, "%d", &prime);
	while (!feof(f)) {
		if (mpz_divisible_ui_p (N, prime)) {				// if N (mod p) = 0
			fclose(f);
			return prime;
		}
		fscanf (f, "%d", &prime);
	}
	fclose(f);
	return 0;
}

/*	Trial Division to catch small primes to add to files
Parameters:
	integers values as given in RWG section 2 primality tests (2.8 / 2.10)
	the following relationship between the inputs should hold: n < (log(15485863-y)-log(A))/log(r)
Returns:
	1 => prime
	0 => not prime
*/
int trail_div_small(int A, int r, int n, int y) {
	mpz_t N; mpz_init (N);
	mpz_ui_pow_ui (N, r, n);
	mpz_mul_ui (N, N, A);
	if (y == 1) {
		mpz_add_ui (N, N, 1);
	}
	else {
		mpz_sub_ui (N, N, 1);
	}						// N = A*r^n+y
	int divRes = trial_div(N, 0);
	if (divRes) {			// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
		if (mpz_cmpabs_ui (N, divRes) == 0) {
			mpz_clear(N);
			return 1;
		}
	}
	mpz_clear(N);
	return 0;
}

/*			compute and return gcd(x,y)
Parameters:
	x,y		integers to find the gcd of
Returns:
	integer		gcd(x,y)
	
For use with small integers! (large integers should instead use the mpz version)
Runtime:	O(1)
	only called with -100 < x, y < 100
*/
int gcd(int x, int y) {
	int temp;
	while (y != 0) {
		temp = x % y;
		x = y;
		y = temp;
	}
	return x;
}

_Bool isSquare(int x) {
	float root_x = sqrt(x);
	if (((int) root_x) == root_x && x >= 0) {
		return 1;
	}
	return 0;
}

//	Comparison of two integers for use with qsort (equivalent to <)
int cmp (const void * a, const void * b) {
	return (*(int *)a - *(int *)b);
}

/*	Sort integer arrays and remove duplicates
Parameters:
	arr			the array to be sorted
	init_size	initial size of integer array
Returns:
	arr is now sorted and contains no duplicates (and has minimal size)
	new_size of arr is returned
Runtime:	O(init_size)
*/
int remove_dup_array(int * arr, int init_size) {
	qsort(arr, init_size, sizeof(int), cmp);
	int current = *arr;
	int new_size = 1;
	for (int i = 1; i < init_size; i++) {
		if (current != *(arr + i)) {
			current = *(arr + i);
			*(arr + new_size) = current;
			new_size++;
		}
	}
	arr = (int *) realloc (arr, sizeof(int) * new_size);
	return new_size;
}

/*	Get information about which groups of n will always give composite numbers for N=Ar^n+y
Parameters:
	A,r,n,y		are positive integers as given in RWG section 2 primality tests (theorems 2.8 / 2.10)
Returns:
	not_to_run			contains triples of p, x, z such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
	size_not_to_run		
Runtime:	O(p_n * M(A * r^(p_n))) 
	where p_n is the nth prime number

Note that initial size of not_to_run should be 2*n
*/
int get_not_to_run_2_8_10 (int * not_to_run, mpz_t A, int r, int n, int y) {
	FILE *f = fopen ("primessmall.txt", "r");
	if (f == NULL) {printf("Error opening file.\n"); return -1;}
	int size_not_to_run = 0;
	int p;
	mpz_t test_N; mpz_init (test_N);
	mpz_t rEXPj; mpz_init (rEXPj);
	_Bool found_x;
	_Bool found_z;
	for (int i = 0; i < n; i++) {
		fscanf (f, "%d", &p);														// get prime from primesSmall
		if (p >= n) {
			break;
		}
		mpz_tdiv_r_ui (test_N, A, p);
		if (y == 1) {
			mpz_add_ui (test_N, test_N, 1);
		}
		else {
			mpz_sub_ui (test_N, test_N, 1);
		}						// test_N = A+y
		mpz_set_ui (rEXPj, 1);
		found_x = 0;
		found_z = 0;
		for (int j = 1; j < p; j++) {					// find x, z such that A*r^x+z = 0 (mod p) and r^z = 1 (mod p)
			if (y == 1) {
				mpz_sub_ui (test_N, test_N, 1);
				mpz_mul_ui (test_N, test_N, r);
				mpz_add_ui (test_N, test_N, 1);
			}
			else {
				mpz_add_ui (test_N, test_N, 1);
				mpz_mul_ui (test_N, test_N, r);
				mpz_sub_ui (test_N, test_N, 1);
			}						// test_N = A*(r^j) + y
			if (!found_x && mpz_divisible_ui_p (test_N, p)) {
				*(not_to_run + 2 * size_not_to_run) = j;
				found_x = 1;
			}
			mpz_mul_ui (rEXPj, rEXPj, r);
			if (!found_z && mpz_congruent_ui_p(rEXPj, 1, p)) {
				*(not_to_run + 2 * size_not_to_run + 1) = j;
				found_z = 1;
			}
			if (found_x && found_z) {					// not_to_run[size] = (x, z) such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
				size_not_to_run++;
				break;
			}
		}
	}
	mpz_clear (test_N);
	mpz_clear (rEXPj);
	fclose(f);
	not_to_run = (int *) realloc (not_to_run, sizeof(int) * 2 * size_not_to_run);
	return size_not_to_run;
}

int get_n_to_run_2_8_10 (int * to_run, mpz_t A, int r, int n, int nf, int y) {
	int *not_to_run = (int *) malloc(sizeof(int) * 2 * nf);						// get p, x, z such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
	int size_not_to_run = get_not_to_run_2_8_10 (not_to_run, A, r, nf, y);	// no n = x+kz for any integer k should be tested (all divisble by p)
	int k = 0;
	int exp;
	for (int i = 0; i < size_not_to_run; i++) {		// for each small prime tested
		exp = (int) ceil ((n - *(not_to_run + 2 * i)) / (double) *(not_to_run + 2 * i + 1)) * (*(not_to_run + 2 * i + 1)) + *(not_to_run + 2 * i) - n;		// get initial exponent's index (x - n + z*ceil( (n-x)/z) )
		while (exp < nf - n + 1) {
			*(to_run + exp) = 1;					// this N with exponent (exp+n) will be divisible by some small prime
			exp += *(not_to_run + 2 * i + 1);		// add z (next iteration remove new exponent from the list if still in range)
		}
	}
	int size_to_run = 0;
	for (int i = 0; i < nf - n + 1; i++) {
		if (!*(to_run + i)) {						// has not been marked as divisible by some x + kz above
			*(to_run + size_to_run) = i + n;		// this integer should be run
			size_to_run++;
		}
	}
	to_run = (int *) realloc (to_run, sizeof(int) * size_to_run);
	return size_to_run;
}

int LARGEST_IN_SMALL_PRIMES = 15485863;
/*
Parameters:
	fileStr			file to output results to
	test_num		2 => do tests from section 2 of RWG, 7 => tests from section 7
	A, r, y_eta		integers as described in RWG Theorems 2.8, 2.10, 7.2, 7.4, 7.5
	n, nf			test for primes using above constants and exponent n, n+1, n+2, ... , nf-2, nf-1, nf
*/
void OEIS_testing (char * fileStr, int test_num, int A, int r, int n, int nf, int y_eta) {
	if (test_num == 2) {															// do the tests from section 2
		mpz_t N; mpz_init (N);
		while (n < (log(LARGEST_IN_SMALL_PRIMES - y_eta)-log(A))/log(r)) {				// if the potential prime is <= this big it will be in smallPrimes.txt
			mpz_ui_pow_ui (N, r, n);
//gmp_printf("here");
			mpz_mul_ui (N, N, A);
//gmp_printf("here");
			if (y_eta == 1) {
				mpz_add_ui (N, N, 1);
			}
			else {
				mpz_sub_ui (N, N, 1);
			}						// N = A*r^n+y
			int divRes = trial_div(N, 0);
			if (divRes) {			// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
				if (mpz_cmpabs_ui (N, divRes) == 0) {
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file_OEIS(n, fileStr);			// add to file
				}
			}
			printf("%d,", n);
			fflush(stdout);
			n++;
		}
		mpz_clear(N);
		mpz_t Az; mpz_init (Az);
		mpz_set_ui (Az, A);
		int *to_run = (int *) calloc((nf - n + 1), sizeof(int));
		int size_to_run = get_n_to_run_2_8_10 (to_run, Az, r, n, nf, y_eta);
		if (r == 2) {
			for (int i = 0; i < size_to_run; i++) {
				if (primality_test_2_10 (Az, to_run[i], y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", to_run[i]);
					add_prime_file_OEIS(to_run[i], fileStr);			// add to file
				}
				printf("%d,", to_run[i]);
				fflush(stdout);
			}
		}
		else {
			for (int i = 0; i < size_to_run; i++) {
				if (primality_test_2_8 (Az, r, to_run[i], y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", to_run[i]);
					add_prime_file_OEIS(to_run[i], fileStr);			// add to file
				}
				printf("%d,", to_run[i]);
				fflush(stdout);
			}
		}
		mpz_clear (Az);
	}
	else {																// use the tests from section 7
		if (r == 2 || (r%4 != 1)) {
			printf ("Section 7 primality tests do not apply to r = 2 or r = 1 (mod 4)\n");
			return;
		}
		if (r == 5) {
			while (n <= nf) {
				if (primality_test_7_5 (A, n, y_eta) == 1) {			// specific to r = 5
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file_OEIS(n, fileStr);			// add to file						}
				}
				else {
					printf("%d,", n);
				}
				fflush(stdout);
				n++;
			}
		}
		else {			// use primality testing from both Theorem 7.2 and 7.4 as they require computing the same sequences are have the same requirements
			while (n <= nf) {
				if (primality_test_7_2_4 (A, r, n, y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file_OEIS(n, fileStr);			// add to file
				}
				else {
					printf("%d,", n);
				}
				fflush(stdout);
				n++;
			}
		}
	}
	printf("\n");
	return;
}

void SIDH_testing (char * fileStr, int bit_length, int r) {
	mpz_t A; mpz_init (A);
	for (int x = (bit_length / 2) - 20; x < bit_length / 2; x++) {
		int y;
		for (y = (int) floor((bit_length / 2) / log2(r)) - 20; y < (x-1)/log2(r) + 2; y++) {
			if (y < 0) {
				continue;
			}
			for (int f = 1; f < 100; f++) {
				if (f%r == 0) {
					break;
				}
				mpz_set_ui (A, f);
				mpz_mul_2exp (A, A, x);
				if (primality_test_2_8 (A, r, y, -1) == 1) {
					printf("\n\nPrime for values: f = %d, x = %d, y = %d\n\n", f, x, y);
					add_prime_file_SIDH(f, x, y, fileStr);				// add to file
				}
				printf("%d/%d/%d,", f, x, y);
				fflush(stdout);
			}
		}
		for ( ; y < bit_length / 2; y++) {
			for (int f = 1; f < 100; f+=2) {
				mpz_ui_pow_ui (A, r, y);
				mpz_mul_ui (A, A, f);
				if (primality_test_2_10 (A, x, -1) == 1) {
					printf("\n\nPrime for values: f = %d, x = %d, y = %d\n\n", f, x, y);
					add_prime_file_SIDH(f, x, y, fileStr);				// add to file
				}
				printf("%d/%d/%d,", f, x, y);
				fflush(stdout);
			}
		}
	}
	// select A = f * 2^x, r^n such that r^n has about x bits (>x => use 2.8 else use 2.10)
	// can use small f, want 
	mpz_clear (A);
	printf("\n");
	return;
}


















