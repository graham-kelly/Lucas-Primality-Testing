#include <gmp.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"

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

/*			compute and return gcd(x,y)
Parameters:
	x		(int) integer to to check if it is square
Returns:
	bool	x = k*k for some integer k
	
For use with small integers!
Runtime:	O(1)
*/
_Bool isSquare(int x) {
	float root_x = sqrt(x);
	if (((int) root_x) == root_x && x >= 0) {
		return 1;
	}
	return 0;
}


















