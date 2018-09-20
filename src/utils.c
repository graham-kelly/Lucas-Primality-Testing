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
*/
void mpzrz_mul (mpz_t z[3], mpz_t x[3], mpz_t y[3], mpz_t p) {
	mpz_t constant1; mpz_init(constant1);
	mpz_t constant2; mpz_init(constant2);
	mpz_t result; mpz_init(result);
	mpz_mul (constant1, x[0], y[0]);			// x_a * y_a
	mpz_mul (constant2, x[1], y[1]);			// x_b * y_b
	mpz_mul (constant2, constant2, x[2]);		// x_b * y_b * x/y_c
	mpz_add (result, constant1, constant2);
	mpz_mul (constant1, x[0], y[1]);			// x_a * y_b
	mpz_mul (constant2, x[1], y[0]);			// x_b * y_a
	
	mpz_set (z[0], result);
	mpz_add (z[1], constant1, constant2);
	mpz_set (z[2], x[2]);						// z = (x_a * y_a + x_b * y_b * x/y_c) + (x_a * y_b + x_b * y_a) * sqrt(x/y_c)
	
	mpz_clear(constant1);
	mpz_clear(constant2);
	mpz_clear(result);
	
	mpz_mod (z[0], z[0], p);
	mpz_mod (z[1], z[1], p);
}

/*			squaring for numbers a + b * sqrt(c)
Parameters:
	x	an array of 3 mpz_t representing the number x[0] + x[1]*sqrt(x[2])
	p	modulus for arithmetic
Returns:
	z	= x*x (mod p)
*/
void mpzrz_sqr (mpz_t z[3], mpz_t x[3], mpz_t p) {
	mpz_t constant; mpz_init(constant);
	mpz_t result; mpz_init(result);
	
	mpz_mul (constant, x[0], x[0]);				// = x_a^2
	mpz_mul (result, x[1], x[1]);				// = x_b^2
	mpz_mul (result, result, x[2]);				// = x_b^2 * x_c
	mpz_add (result, constant, result);
	
	mpz_mul (constant, x[0], x[1]);				// = x_a * x_b
	
	mpz_set (z[0], result);
	mpz_mul_ui (z[1], constant, 2);
	mpz_set (z[2], x[2]);						// z = (x_a^2 + x_b^2 * x_c) + (2 * x_a * x_b) * sqrt(x_c)
	
	mpz_clear(constant);
	mpz_clear(result);
	
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
*/
void mpzrz_exp (mpz_t int_result, mpz_t root_result, mpz_t base[3], mpz_t exp, mpz_t mod) {
	mpz_t cexp; mpz_init(cexp);
	mpz_t term[3]; mpz_init(term[0]); mpz_init(term[1]); mpz_init(term[2]);
	mpz_t product[3]; mpz_init(product[0]); mpz_init(product[1]); mpz_init(product[2]);
	mpz_set (cexp, exp);									// need a copy of the exponent to work with
	mpz_set (term[0], base[0]);								// initialize term
	mpz_set (term[1], base[1]);
	mpz_set (term[2], base[2]);
	while (! mpz_tstbit (cexp, 0)) {						// before first set bit of exponent
		mpzrz_sqr (term, term, mod);						// square the term to multiply by
		mpz_tdiv_q_2exp (cexp, cexp, 1);					// right shift exponent by 1
	}														// now at first set bit of exponent => can initialize product
	mpz_set (product[0], term[0]);							// initialize product
	mpz_set (product[1], term[1]);
	mpz_set (product[2], term[2]);
	mpz_tdiv_q_2exp (cexp, cexp, 1);						// right shift exponent by 1
	while (mpz_cmp_ui (cexp, 0) != 0) {						// while exponent has more bits
		mpzrz_sqr (term, term, mod);						// square the term to multiply by
		if (mpz_tstbit (cexp, 0)) {							// if rightmost bit is 1
			mpzrz_mul (product, product, term, mod);		// do multiplication
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
*/
int trial_div (mpz_t N, int ndiv) {
	if (ndiv > 1000000 || ndiv < 1) {
		ndiv = 1000000;
	}
	FILE *f = fopen ("primessmall.txt", "r");
	if (f == NULL) {printf("Error opening primessmall.txt.\n"); return -1;}
	int i = 0;
	int prime;
	int rprime = 0;
	int character;
	mpz_t r; mpz_init(r);
	while (i++ < ndiv) {
		prime = 0;
		while ((character = getc(f)) != EOF) {
			if (character == 13 || character == 10) {
				break;
			}
			prime *= 10;
			prime += character - 48;
		}
		if (! mpz_tdiv_r_ui(r, N, prime) && !mpz_cmp_ui(N, prime) == 0) {				// if N (mod p) = 0
			rprime = prime;
			break;
		}
		//if (mpz_cmp_ui (N, prime * prime) < 0) {break;}				// doesn't occur for reasonably large primes so no need to test
	}
	mpz_clear(r);
	fclose(f);
	return rprime;
}

/*			compute and return gcd(x,y)
Parameters:
	x,y		integers to find the gcd of
Returns:
	integer		gcd(x,y)
	
For use with small integers
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

/*			check if n is a square integer
Parameters:
	n		integer to check if there exists an integer a such that a^2 = n
Returns:
	_Bool	true or false
*/
_Bool isSquare(int n) {
	return !(n - (int) pow(floor(sqrt(n)), 2));
}