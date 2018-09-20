#include <gmp.h>
#include <math.h>
#include "utils.h"
#include "RWG2.h"
#include "RWG2S.h"


/*			Get the coefficients of the g_r(x) polynomial from RWG Theorem 2.4
Parameters:
	r	the small (odd) prime used to compute N
	d	delta value based on p and q (which were selected based on N)
	N	N = A*r^n+gamma which is hoped to be prime

Returns:
	constants	will contain the coefficients of g_r(x)
	bool		if inversion fails (mod N) then N cannot be prime and 0 will be returned, otherwise returns 1 (doesn't guarentee N is prime)

takes an array of mpz_t elements and fills them with the coefficients of g_r(x)
g_r(x) = sum(j = 0...k) [ (r/2j+1) * binomial_coef(k+j, k-j) * delta^j * x^(2j+1) ]
where k = (r-1)/2 computing the k+1 coefficients and storing them in constants
*/
_Bool get_g_r_constants (mpz_t constants[], long int r, int d, mpz_t N) {
	int k = (r-1)/2;											// note r is an odd prime (this is exact)
	mpz_t term; mpz_init(term);
	mpz_t constant; mpz_init(constant);
	_Bool prime = 1;											// true if inversion (mod N) has not failed
	int j = 0;
	while (j < k && prime) {
		mpz_ui_pow_ui (term, d, j);								// d^j
		mpz_mul_ui (term, term, r);								// r * d^j
		mpz_set_ui (constant, k - j);
		prime = mpz_invert(constant, constant, N);				// bool return wll be non-zero if inversion was sucessful every time
		mpz_mul (term, term, constant);							// = r/(k-j) * d^j
		mpz_bin_uiui (constant, k + j, k - j - 1);				
		mpz_mul (term, term, constant);							// = r/(2j+1) binom(k+j, k-j) delta^j
		mpz_mod (term, term, N);								// gives term (mod N)
		mpz_set (constants[j], term);
		j++;
	}
	mpz_ui_pow_ui (constants[j], d, j);							// d^j
	mpz_mod (constants[j], constants[j], N);					// gives term (mod N)
	mpz_clear (term);
	mpz_clear (constant);
	return prime;
}

/*			Get appropreate values of p and q for primality tests from RWG Theorem 2.8/2.10
Parameters:
	N	N=Ar^n+gamma
	y	gamma
	
Returns:
	arr[3]	integer array containing p, q and delta values
	
pass arr[p, q, d] as parameter; this array will be modified to contain correct values with respect to N, y (other parameters)
all array entries will be set to 0 if not successful

p, q coprime, delta = p^2 - 4q, gcd (N, q) == 1 and Jacobi(delta,N) == gamma will be true after this method (required for primality tests in RWG section 2)
*/
void old_find_p_q (int arr[3], mpz_t N, short y) {
	mpz_t constant; mpz_init (constant);
	mpz_t constant2; mpz_init (constant2);
	int p = 1, q = 1, d;												// initial values
	int max_val = 100;													// set max absolute value for p and q (method will fail if no such values are found)
	_Bool notFound = 1;
	while ((notFound) && p <= max_val) {
		d = p*p - 4 * q;
		while (q<=max_val && (notFound) && d>0) {						// q = 1, -1, 2, -2, 3, -3, ...
			mpz_set_si (constant, d);
			mpz_set_si (constant2, q);
			if (y == mpz_jacobi (constant, N) && 1 == mpz_jacobi(constant2, N)) {						// gamma = jacobi(d,N) this doesn't make sense...
				notFound = 0;
			}
			else {
				if (q > 0) {
					q = -q;												// q = 1, -1, 2, -2, 3, -3, ...
				}
				else {
					q = -q + 1;											// q = 1, -1, 2, -2, 3, -3, ...
				}
				d = p*p - 4 * q;
			}
		}
		if (notFound) {													// try next value of p
			q = 1;
			p = p + 1;
		}
	}
	if (notFound) {
		if (p == max_val) {
			gmp_printf("WARNING! p and q were NOT FOUND for N = %Zd!\n", N);
		}
		p = 0;
		q = 0;
		d = 0;
	}
	arr[0] = p;															//set up return values
	arr[1] = q;
	arr[2] = d;
	mpz_clear(constant);
	mpz_clear(constant2);
}

void find_p_q (int arr[3], mpz_t N, short y) {
	mpz_t constant; mpz_init (constant);
	mpz_t constant2; mpz_init (constant2);
	int p = 1, q = 0, d;												// initial values
	int max_val = 100;													// set max absolute value for p and q (method will fail if no such values are found)
	_Bool notFound = 1;
	while (notFound && p <= max_val) {
		while (notFound && q<=max_val) {
			if (q > 0) {
				q = -q;												// q = 1, -1, 2, -2, 3, -3, ...
			}
			else {
				q = -q + 1;											// q = 1, -1, 2, -2, 3, -3, ...
			}
			d = p*p - 4 * q;
			while (isSquare(d)) {
				if (q > 0) {
					q = -q;												// q = 1, -1, 2, -2, 3, -3, ...
				}
				else {
					q = -q + 1;											// q = 1, -1, 2, -2, 3, -3, ...
				}
				d = p*p - 4 * q;
			}
			mpz_set_si (constant, d);
			mpz_set_si (constant2, q);
			if (y == mpz_jacobi (constant, N) && 1 == mpz_jacobi(constant2, N)) {						// gamma = jacobi(q,N) (typo in paper?)
				notFound = 0;
			}
		}
		if (notFound) {													// try next value of p
			q = 1;
			p = p + 1;
		}
	}
	if (notFound) {
		if (p == max_val) {
			gmp_printf("WARNING! p and q were NOT FOUND for N = %Zd!\n", N);
		}
		p = 0;
		q = 0;
		d = 0;
	}
	arr[0] = p;															//set up return values
	arr[1] = q;
	arr[2] = d;
	mpz_clear(constant);
	mpz_clear(constant2);
}

/*			Theorem 2.8 Primality Test from RWG
Parameters:
	A	an even positive integer not divisible by r
	r	a (small) odd prime
	n	a positive integer
	y	= +/-1

Returns:
	integer type	1 => N = Ar^n+y is prime
					0 => N is composite
					-1 => the primality test doesn't apply (gamma^2 != 1 or A not even)
					2 => N is possibly prime
	prints to console result of primality testing if information other than N is composite was obtained

Let r be an odd prime and N = Ar^n + y where y^2 = 1,
2|A, A<=2r^n, (N,q) = 1 and (d/N) = y. Put s_0 = u_a/q^(A/2) (mod N)
and define s_i+1 = g_r(s_i) (mod N) (i = 0, 1, ... , n-1)
If not N|s_n then N is composite. If alpha (<=n) is the least value of 
i such that not N|s_alpha-1 and N|s_alpha then a prime devisor of N
must have the form 2kr^alpha + y. Also if A<= 2r^(2alpha-n)
then N is prime
*/
int primality_test_2_8 (long int A, long int r, long int n, short y) {
	int result = 1;													//assume result will be prime (if the number is composite then result will be changed)
	mpz_t N; mpz_init(N);
	mpz_t constant; mpz_init(constant);
	mpz_ui_pow_ui (constant, r, n);
	mpz_mul_ui (N, constant, A);
	if (y == 1) {
		mpz_add_ui (N, N, 1);										// N = A * r^n + y
	}
	else {
		mpz_sub_ui (N, N, 1);										// N = A * r^n + y
	}
	if (mpz_cmp_ui (constant, A/2) < 0) {							// if A > 2r^n
		result = -1;												// test doesn't apply
	}
	if (y*y != 1 || A%2 != 0) {										// if gamma^2 != 1 or A is odd
		result = -1;
	}
	if (result == 1 && trial_div(N, 0)) {							// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
		result = 0;
	}
	long int a = -1;
/*	if (log2(n) > r * log2(r)) {									// this commented out section finds primality as well as possible prime divisors if N is not found to be prime
		long int k = (r-1)/2;										// the k for calculating g_r(x)
		mpz_t g_r_constants[k + 1];									// has indexes 0 - k available
		while (k >= 0) {
			mpz_init (g_r_constants[k]);							// initialize constants for calculating g_r(x)
			k--;
		}
		int d = get_s0 (constant, A, y, N);							// get s_0 and store in constant (and the associated delta value used)
		get_g_r_constants (g_r_constants, r, d, N);
		if (d == 0 && result == 1) {
			result = 0;												// integer inversion failed so N is not prime
		}
		long int i = 0;
		mpz_t constants[2]; mpz_init (constants[0]); mpz_init (constants[1]);	// used in each iteration of get_next_s_i method (don't want to initialize n times)
		while ((0 != mpz_cmp_ui (constant, 0)) && i++ < n) {		// while s_i != 0
			get_next_s_i (constant, r, constants, g_r_constants, N);			// finds s_i+1 and stores it in s_i
		}
		mpz_clear (constants[0]); mpz_clear (constants[1]);
		if (0 == mpz_cmp_ui (constant, 0) && i <= n) {
			a = i;													// first value of i such that N !| s_i-1 and N | s_i (i has already been incremented)
		}
		k = (r-1)/2;												// the k for calculating g_r(x)
		while (k >= 0) {
			mpz_clear (g_r_constants[k]);							// clear constants for calculating g_r(x)
			k--;
		}
	}
	else {
		int pqd[3];
		find_p_q (pqd, N, y);
		a = bin_search_s_i (A, r, n, pqd, N);						// note that binary search is not necessarily faster *** CAN PROBABLY BE OPTIMIZED USING SOME CUTOFF ***
	}
	if (a == -1 && result == 1) {									// if N !| s_a (s_n if no suitable a could be found)
		result = 0;													// N is not prime
	}
	if (result == 1) {
		if (2 * a - n < 0 || log(A/2)/log(r) >= 2*a - n) {
			result = 2;												// primality status unknown
		}
	}
*/
	if (result == 1) {												// this conditional section gives a test for primality but does not attempt to find prime divisors if this is not the case
		long int min_alpha = (long int) floor((log(A/2) / log(r) + n) / 2);		// minimum alpha which removes ability to confirm primality
		int pqd[3];
		find_p_q (pqd, N, y);
		if (! get_s_i (constant, A, r, min_alpha, pqd, N)) {
			result = 0;
		}
		else {
			if (mpz_cmp_ui (constant, 0) == 0) {
				result = 2;											// can't confirm primality (doesn't do compositeness test (if N !| s_n) which is the most expensive check)
			}
			else {
				if (! get_s_i (constant, A, r, n, pqd, N)) {		//compute s_n and check that it is 0 (else N is composite)
					result = 0;
				}
				else {
					if (mpz_cmp_ui (constant, 0) != 0) {
						result = 0;
					}
				}
			}
		}
	}
	
	mpz_clear (N);
	mpz_clear (constant);
	
	if (result == 1) {
		gmp_printf("A = %d, r = %d, n = %d, y = %d is prime!\n", A, r, n, y);
	}
 	if (result == 0 && a != -1) {
		gmp_printf("A = %d, r = %d, n = %d, y = %d is composite with prime divisors of the form 2kr^(%d) +/- 1.\n", A, r, n, y, a);
	}
	if (result == -1) {
		gmp_printf("The test doesn't apply to A = %d, r = %d, n = %d, y = %d.\n", A, r, n, y);
	}
	if (result == 2 && a == -1) {
		gmp_printf("A = %d, r = %d, n = %d, y = %d may be prime.\n", A, r, n, y);
	}
	if (result == 2 && a != -1) {
		gmp_printf("A = %d, r = %d, n = %d, y = %d may be prime. Prime divisors would have the form 2kr^(%d) +/- 1.\n", A, r, n, y, a);
	}
	
	return result;
}

/*			Theorem 2.10 Primality Test from RWG
Parameters:
	A	an odd positive integer
	n	a positive integer
	y	= +/-1

Returns:
	integer type	1 => N = A*2^n+y is prime
					0 => N is composite
					-1 => the primality test doesn't apply (gamma^2 != 1)
					2 => N is possibly prime
	prints to console result of primality testing if information other than N is composite was obtained

Let N = A*2^n+y, where n >=2, 2!|A, y^2 = 1, A<2^n, (d/N) = y.
Set r_0 = u_A(p,q)/q^(A/2), r_1 = v_A(p,q)/q^(A/2) (mod N) and define
r_(i+1) = r_i^2 - 2 (mod N)
Then N is not a prime if N !| r_i (i = 0, 1, ... , n). If a is the least
value of i such that N|r_i then any prime devisor of N must have the form
k*2^a +/- 1. Also if A < 2^(2a-n)-1, then N is a prime
*/
short primality_test_2_10 (long int A, long int n, short y) {
	short result = 1;												// assume the number is prime until found composite
	if (A % 2 == 0 || y*y != 1) {
		result = -1;												// test doesn't apply
	}
	mpz_t N; mpz_init(N);
	mpz_t constant; mpz_init(constant);
	mpz_ui_pow_ui (constant, 2, n);									// constant = 2^n
	if (mpz_cmp_ui (constant, A) < 0) {								// if A >= 2^n
		result = -1;												// test doesn't apply
	}
	mpz_mul_ui (N, constant, A);
	if (y == 1) {
		mpz_add_ui (N, N, 1);										// N = A * r^n + y
	}
	else {
		mpz_sub_ui (N, N, 1);										// N = A * r^n + y
	}
	if (result == 1 && trial_div(N, 0)) {							// trial division by first million primes
		result = 0;
	}
	mpz_t r_i; mpz_init (r_i);
	if (result == 1) {
		if (get_r0_r1 (constant, r_i, A, y, N) == 0) {				// constant will (temporarily) store r0, and r_i will store r1
			result = 0;												// if a multiplicative inverse doesn't exist mod N then N is not prime
		}
	}
	long int a = -1;
	if (result == 1) {
		if (0 == mpz_cmp_ui (constant, 0)) {							// if r0 = 0
			a = 0;
		}
		if (0 == mpz_cmp_ui (r_i, 0)) {									// if r1 = 0
			a = 1;
		}
		long int i = 1;
		while (i++ < n) {												// while a value hasn't been found yet
			mpz_powm_ui (r_i, r_i, 2, N);
			mpz_sub_ui (r_i, r_i, 2);									// r_(i+1) = r_i^2 - 2 (mod N)
			if (mpz_cmp_ui (r_i, 0) == 0) {								// if N | r_i
				a = i;													// a value found
				break;
			}
		}
		if (a == -1) {													// if N !| r_i for i = 0, 1, ... , n
			result = 0;													// then N is composite
		}
		if (result == 1 && (2*a <= n || log(A)/log(2) >= 2*a - n)) {
			result = 2;													// potentially prime
		}
	}
	mpz_clear (N);
	mpz_clear (constant);
	mpz_clear (r_i);
	
	if (result == 1) {
		gmp_printf("A = %d, r = 2, n = %d, y = %d is prime!\n", A, n, y);
	}
 	if (result == 0 && a != -1) {
		gmp_printf("A = %d, r = 2, n = %d, y = %d is composite with prime divisors of the form k 2^(%d) +/- 1.\n", A, n, y, a);
	}
	if (result == -1) {
		gmp_printf("The test doesn't apply to A = %d, r = 2, n = %d, y = %d.\n", A, n, y);
	}
	if (result == 2 && a == -1) {
		gmp_printf("A = %d, r = 2, n = %d, y = %d may be prime.\n", A, n, y);
	}
	if (result == 2 && a != -1) {
		gmp_printf("A = %d, r = 2, n = %d, y = %d may be prime. Prime divisors would have the form k 2^(%d) +/- 1.\n", A, n, y, a);
	}
	
	return result;
}

