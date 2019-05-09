#include <gmp.h>
#include <math.h>
#include "utils.h"
#include "RWG2.h"
#include "RWG2S.h"

#define MAX_ABSOLUTE_P_Q_SIZE 100;

/*			Get the coefficients of the g_r(x) polynomial from RWG Theorem 2.4
Parameters:
	r	the small (odd) prime used to compute N
	d	delta value based on p and q (which were selected based on N)
	N	N = A*r^n+gamma which is hoped to be prime

Returns:
	coefficients	will contain the coefficients of g_r(x)
	bool		if inversion fails (mod N) then N cannot be prime and 0 will be returned, otherwise returns 1 (doesn't guarentee N is prime)

takes an array of mpz_t elements and fills them with the coefficients of g_r(x)
g_r(x) = sum(j = 0...k) [ (r/2j+1) * binomial_coef(k+j, k-j) * delta^j * x^(2j+1) ]
where k = (r-1)/2 computing the k+1 coefficients and storing them in coefficients
Runtime:	O(r*log(N)^2)
	r inversions (mod N)
*/
_Bool get_g_r_coef (mpz_t coefficients[], long int r, int d, mpz_t N) {
	int k = (r-1)/2;												// note r is an odd prime (this is exact)
	mpz_t term; mpz_init(term);
	mpz_t tmp_val; mpz_init(tmp_val);
	int j = 0;
	while (j <= k) {
		mpz_ui_pow_ui (term, d, j);									// d^j
		if (k != j) {												// if k == j then coefficient * binomial coefficient = 1
			mpz_mul_ui (term, term, r);								// r * d^j
			mpz_set_ui (tmp_val, k - j);
			if (!mpz_invert(tmp_val, tmp_val, N)) {					// if inversion fails then N is not prime
				mpz_clear (term);
				mpz_clear (tmp_val);
				return 0;
			}
			mpz_mul (term, term, tmp_val);							// = r/(k-j) * d^j
			mpz_bin_uiui (tmp_val, k + j, k - j - 1);				
			mpz_mul (term, term, tmp_val);							// = r/(2j+1) binom(k+j, k-j) delta^j
		}
		mpz_mod (term, term, N);									// gives term (mod N)
		mpz_set (coefficients[j], term);
		j++;
	}
	mpz_clear (term);
	mpz_clear (tmp_val);
	return 1;
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
Runtime:	O(1)
	runs jacobi on int with absolute size < 100 (MAX_ABSOLUTE_P_Q_SIZE)
*/
void find_p_q (int arr[3], mpz_t N, short y) {
	mpz_t tmp_val[2]; mpz_init (tmp_val[0]); mpz_init (tmp_val[1]);
	int p = 1, q = 1, d;												// initial values
	int max_val = MAX_ABSOLUTE_P_Q_SIZE;								// set max absolute value for p and q (method will fail if no such values are found)
	while (p <= max_val) {
		d = p * p - 4 * q;
		while (q <= max_val && d > 0) {									// q = 1, -1, 2, -2, 3, -3, ...
			mpz_set_si (tmp_val[0], d);
			mpz_set_si (tmp_val[1], q);
			if (y == mpz_jacobi (tmp_val[0], N) && 1 == mpz_jacobi(tmp_val[1], N)) {						// gamma = jacobi(d,N) this doesn't make sense...
				arr[0] = p;												//set up return values
				arr[1] = q;
				arr[2] = d;
				mpz_clear(tmp_val[0]);
				mpz_clear(tmp_val[1]);
				return;
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
		q = 1;															// try next value of p
		p = p + 1;
	}
	gmp_printf("WARNING! p and q were NOT FOUND for N = %Zd!\n", N);	// if the loop is exited then failure has occured
	arr[0] = 0;															// 0s => failure
	arr[1] = 0;
	arr[2] = 0;
	mpz_clear(tmp_val[0]);
	mpz_clear(tmp_val[1]);
	return;
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
Runtime:	O(min(r * log(N)^2, log(log(N)) * log(N)^2))
*/
int primality_test_2_8 (long int A, long int r, long int n, short y) {
	mpz_t N; mpz_init(N);
	mpz_t tmp_val; mpz_init(tmp_val);
	mpz_ui_pow_ui (tmp_val, r, n);
	mpz_mul_ui (N, tmp_val, A);
	if (y == 1) {
		mpz_add_ui (N, N, 1);										// N = A * r^n + y
	}
	else {
		mpz_sub_ui (N, N, 1);										// N = A * r^n + y
	}
	if (mpz_cmp_ui (tmp_val, A/2) < 0 || y*y != 1 || A%2 != 0) {	// if A > 2r^n or gamma^2 != 1 or A is odd
		mpz_clear(N);
		mpz_clear(tmp_val);
		return -1;													// test doesn't apply
	}
	int divRes = trial_div(N, 0);
	if (divRes) {													// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
		if (mpz_cmpabs_ui (N, divRes) == 0) {
			mpz_clear(N);
			mpz_clear(tmp_val);
			return 1;
		}
		else {
			mpz_clear(N);
			mpz_clear(tmp_val);
			return 0;
		}
	}
	long int a = -1;
	if (log2(n) > r * log2(r)) {									// if log(n) > r(log(r))
		long int k = (r-1)/2;										// the k for calculating g_r(x)
		mpz_t g_r_coef[k + 1];										// has indexes 0 - k available
		while (k >= 0) {
			mpz_init (g_r_coef[k]);									// initialize coefficents for calculating g_r(x)
			k--;
		}
		int d = get_s0 (tmp_val, A, y, N);							// get s_0 and store in tmp_val (and the associated delta value used)
		if (d == 0 || !get_g_r_coef (g_r_coef, r, d, N)) {
			while (k >= 0) {
				mpz_clear (g_r_coef[k]);							// clear tmp_vals for calculating g_r(x)
				k--;
			}
			mpz_clear(N);
			mpz_clear(tmp_val);
			return 0;
		}
		long int i = 0;
		while ((0 != mpz_cmp_ui (tmp_val, 0)) && i++ < n) {			// while s_i != 0
			get_next_s_i (tmp_val, r, g_r_coef, N);					// finds s_i+1 and stores it in tmp_val
		}
		if (0 == mpz_cmp_ui (tmp_val, 0) && i <= n) {
			a = i;													// first value of i such that N !| s_i-1 and N | s_i (i has already been incremented)
		}
		k = (r-1)/2;												// the k for calculating g_r(x)
		while (k >= 0) {
			mpz_clear (g_r_coef[k]);								// clear tmp_vals for calculating g_r(x)
			k--;
		}
	}
	else {
		int pqd[3];
		find_p_q (pqd, N, y);
		a = bin_search_s_i (A, r, n, pqd, N);						// note that binary search is not necessarily faster (iterative method if log(n) > r*log(r))
	}
	if (a == -1) {													// if N !| s_a (s_n if no suitable a could be found)
		mpz_clear (N);
		mpz_clear (tmp_val);
		return 0;													// N is not prime
	}
	if (2 * a - n < 0) {
		//gmp_printf("%d * %d^%d + %d may be prime. Prime divisors would have the form 2kr^(%d) +/- 1.\n", A, r, n, y, a);
		mpz_clear (N);
		mpz_clear (tmp_val);
		return 2;
	}
	else {
		mpz_ui_pow_ui (tmp_val, r, 2 * a - n);
		mpz_mul_ui (tmp_val, tmp_val, 2);							//tmp_val = 2r^(2a-n)
		if (mpz_cmp_ui (tmp_val, A) <= 0) {							//if (A > 2r^(2i-n))
			//gmp_printf("%d * %d^%d + %d may be prime. Prime divisors would have the form %d * %d * %d^(%d) +/- 1.\n", A, r, n, y, r, (r-1)/2, r, a);
			mpz_clear (N);
			mpz_clear (tmp_val);
			return 2;												//N has form rk(r^a)+/-1
		}
		else {														// N is prime
			mpz_clear (N);
			mpz_clear (tmp_val);
			return 1;
		}
	}
	//gmp_printf("%d * %d^%d + %d is composite with prime divisors of the form 2kr^(%d) +/- 1.\n", A, r, n, y, a);
	mpz_clear (N);
	mpz_clear (tmp_val);
	return 0;
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
Runtime:	O(log(N)^2)
	2 inversions in get_r0_r1 dominate
*/
short primality_test_2_10 (long int A, long int n, short y) {
	short result = 1;												// assume the number is prime until found composite
	if (A % 2 == 0 || y*y != 1) {
		return -1;													// test doesn't apply
	}
	mpz_t N; mpz_init(N);
	mpz_t tmp_val; mpz_init(tmp_val);
	mpz_ui_pow_ui (tmp_val, 2, n);									// tmp_val = 2^n
	if (mpz_cmp_ui (tmp_val, A) < 0) {								// if A >= 2^n
		return -1;													// test doesn't apply
	}
	mpz_mul_ui (N, tmp_val, A);
	if (y == 1) {
		mpz_add_ui (N, N, 1);										// N = A * r^n + y
	}
	else {
		mpz_sub_ui (N, N, 1);										// N = A * r^n + y
	}
	int divRes = trial_div(N, 0);
	if (divRes) {													// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
		if (mpz_cmpabs_ui (N, divRes) == 0) {
		mpz_clear(N);
		mpz_clear(tmp_val);
			return 1;
		}
		else {
		mpz_clear(N);
		mpz_clear(tmp_val);
			return 0;
		}
	}
	mpz_t r_i; mpz_init (r_i);
	if (!get_r0_r1 (tmp_val, r_i, A, y, N)) {						// tmp_val will (temporarily) store r0, and r_i will store r1
		mpz_clear (r_i);
		mpz_clear (N);
		mpz_clear (tmp_val);
		return 0;													// if a multiplicative inverse doesn't exist mod N then N is not prime
	}
	long int a = -1;
	if (0 == mpz_cmp_ui (tmp_val, 0)) {								// if r0 = 0
		a = 0;
	}
	if (0 == mpz_cmp_ui (r_i, 0)) {									// if r1 = 0
		a = 1;
	}
	long int i = 1;
	while (i < n) {													// while a value hasn't been found yet
		mpz_mul (r_i, r_i, r_i);
		mpz_sub_ui (r_i, r_i, 2);									// r_(i+1) = r_i^2 - 2
		mpz_mod (r_i, r_i, N);
		i++;
		if (mpz_cmp_ui (r_i, 0) == 0) {								// if N | r_i
			a = i;													// a value found
			break;
		}
	}
	mpz_clear (r_i);
	if (a == -1) {													// if N !| r_i for i = 0, 1, ... , n
		mpz_clear (N);
		mpz_clear (tmp_val);
		return 0;													// then N is composite
	}
	mpz_ui_pow_ui (tmp_val, 2, 2*a-n);
	if (mpz_cmp_ui (tmp_val, A + 1) < 0) {
		//gmp_printf("%d * 2^%d + %d may be prime. Prime divisors would have the form k 2^(%d) +/- 1.\n", A, n, y, a);
		mpz_clear (N);
		mpz_clear (tmp_val);
		return 2;													// potentially prime
	}
	else {
		mpz_clear (N);
		mpz_clear (tmp_val);
		return 1;													// N is prime
	}
	//gmp_printf("%d * 2^%d + %d is composite with prime divisors of the form k 2^(%d) +/- 1.\n", A, n, y, a);
	mpz_clear (N);
	mpz_clear (tmp_val);
	return 0;
}







