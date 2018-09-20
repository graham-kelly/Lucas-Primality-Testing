#include <gmp.h>
#include <math.h>
#include "utils.h"
#include "roots.h"
#include "RWG2.h"
#include "RWG2S.h"

/*			Get S_i from RWG Section 2 (arbitrary i in 0...n)
Parameters:
	A	even positive integer coprime with r
	r	small odd prime integer
	i	integer in range (0...n)
	pqd	values of p, q, and d selected to correspond to N
	N	N = Ar^n+gamma

Returns:
	result	= s_i as defined in RWG Section 2 (just after Theorem 2.4)
*/
_Bool get_s_i (mpz_t result, long int A, long int r, long int i, int pqd[3], mpz_t N) {
	_Bool possibly_prime = 1;									// will be set to 0 if inversion fails (can occur if N is not prime)
	int p = 0, q = 1, d = 2;									// indices for the pqd array
	mpz_t constant; mpz_init(constant);
	mpz_t constant2; mpz_init(constant2);
	mpz_ui_pow_ui (constant, r, i);
	mpz_mul_ui (constant, constant, A);							// = Ar^i
	if (!get_u_i (result, constant, pqd[p], pqd[d], N)) {		// result = u_(Ar^i)
		possibly_prime = 0;
	}
	mpz_cdiv_q_ui (constant2, N, 2);							// constant2 = 1/2
	mpz_mul (constant, constant, constant2);					// constant = (A* r^i) / 2
	mpz_mod (constant, constant, N);
	mpz_set_ui (constant2, pqd[q]);
	possibly_prime = possibly_prime && mpz_invert (constant2, constant2, N);		// 1/q
	mpz_powm (constant, constant2, constant, N);				// constant = q ^ (-(A* r^i) / 2)
	mpz_mul (result, result, constant);							// result = u_(Ar^i) * q ^ (-(A* r^i) / 2)
	mpz_mod (result, result, N);								// result = s_i
	mpz_clear(constant);
	mpz_clear(constant2);
	return possibly_prime;
}

/*			Find smallest a such that N !| S_(a-1) and N | S_a
Parameters:
	A	even positive integer coprime with r
	r	small odd prime integer
	n	positive integer
	pqd	values of p, q, and d selected to correspond to N
	N	N = Ar^n+gamma

Returns:
	integer		smallest a such that N !| S_(a-1) and N | S_a (based on RWG Theorem 2.8)
				will be -1 if no value in 0 ... n is found

Runtime:		should be filled in probably
*/
long int bin_search_s_i (long int A, long int r, long int n, int pqd[3], mpz_t N) {
	mpz_t s_i; mpz_init(s_i);
	if (!get_s_i (s_i, A, r, n, pqd, N)) {						// if s_n is not 0 (mod N) then N is not prime (don't need anything else)
		mpz_clear (s_i);
		return -1;
	}
	if (mpz_cmp_ui (s_i, 0) != 0) {								// ensure that alpha in [0,n] exists
		mpz_clear (s_i);
		return -1;
	}
	else {
		long int lep = 0, rep = n, mpt;							// left and right endpoints (+ midpoint) for binary search
		while (rep - lep > 1) {
			mpt = (lep + rep - 1) / 2 + 1;						// mpt = (lep+rep)/2 if (lep-rep) is even otherwise mpt = (lep+rep)/2 + 1
			if (!get_s_i (s_i, A, r, mpt, pqd, N)) {			// s_mpt
				mpz_clear (s_i);
				return -1;
			}
			if (mpz_cmp_ui (s_i, 0) == 0) {						// if s_mpt == 0 then alpha is in left part
				rep = mpt;
			}
			else {												// if s_mpt != 0 then alpha is in right part
				lep = mpt;
			}
		}
		get_s_i (s_i, A, r, lep, pqd, N);						// done binary search, either lep or lep+1 = alpha
		if (mpz_cmp_ui (s_i, 0) == 0) {
			mpz_clear (s_i);
			return lep;
		}
		else {
			mpz_clear (s_i);
			return lep + 1;
		}
	}
	mpz_clear (s_i);
	return -1;
}

/*			Get S_0 from RWG Section 2
Parameters:
	A	even positive integer coprime with r
	y	gamma = +/-1
	N	N = Ar^n+gamma

Returns:
	s0	first value in S_i sequence defined in RWG Theorem 2.8 (if method fails this will be 0)
	integer 	delta value used to generate s0 (if method fails this will be 0)
*/
int get_s0 (mpz_t s0, long int A, short y, mpz_t N) {
	int arr[3] = {0, 0, 0};											// used to store p, q, d
	find_p_q (arr, N, y);
	int p = arr[0];
	int q = arr[1];
	int d = arr[2];
	if (d == 0) {													// if d=0 then no suitable p,q were found => this method fails as well
		mpz_set_ui (s0, 0);
		return d;
	}
	mpz_t constant; mpz_init (constant);
	mpz_set_ui (constant, A);
	if (! get_u_i (s0, constant, p, d, N)) {						// = u_i
		d = 0;
		mpz_set_ui (s0, 0);
	}
	if (q < 0) {
		q = -q;														// q = |q|
		if ((A/2) % 2) {											// if q^(A/2) < 0
			mpz_sub (s0, N, s0);									// set s0 = - s0 (mod N)
		}
	}
	mpz_ui_pow_ui (constant, q, A/2);								// note A is even
	mpz_mod (constant, constant, N);
	if (! mpz_invert (constant, constant, N)) {						// if not invertable => not prime (note: not prime =/> not invertable)
		d = 0;
		mpz_set_ui (s0, 0);
	}
	mpz_mul (s0, s0, constant);										// s0 = u_A / q^(A/2)
	mpz_mod (s0, s0, N);											// s0 (mod N)
	mpz_clear(constant);
	return d;
}

/*			Get S_(i+1) from RWG Section 2
Parameters:
	r	small odd prime integer
	constants		array of (r-1)/2 + 1 coefficients for the g_r(x) function
	N	N = Ar^n+gamma

Returns:
	s_(i+1)		= g_r(s_i), as defined in RWG section 2
*/
void get_next_s_i(mpz_t s_i, long int r, mpz_t constants[2], mpz_t g_r_constants[], mpz_t N) {
	mpz_set_ui (constants[0], 0);											// represents sum, (constants[1] represents current term)
	long int k = (r-1)/2;
	long int j = 0;
	while (j <= k) {
		mpz_powm_ui (constants[1], s_i, 2 * j + 1, N);						// term = s_i^(2j+1) (mod N)
		mpz_mul (constants[1], constants[1], g_r_constants[j]);
		mpz_add (constants[0], constants[0], constants[1]);									// max value for sum is k*N^2 => only use oen mod operation at the end
		j++;
	}
	mpz_mod (s_i, constants[0], N);											// set s_i = sum (mod N)
}

/*			Get corresponding lucas sequence term (RWG section 1)
Parameters:
	i	which term of the lucas sequence to get
	p	from u_i(p, q), coprime to q
	d	discriminant delta = p^2 - 4q (doesn't have to be >= 0)

Returns:
	u_i		ith lucas sequence term
*/
 _Bool get_u_i (mpz_t u_i, mpz_t i, int p, int d, mpz_t N) {
	_Bool possibly_prime = 1;
	if (d >= 0 && pow(floor(sqrt(d)), 2) == d) {			// if d is a perfect square
		long int rd = (long int) sqrt(d);
		mpz_t a; mpz_init(a);
		mpz_t b; mpz_init(b);
		mpz_cdiv_q_ui (b, N, 2);							// b = 1/2
		mpz_mul_ui (a, b, p + rd);							// = alpha
		mpz_mul_ui (b, b, p - rd);							// = beta
		mpz_powm (a, a, i, N);								// a = alpha^i
		mpz_powm (b, b, i, N);								// b = beta^i
		mpz_sub (a, a, b);									// = sqrt(d) * u_i
		mpz_set_ui (b, rd);
		possibly_prime = mpz_invert (b, b, N);				// b = (sqrt(d))^(-1)
		mpz_mul (u_i, a, b);
		mpz_mod (u_i, u_i, N);
		mpz_clear(a);
		mpz_clear(b);
	}
	else {
		mpz_t alpha[3]; mpz_init(alpha[0]); mpz_init(alpha[1]); mpz_init(alpha[2]);
		mpz_t constant; mpz_init(constant);
		mpz_cdiv_q_ui (constant, N, 2);						// constant = 1/2
		mpz_mul_ui (alpha[0], constant, p);
		mpz_set (alpha[1], constant);
		mpz_set_si (alpha[2], d);							// = alpha
		mpzrz_exp (constant, u_i, alpha, i, N);				// u_i = alpha^i, discard constant (= integer part of alpha^i which cancels with integer part of beta^i)
		mpz_mul_ui (u_i, u_i, 2);							// = u_i = 2* non-integer part of alpha^i (sqrt(delta) cancels with denominator)
		mpz_mod (u_i, u_i, N);
		mpz_clear(alpha[0]);
		mpz_clear(alpha[1]);
		mpz_clear(alpha[2]);
		mpz_clear(constant);
	}
	return possibly_prime;
}

/*			Get corresponding lucas sequence term (RWG section 1)
Parameters:
	i	which term of the lucas sequence to get
	p	from u_i(p, q), coprime to q
	d	discriminant delta = p^2 - 4q (doesn't have to be >= 0)

Returns:
	u_i		ith lucas sequence term
	v_i		ith lucas sequence term

Note takes much less time than calling both get_u_i and get_v_i separately
*/
_Bool get_uv_i (mpz_t u_i, mpz_t v_i, mpz_t i, int p, int d, mpz_t N) {
	_Bool possibly_prime = 1;
	if (d >= 0 && pow(floor(sqrt(d)), 2) == d) {			// if d is a perfect square
		long int rd = (long int) sqrt(d);
		mpz_t a; mpz_init(a);
		mpz_t b; mpz_init(b);
		mpz_cdiv_q_ui (b, N, 2);							// b = 1/2
		mpz_mul_ui (a, b, p + rd);							// = alpha
		mpz_mul_ui (b, b, p - rd);							// = beta
		mpz_powm (a, a, i, N);								// a = alpha^i
		mpz_powm (b, b, i, N);								// b = beta^i
		mpz_add (v_i, a, b);								// v_i = alpha^i + beta^i
		mpz_sub (a, a, b);									// = sqrt(d) * u_i
		mpz_set_ui (b, rd);
		possibly_prime = mpz_invert (b, b, N);				// b = (sqrt(d))^(-1)
		mpz_mul (u_i, a, b);
		mpz_mod (u_i, u_i, N);
		mpz_clear(a);
		mpz_clear(b);
	}
	else {
		mpz_t alpha[3]; mpz_init(alpha[0]); mpz_init(alpha[1]); mpz_init(alpha[2]);
		mpz_t constant; mpz_init(constant);
		mpz_cdiv_q_ui (constant, N, 2);						// constant = 1/2
		mpz_mul_ui (alpha[0], constant, p);
		mpz_set (alpha[1], constant);
		mpz_set_si (alpha[2], d);							// alpha = alpha
		mpzrz_exp (v_i, u_i, alpha, i, N);					// u_i / 2 = irrational part of alpha^i, v_i / 2 = integer part of alpha^i
		mpz_mul_ui (u_i, u_i, 2);							// = u_i
		mpz_mul_ui (v_i, v_i, 2);							// = v_i
		mpz_mod (u_i, u_i, N);
		mpz_mod (v_i, v_i, N);
		mpz_clear(alpha[0]);
		mpz_clear(alpha[1]);
		mpz_clear(alpha[2]);
		mpz_clear(constant);
	}
	return possibly_prime;
}

/*			Get corresponding lucas sequence term (RWG section 1)
Parameters:
	i	which term of the lucas sequence to get
	p	from u_i(p, q), coprime to q
	d	discriminant delta = p^2 - 4q (doesn't have to be >= 0)

Returns:
	v_i		ith lucas sequence term
*/
_Bool get_v_i (mpz_t v_i, mpz_t i, int p, int d, mpz_t N) {
	_Bool possibly_prime = 1;
	if (d >= 0 && pow(floor(sqrt(d)), 2) == d) {			// if d is a perfect square
		long int rd = (long int) sqrt(d);
		mpz_t constant; mpz_init(constant);
		mpz_set_ui (v_i, p + rd);							// = 2 * alpha
		mpz_cdiv_q_ui (constant, N, 2);						// constant = 1/2
		mpz_mul (v_i, v_i, constant);						// = alpha
		mpz_powm (v_i, v_i, i, N);							// = (alpha)^i
		mpz_mul_ui (v_i, v_i, 2);							// v_i = (alpha^i + beta^i)
		mpz_mod (v_i, v_i, N);
		mpz_clear(constant);
	}
	else {
		mpz_t alpha[3]; mpz_init(alpha[0]); mpz_init(alpha[1]); mpz_init(alpha[2]);
		mpz_t constant; mpz_init(constant);
		mpz_cdiv_q_ui (constant, N, 2);						// constant = 1/2
		mpz_mul_ui (alpha[0], constant, p);
		mpz_set (alpha[1], constant);
		mpz_set_si (alpha[2], d);							// = alpha
		mpzrz_exp (v_i, constant, alpha, i, N);				// v_i = alpha^i, constant = root part of v_i (discarded as it cancels with root part of beta^i)
		mpz_mul_ui (v_i, v_i, 2);							// alpha^i + beta^i
		mpz_mod (v_i, v_i, N);
		mpz_clear(alpha[0]);
		mpz_clear(alpha[1]);
		mpz_clear(alpha[2]);
		mpz_clear(constant);
	}
	return possibly_prime;
}

/*			Get r0 and r1 as defined in RWG Theorem 2.10
Parameters:
	A	even positive integer coprime with r
	y	gamma = +/-1
	N	N = Ar^n+gamma

Returns:
	r0, r1 as defined by RWG Theorem 2.10
	integer		delta value used to calculate r0, r1
*/
int get_r0_r1 (mpz_t r0, mpz_t r1, long int A, short y, mpz_t N) {
	int arr[3] = {0, 0, 0};											//used to store p, q, d
	find_p_q (arr, N, y);
	int p = arr[0];
	int q = arr[1];
	int d = arr[2];
	if (d == 0) {													// failure to find appropreate p, q results in failure of this method as well
		mpz_set_ui (r0, 0);
		mpz_set_ui (r1, 0);
		return d;
	}
	mpz_t constant; mpz_init (constant);
	mpz_set_ui (constant, A);
	if (! get_uv_i (r0, r1, constant, p, d, N)) {					// get u_A, v_A
		d = 0;
		mpz_set_ui (r0, 0);
		mpz_set_ui (r1, 0);
	}																//now multiply by q^(-A/2)
	mpz_t constant2; mpz_init (constant2);
	mpz_set_si (constant, q);
	if (! get_sqrt(constant2, constant, N)) {						//constant = sqrt(q) (mod N) (note: picked q so that Jacobi (q,N) == 1)
		d = 0;
		mpz_set_ui (r0, 0);
		mpz_set_ui (r1, 0);
	}
	if (! mpz_invert (constant, constant2, N)) {					//if not invertable => not prime (note: not prime =/> not invertable)
		d = 0;
		mpz_set_ui (r0, 0);
		mpz_set_ui (r1, 0);
	}																//constant = 1/sqrt(q)
	mpz_pow_ui (constant, constant, A);								//constant = 1/(q^(A/2))
	mpz_mul (r0, r0, constant);										//r0 = u_A / q^(A/2)
	mpz_mod (r0, r0, N);											//r0 (mod N)
	mpz_mul (r1, r1, constant);										//r1 = v_A / q^(A/2)
	mpz_mod (r1, r1, N);											//r1 (mod N)
	mpz_clear(constant);
	mpz_clear(constant2);
	return d;
}

