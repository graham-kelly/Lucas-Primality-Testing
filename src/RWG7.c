#include <gmp.h>
#include <math.h>
#include "utils.h"
#include "roots.h"
#include "RWG7.h"
#include "RWG7S.h"


/*			Evaluate the polynomial H_k(X,Y) (mod N)
Parameters:
	X	the first argument to the H_k polynomial defined in section 4 of RWG
	Y	the second argument to the H_k polynomial defined in section 4 of RWG
	k	index of H_k function defined in section 4 of RWG, called with (r-1)/2 in the primality tests 7.1-4
	N	number to evaluate H_k modulo, N = Ar^n + eta*gamma_n(r)
	
Returns:
	rop		the value of H_k(X, Y) (mod N)
*/
void get_H_k (mpz_t rop, mpz_t X, mpz_t Y, long int k, mpz_t N) {
	long int i;
	long int j = 0;
	mpz_t outer_sum; mpz_init (outer_sum);
	mpz_t inner_sum; mpz_init (inner_sum);
	mpz_t constant; mpz_init (constant);
	mpz_t constant2; mpz_init (constant2);
	mpz_set_ui (outer_sum, 0);
	while (j <= k) {
		i = 0;
		mpz_set_ui (inner_sum, 0);
		while (i <= j) {													// sum(i=0 -> j) binomial_coef(2j+1, 2i) * X^i * Y^(j-i)
			mpz_powm_ui (constant, X, i, N);								// X^i
			mpz_powm_ui (constant2, Y, j-i, N);								// Y^(j-i)
			mpz_mul (constant, constant, constant2);						// X^i * Y^(j-i)
			mpz_bin_uiui (constant2, 2*j+1, 2*i);							// binom (2j+1, 2i)
			mpz_mul (constant, constant, constant2);
			mpz_add (inner_sum, inner_sum, constant);						// inner_sum = (sum(i=0 -> j) binomial_coef (2j+1, 2i) * X^i * Y^(j-i)) / (2j+1)
			i++;
		}
		if ((k+j)%2 == 1) {													// *= (-1)^k+j
			mpz_neg (inner_sum, inner_sum);
		}
		mpz_bin_uiui (constant, k+j, k-j);									// binomial_coef (k+j, k-j)
		mpz_mul (inner_sum, inner_sum, constant);
		mpz_mul_ui (inner_sum, inner_sum, 2*k+1);							// *= 2k+1
		mpz_divexact_ui (inner_sum, inner_sum, 2*j+1);
		mpz_add (outer_sum, outer_sum, inner_sum);
		j++;
	}
	mpz_mod (rop, outer_sum, N);
	mpz_clear (outer_sum);
	mpz_clear (inner_sum);
	mpz_clear (constant);
	mpz_clear (constant2);
}

/*			Get appropreate values of Q, P1 and P2 for theorem 7.2/4
Parameters:
	N	N = Ar^n + eta*gamma_n(r). Selection must result in Jacobi(Delta, N) = Jacobi(E, N) = -1 where Delta = P1^2 - 4P2 and E = (P2+4Q)^2 - 4QP1^2
	
Returns:
	QPP		integer array [Q, P1, P2] such that Jacobi(Delta, N) = Jacobi(E, N) = -1 where Delta = P1^2 - 4P2 and E = (P2+4Q)^2 - 4QP1^2
*/
_Bool find_QPP (int QPP[3], mpz_t N) {
	QPP[0] = 1; QPP[1] = 1; QPP[2] = 1;
	int MAX_VAL = 100;
	mpz_t Delta; mpz_init (Delta);		// Delta = P1^2 - 4P2
	mpz_t E; mpz_init (E);				// E = (P2+4Q)^2 - 4QP1^2
	while (QPP[0] < MAX_VAL) {
		QPP[1] = 1;
		while (QPP[1] < MAX_VAL) {
			QPP[2] = 1;
			while (QPP[2] < MAX_VAL) {
				if (gcd (QPP[0], QPP[1]) == 1 && gcd (QPP[0], QPP[2]) == 1 && gcd (QPP[1], QPP[2]) == 1) {	// gcd (QPP) = 1
					mpz_set_si (Delta, QPP[1]*QPP[1] - 4*QPP[2]);											// Delta = P1^2 - 4P2
					mpz_set_si (E, (QPP[2] + 4*QPP[0]) * (QPP[2] + 4*QPP[0]) - 4*QPP[0]*QPP[1]*QPP[1]);		// E = (P2+4Q)^2 - 4QP1^2
					if ((mpz_jacobi(Delta, N) == -1) && (mpz_jacobi(E, N) == -1)) {							// if Jacobi (Delta/N) = Jacobi (E/N) = -1
						mpz_clear (Delta);
						mpz_clear (E);
						return 1;
					}
				}
				if (QPP[2] > 0) {
					QPP[2] = -QPP[2];				// P2 = 1 , -1, 2, -2, 3, ...
				}
				else {
					QPP[2] = -QPP[2] + 1;
				}
			}
			if (QPP[1] > 0) {
				QPP[1] = -QPP[1];					// P1 = 1 , -1, 2, -2, 3, ...
			}
			else {
				QPP[1] = -QPP[1] + 1;
			}
		}
		if (QPP[0] > 0) {
			QPP[0] = -QPP[0];						// Q = 1 , -1, 2, -2, 3, ...
		}
		else {
			QPP[0] = -QPP[0] + 1;
		}
	}
	mpz_clear (Delta);
	mpz_clear (E);
	return 0;
}

/*			Theorem 7.1 Primality Test from RWG
Parameters:
	A	an even positive integer not divisible by r
	r	a (small) odd prime
	n	a positive integer
	eta	= +/-1

Returns:
	integer		based on N = Ar^n + eta*gamma_n(r)
				-2 if no such number exists (gamma_n(r) cannot be computed (sqrt(-1) doesn't exist for r = 3 (mod 4)) )
				-1 if primality is unknown
				0 if composite
				1 if prime
	
Let m, n be positive integers with n >= m and let N be given by 6.3 (of RWG) where r !| A and gamma_m(r) !| N.
If A < 2r^(2m-n) and a prime divisor of p of N satisfies:		p^4 = 1 (mod r^m),			then N is a prime.
*/
int primality_test_7_1 (long int A, long int r, long int n, int eta) {
	int prime = -1;
	if (A % 2 != 0 || eta*eta != 1 || A % r == 0) {
		return prime;
	}
	mpz_t RST_m[3]; mpz_init (RST_m[0]); mpz_init (RST_m[1]); mpz_init (RST_m[2]);
	mpz_t constants[2]; mpz_init (constants[0]); mpz_init (constants[1]);
	mpz_t N; mpz_init (N);
	mpz_t gamma_n_r; mpz_init (gamma_n_r);
	mpz_t gamma_m_r; mpz_init (gamma_m_r);
	mpz_t rEXPn; mpz_init (rEXPn);
	mpz_t rEXPm; mpz_init (rEXPm);
	mpz_ui_pow_ui (rEXPn, r, n);															// r^n
	mpz_set_si (constants[0], -1);
	mpz_set_ui (constants[1], r);
	if (! h_lift_root (gamma_n_r, constants[0], constants[1], n)) {							// = sqrt(-1) if r^n = 1 (mod 4)
		prime = -2;
	}
	else {
		if (mpz_tstbit(gamma_n_r, 0) == 0) {												// if root found is even
			mpz_neg (gamma_n_r, gamma_n_r);
			mpz_add (gamma_n_r, gamma_n_r, rEXPn);											// to get odd sqrt(-1)
		}																					// gamma_n(r) = sqrt(-1) (mod r^n)
		mpz_mul_ui (N, rEXPn, A);															// = Ar^n
		mpz_mul_si (constants[0], gamma_n_r, eta);
		mpz_add (N, N, constants[0]);														// computed N
//if (mpz_probab_prime_p (N, 20)) {gmp_printf("n = %d gives a likely prime\n", n);}
		if (!trial_div (N, 0)) {															// trial division by first million primes
			long int m = (long int) (log(A/2)/log(r) + n)/2;								// start with m = ceil ((log_r(A/2) + n)/2)
			int QPP[3] = {3,2,1};															// pick QPP, gcd(N,Q) = 1, gcd (P1, P2, Q) = 1
			if (get_RST_i (RST_m, m, QPP, r, rEXPn, N)) {
				m++;
				while (m <= n) {
					mpz_ui_pow_ui (rEXPm, r, m);											// do using binomial theorem ?
					mpz_set_si (constants[0], -1);
					mpz_set_ui (constants[1], r);
					if (h_lift_root (gamma_m_r, constants[0], constants[1], m)) {			// = sqrt(-1) mod r^m
						if (mpz_tstbit(gamma_m_r, 0) == 0) {								// if root found is even
							mpz_neg (gamma_m_r, gamma_m_r);
							mpz_add (gamma_m_r, gamma_m_r, rEXPm);							// to get odd sqrt(-1)
						}																	// gamma_n(r) = sqrt(-1) (mod r^n)
						if (mpz_cmp_ui (gamma_m_r, 1) != 0 && mpz_divisible_p(N, gamma_m_r)) {		// check that gamma_m(r) !| N
							prime = 0;
							break;
						}
					}
					get_next_RST_i (RST_m, RST_m, constants, QPP, r, N);
					mpz_mul (constants[0], RST_m[2], RST_m[2]);								// = (T_m)^2 (if N is prime will be 4 (mod N))
					mpz_set_ui (constants[1], 4);
					if (mpz_congruent_p(constants[0], constants[1], N) && mpz_divisible_p(RST_m[0], N)){	// if (T_m)^2 = 4 (mod N) && N | R_m
						prime = 1;
						break;																// then N is prime			
					}
					m++;
				}
			}
			else {
				prime = 0;
			}
		}
		else {
			prime = 0;
		}
	}
	mpz_clear (RST_m[0]); mpz_clear (RST_m[1]); mpz_clear (RST_m[2]);
	mpz_clear (constants[0]); mpz_clear (constants[1]);
	mpz_clear (N);
	mpz_clear (gamma_n_r);
	mpz_clear (gamma_m_r);
	mpz_clear (rEXPn);
	mpz_clear (rEXPm);
	return prime;																			// 1 if prime, 0 if composite, -1 if unknown, -2 if N doesn't exist
}

/*			Theorem 7.2 and 7.4 Primality Tests from RWG
Parameters:
	A	an even positive integer not divisible by r
	r	a (small) odd prime
	n	a positive integer
	eta	= +/-1

Returns:
	integer		based on N = Ar^n + eta*gamma_n(r)
				-2 if no such number exists (gamma_n(r) cannot be computed (sqrt(-1) doesn't exist for r = 3 (mod 4)) )
				-1 if primality is unknown
				0 if composite
				1 if prime

7.2		Let r be an odd prime and N = Ar^n + eta*gamma_n(r) where eta^2 = 1, 2|A, A<2r^n, and gcd(N,Q) = 1.
		Suppose that Jacobi(Delta, N) = Jacobi(E, N) = -1. Define the sequences {R_i}, {S_i}, {T_i} by 6.4 - 6.7
		of RWG. If N !| (T_n^2 - 4, S_n), then N is composite. If alpha (<= n) is the least positive value
		of i such that N | (T_i^2 - 4, S_i) and N !| S_(i-1), then a prime divisor p of N must satisfy:
		p^4 = 1 (mod r^alpha). Also if A <= 2r^(2alpha-n) then n is prime.

7.4		Let r be an odd prime and N = Ar^n + eta*gamma_n(r) where eta^2 = 1, 2|A, A<2r^n, and gcd(N,Q) = 1.
		Suppose that Jacobi(Delta, N) = Jacobi(E, N) = -1. Define the sequences {R_i}, {S_i}, {T_i} by 6.4 - 6.7
		of RWG. If either N !| S_n or N !| T_i^2 - 4 whenever N | R_i (i in 0, ... n) then N is composite.
		If alpha is the least value of i such that N | T_i^2 - 4 and N | R_i then any prime divisor p of N
		must satisfy p^4 = 1 (mod r^alpha). Also if A <= 2r^(2alpha - n) then N is prime
		
*/
int primality_test_7_2_4 (long int A, long int r, long int n, int eta) {
	int prime = -1;
	if (A % 2 != 0 || eta*eta != 1 || A % r == 0) {
		return prime;
	}
	mpz_t N; mpz_init (N);
	mpz_t rEXPn; mpz_init (rEXPn);
	mpz_ui_pow_ui (rEXPn, r, n);															// r^n
	if (mpz_cmp_ui (rEXPn, A/2) <= 0) {
		return prime;
	}
	mpz_t gamma_n_r; mpz_init (gamma_n_r);
	mpz_t constants[2]; mpz_init (constants[0]); mpz_init (constants[1]);
	mpz_set_si (constants[0], -1);
	mpz_set_ui (constants[1], r);
	if (! h_lift_root (gamma_n_r, constants[0], constants[1], n)) {							// = sqrt(-1) if r^n = 1 (mod 4)
		prime = -2;
	}
	else {
		if (mpz_tstbit(gamma_n_r, 0) == 0) {												// if root found is even
			mpz_neg (gamma_n_r, gamma_n_r);
			mpz_add (gamma_n_r, gamma_n_r, rEXPn);											// to get odd sqrt(-1)
		}																					// gamma_n(r) = sqrt(-1) (mod r^n)
		mpz_mul_ui (N, rEXPn, A);															// = Ar^n
		mpz_mul_si (constants[0], gamma_n_r, eta);
		mpz_add (N, N, constants[0]);														// computed N
//************************************************************************************************************************************************************************************
//if (mpz_probab_prime_p (N, 20)) {gmp_printf("n = %d gives a likely prime\n", n);}			// used for testing!
//************************************************************************************************************************************************************************************
		if (!trial_div (N, 0)) {															// trial division by first million primes
			int QPP[3];
			if (find_QPP (QPP, N)) {
				mpz_t RST_i[3]; mpz_init (RST_i[0]); mpz_init (RST_i[1]); mpz_init (RST_i[2]);
				mpz_t S_im1; mpz_init (S_im1);
				if (get_RST_i (RST_i, 0, QPP, r, rEXPn, N)) {
					long int alpha = 0;
					while (alpha++ < n) {														// alpha = 1 ... n
						mpz_set (S_im1, RST_i[1]);
						get_next_RST_i (RST_i, RST_i, constants, QPP, r, N);
						mpz_mul (constants[0], RST_i[2], RST_i[2]);								// = (T_i)^2 (if N is prime will be 4 (mod N))
						mpz_sub_ui (constants[0], constants[0], 4);
						if (mpz_divisible_p(RST_i[0], N) && !mpz_divisible_p(constants[0], N)) {	// based on Theorem 7.4
							prime = 0;
							break;
						}
						if (mpz_divisible_p(constants[0], N)) {									// if N | (T_i)^2 - 4
							if (mpz_divisible_p(RST_i[1], N) && !mpz_divisible_p(S_im1, N) || mpz_divisible_p(RST_i[0], N)) {	// if N | S_i && N !| S_i-1 or N | R[i]
								break;
							}
						}
					}
					if (log(A/2)/log(r) < 2*alpha - n) {									// N is prime (if N | T_n^2 -4 and S_n)
						prime = 1;
					}
					else {																	// may be prime... prime divisor must satisfy p^4 = 1 (mod r^alpha)
						gmp_printf ("A prime divisor of N = %d*%d^%d + (%d)gamma_n(r) must satisfy p^4 = 1 (mod %d^%d\n", A, r, n, eta, r, alpha);
					}
					if (prime != 0) {
						while (alpha++ < n) {
							get_next_RST_i (RST_i, RST_i, constants, QPP, r, N);
							mpz_mul (constants[0], RST_i[2], RST_i[2]);								// = (T_i)^2 (if N is prime will be 4 (mod N) if N | R_i)
							mpz_sub_ui (constants[0], constants[0], 4);
							if (mpz_divisible_p(RST_i[0], N) && !mpz_divisible_p(constants[0], N)) {	// based on Theorem 7.4
								prime = 0;
								break;
							}
						}
						mpz_mul (constants[0], RST_i[2], RST_i[2]);											// = (T_i)^2 (if N is prime will be 4 (mod N))
						mpz_sub_ui (constants[0], constants[0], 4);
						if (!mpz_divisible_p(constants[0], N) || !mpz_divisible_p(RST_i[1], N)) {				// check that N | T_n^2 -4 and S_n
							prime = 0;
						}
					}
				}
				else {
					prime = 0;
				}
				mpz_clear (RST_i[0]); mpz_clear (RST_i[1]); mpz_clear (RST_i[2]);
				mpz_clear (S_im1);
			}
			else {
				prime = 0;
			}
		}
		else {
			prime = 0;
		}
	}	
	return prime;
	
	
	
	
}

/*			Theorem 7.5 Primality Test from RWG
Parameters:
	A	an even positive integer not divisible by r
	n	a positive integer
	eta	= +/-1

Returns:
	integer		based on N = A*5^n + eta*gamma_n(5)
				-1 if primality is unknown
				0 if composite
				1 if prime

7.5		Let N = A*5^n + eta*gamma_n(5) where eta^2 = 1, 2|A and A<2*5^n. Suppose that
		q is a prime such that q = 1 (mod 5) and N^((q-1)/5) != 1 (mod q) and let P1, 
		P2, Q be defined as in 7.1 of RWG. Put R0 = S0 = U_(N^2 + 1) / 2Q^((N^2+1)/(2*5^n))
		T0 = U_(N^2 + 1) / 2Q^((N^2+1)/(2*5^n)) (mod N) and define
			S_(i+1) = S_i(Delta^2 * S_i^4 + 10*Delta*S_i^2 * T_i^2 + 5T_i^4 - 5*Delta*S_i^2 - 15T_i^2 +5) (mod N)
			T_(i+1) = T_i(T_i^4 + 10*Delta*S_i^2 * T_i^2 + 5*Delta^2 * S_i^4 - 5*T_i^2 - 15*Delta*S_i^2 +5) (mod N)
		for i = 0 , ... n-1, then N is prime iff 
			4*Delta*(S_(n-1))^2 = 5 (mod N) and 2*T_(n-1) = -1 (mod N)

*/
int primality_test_7_5 (long int A, long int n, int eta) {
	return 0;
}




