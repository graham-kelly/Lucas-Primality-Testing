#include <gmp.h>
#include <math.h>
#include "utils.h"
#include "roots.h"
#include "RWG7.h"
#include "RWG7S.h"

static const int P_ipq_table_size = 12;
static const int P_ipq_table[12][3] =	{{11,-89,1199},{31,-409,22289},{41,981,239809},{61,1111,214049},		// format is [q, P(1,5,q), P(2,5,q)]
										{71,101,-1310731},{101,271,-1294921},{131,-4009,3735989},
										{151,596,-4423696},{181,1691,-7254661},{191,1331,-18326641},
										{211,961,-24801151},{241,-3344,1283084}};

/*			Evaluate the polynomial H_k(X,Y), I_k(X,Y) (mod N)
Parameters:
	X	the first argument to the H_k polynomial defined in section 4 of RWG
	Y	the second argument to the H_k polynomial defined in section 4 of RWG
	k	index of H_k function defined in section 4 of RWG, called with (r-1)/2 in the primality tests 7.1-4
	N	number to evaluate H_k modulo, N = Ar^n + eta*gamma_n(r)
	
Returns:
	rop1	the value of H_k(X, Y) (mod N)
	rop2	the value of I_k(X, Y) = H_k(Y,X) (mod N)
Runtime:	O(log(r) * r^2 * M(N))
	k = (r-1)/2 and most expensive operation is exponentiation in nested loop
*/
void get_HI_k (mpz_t rop1, mpz_t rop2, mpz_t X, mpz_t Y, int k, mpz_t N) {
	int i;
	int j = 0;
	mpz_t outer_sum1; mpz_init (outer_sum1);
	mpz_t inner_sum1; mpz_init (inner_sum1);
	mpz_t outer_sum2; mpz_init (outer_sum2);
	mpz_t inner_sum2; mpz_init (inner_sum2);
	mpz_t tmp_val[2]; mpz_init (tmp_val[0]); mpz_init (tmp_val[1]);
	mpz_set_ui (outer_sum1, 0);
	mpz_set_ui (outer_sum2, 0);
	while (j <= k) {
		i = 0;
		mpz_set_ui (inner_sum1, 0);
		mpz_set_ui (inner_sum2, 0);
		while (i <= j) {													// sum(i=0 -> j) binomial_coef(2j+1, 2i) * X^i * Y^(j-i)
			mpz_powm_ui (tmp_val[0], X, i, N);								// X^i
			mpz_powm_ui (tmp_val[1], Y, j-i, N);							// Y^(j-i)
			mpz_mul (tmp_val[0], tmp_val[0], tmp_val[1]);					// X^i * Y^(j-i)
			mpz_bin_uiui (tmp_val[1], 2*j+1, 2*i);							// binom (2j+1, 2i)
			mpz_mul (tmp_val[1], tmp_val[0], tmp_val[1]);
			mpz_add (inner_sum1, inner_sum1, tmp_val[1]);					// inner_sum1 = (sum(i=0 -> j) binomial_coef (2j+1, 2i) * X^i * Y^(j-i)) / (2j+1)
			mpz_bin_uiui (tmp_val[1], 2*j+1, 2*i+1);						// binom (2j+1, 2i+1)
			mpz_mul (tmp_val[1], tmp_val[0], tmp_val[1]);
			mpz_add (inner_sum2, inner_sum2, tmp_val[1]);					// inner_sum2 = (sum(i=0 -> j) binomial_coef (2j+1, 2i+1) * X^i * Y^(j-i)) / (2j+1)
			i++;
		}
		if ((k+j)%2 == 1) {													// *= (-1)^k+j
			mpz_neg (inner_sum1, inner_sum1);
			mpz_neg (inner_sum2, inner_sum2);
		}
		mpz_bin_uiui (tmp_val[0], k+j, k-j);								// binomial_coef (k+j, k-j)
		mpz_mul (inner_sum1, inner_sum1, tmp_val[0]);
		mpz_mul_ui (inner_sum1, inner_sum1, 2*k+1);							// *= 2k+1
		mpz_divexact_ui (inner_sum1, inner_sum1, 2*j+1);
		mpz_add (outer_sum1, outer_sum1, inner_sum1);
		mpz_mul (inner_sum2, inner_sum2, tmp_val[0]);
		mpz_mul_ui (inner_sum2, inner_sum2, 2*k+1);							// *= 2k+1
		mpz_divexact_ui (inner_sum2, inner_sum2, 2*j+1);
		mpz_add (outer_sum2, outer_sum2, inner_sum2);
		j++;
	}
	mpz_mod (rop1, outer_sum1, N);
	mpz_mod (rop2, outer_sum2, N);
	mpz_clear (outer_sum1);
	mpz_clear (inner_sum1);
	mpz_clear (outer_sum2);
	mpz_clear (inner_sum2);
	mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
}

#define MAX_ABSOLUTE_Q_P1_P2_SIZE 100;

/*			Get appropreate values of Q, P1 and P2 for theorem 7.2 / 7.4
Parameters:
	N	N = Ar^n + eta*gamma_n(r). Selection must result in Jacobi(Delta, N) = Jacobi(E, N) = -1 where Delta = P1^2 - 4P2 and E = (P2+4Q)^2 - 4QP1^2
	
Returns:
	QPP		integer array [Q, P1, P2] such that Jacobi(Delta, N) = Jacobi(E, N) = -1 where Delta = P1^2 - 4P2 and E = (P2+4Q)^2 - 4QP1^2

Runtime:	O(1)
	jacobi with small integers (absolute size < 100 each)
*/
_Bool find_QPP_7_2_4 (int QPP[3], mpz_t N) {
	if (try_best_QPP_7_2_4(QPP, N)) {	// if this suceeds can avoid get_KL_m method which is slow (same asymptotic run time but much lower constant)
		return 1;
	}
	int MAX_VAL = MAX_ABSOLUTE_Q_P1_P2_SIZE;
	mpz_t Delta; mpz_init (Delta);		// Delta = P1^2 - 4P2
	mpz_t E; mpz_init (E);				// E = (P2+4Q)^2 - 4QP1^2
	for (QPP[0] = 1; QPP[0] < MAX_VAL;) {
		if (isSquare(QPP[0])) {			// already checked all of these
			continue;
		}
		for (QPP[1] = 1; QPP[1] < MAX_VAL;) {
			for (QPP[2] = 1; QPP[2] < MAX_VAL;) {
				if (gcd (QPP[0], QPP[1]) == 1 && gcd (QPP[0], QPP[2]) == 1 && gcd (QPP[1], QPP[2]) == 1) {	// gcd (QPP) = 1
					mpz_set_si (Delta, QPP[1]*QPP[1] - 4*QPP[2]);											// Delta = P1^2 - 4P2
					mpz_set_si (E, (QPP[2] + 4*QPP[0]) * (QPP[2] + 4*QPP[0]) - 4*QPP[0]*QPP[1]*QPP[1]);		// E = (P2+4Q)^2 - 4QP1^2
					if ((mpz_jacobi(Delta, N) == -1) && (mpz_jacobi(E, N) == -1)) {							// if Jacobi (Delta/N) = Jacobi (E/N) = -1
						mpz_clear (Delta);
						mpz_clear (E);
						return 1;
					}
				}
				QPP[2] = -QPP[2];
				if (QPP[2] > 0) {
					QPP[2] += 1;					// P2 = 1 , -1, 2, -2, 3, ...
				}
			}
			QPP[1] = -QPP[1];
			if (QPP[1] > 0) {
				QPP[1] += 1;						// P1 = 1 , -1, 2, -2, 3, ...
			}
		}
		QPP[0] = -QPP[0];
		if (QPP[0] > 0) {
			QPP[0] += 1;							// Q = 1 , -1, 2, -2, 3, ...
		}
	}
	mpz_clear (Delta);
	mpz_clear (E);
	return 0;
}

/*			Get appropreate values of Q, P1 and P2 for theorem 7.2 / 7.4
Parameters:
	N	N = Ar^n + eta*gamma_n(r). Selection must result in Jacobi(Delta, N) = Jacobi(E, N) = -1 where Delta = P1^2 - 4P2 and E = (P2+4Q)^2 - 4QP1^2
	
Returns:
	QPP		integer array [Q, P1, P2] such that Jacobi(Delta, N) = Jacobi(E, N) = -1 where Delta = P1^2 - 4P2 and E = (P2+4Q)^2 - 4QP1^2

Runtime:	O(1)
	jacobi with small integers (absolute size < 100 each)
*/
_Bool try_best_QPP_7_2_4(int QPP[3], mpz_t N) {
	mpz_t Delta; mpz_init (Delta);		// Delta = P1^2 - 4P2
	mpz_t E; mpz_init (E);				// E = (P2+4Q)^2 - 4QP1^2
	int MAX_VAL = MAX_ABSOLUTE_Q_P1_P2_SIZE;
	for (int i = 0; i <= MAX_VAL; i++) {
		QPP[0] = i * i;
		for (QPP[1] = 1; QPP[1] < MAX_VAL;) {
			for (QPP[2] = 1; QPP[2] < MAX_VAL;) {
				if (gcd (QPP[0], QPP[1]) == 1 && gcd (QPP[0], QPP[2]) == 1 && gcd (QPP[1], QPP[2]) == 1) {	// gcd (QPP) = 1
					mpz_set_si (Delta, QPP[1]*QPP[1] - 4*QPP[2]);											// Delta = P1^2 - 4P2
					mpz_set_si (E, (QPP[2] + 4*QPP[0]) * (QPP[2] + 4*QPP[0]) - 4*QPP[0]*QPP[1]*QPP[1]);		// E = (P2+4Q)^2 - 4QP1^2
					if ((mpz_jacobi(Delta, N) == -1) && (mpz_jacobi(E, N) == -1)) {							// if Jacobi (Delta/N) = Jacobi (E/N) = -1
						mpz_clear (Delta);
						mpz_clear (E);
						return 1;
					}
				}
				QPP[2] = -QPP[2];
				if (QPP[2] > 0) {
					QPP[2] += 1;					// P2 = 1 , -1, 2, -2, 3, ...
				}
			}
			QPP[1] = -QPP[1];
			if (QPP[1] > 0) {
				QPP[1] += 1;						// P1 = 1 , -1, 2, -2, 3, ...
			}
		}
	}
	mpz_clear (Delta);
	mpz_clear (E);
	return 0;
}

/*			Use table lookup to find values of P(1,5,q) and P(2,5,q)
Parameters:
	q		small prime = 1 (mod 5) for which to find P(1,5,q) and P(2,5,q)
	N		prime to be checked, will be used for one check if q is suitable
	
Returns:
	QPP[0]	q^3
	QPP[1]	value of P_1 = P(1,5,q) will be here when function returns true
	QPP[2]	value of P_2 = P(2,5,q) will be here when function returns true
	_Bool	true if P(1,5,q) and P(2,5,q) were found sucessfully
Runtime:
	Constant currently but will fail for q > 241 (chance of failure = 10^-22)
	************************ (can implement algorithm for larger q) ************************
*/
_Bool find_QPP_7_5 (int QPP[3], mpz_t tmp_val, mpz_t N) {
	int q;
	for (int i = 0; i < P_ipq_table_size; i++) {	// check array above
		q = P_ipq_table[i][0];
		mpz_pow_ui (tmp_val, N, (q-1)/5);
		mpz_sub_ui (tmp_val, tmp_val, 1);
		if (!mpz_divisible_ui_p (tmp_val, q)) {
			QPP[0] = q*q*q;					// Q
			QPP[1] = P_ipq_table[i][1];		// P_1
			QPP[2] = P_ipq_table[i][2];		// P_2
			return 1;
		}
	}
	QPP[0] = 0;
	QPP[1] = 0;
	QPP[2] = 0;
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
Runtime:	
*/
/*
int primality_test_7_1 (int A, int r, int n, int eta) {
	int prime = -1;
	if (A % 2 != 0 || eta*eta != 1 || A % r == 0) {
		return prime;
	}
	mpz_t RST_m[3]; mpz_init (RST_m[0]); mpz_init (RST_m[1]); mpz_init (RST_m[2]);
	mpz_t tmp_val[2]; mpz_init (tmp_val[0]); mpz_init (tmp_val[1]);
	mpz_t N; mpz_init (N);
	mpz_t gamma_n_r; mpz_init (gamma_n_r);
	mpz_t gamma_m_r; mpz_init (gamma_m_r);
	mpz_t rEXPn; mpz_init (rEXPn);
	mpz_t rEXPm; mpz_init (rEXPm);
	mpz_ui_pow_ui (rEXPn, r, n);															// r^n
	mpz_set_si (tmp_val[0], -1);
	mpz_set_ui (tmp_val[1], r);
	if (! h_lift_root (gamma_n_r, tmp_val[0], tmp_val[1], n)) {							// = sqrt(-1) if r^n = 1 (mod 4)
		prime = -2;
	}
	else {
		if (mpz_tstbit(gamma_n_r, 0) == 0) {												// if root found is even
			mpz_neg (gamma_n_r, gamma_n_r);
			mpz_add (gamma_n_r, gamma_n_r, rEXPn);											// to get odd sqrt(-1)
		}																					// gamma_n(r) = sqrt(-1) (mod r^n)
		mpz_mul_ui (N, rEXPn, A);															// = Ar^n
		mpz_mul_si (tmp_val[0], gamma_n_r, eta);
		mpz_add (N, N, tmp_val[0]);														// computed N
//if (mpz_probab_prime_p (N, 20)) {gmp_printf("n = %d gives a likely prime\n", n);}
		if (!trial_div (N, 0)) {															// trial division by first million primes
			int m = (int) (log(A/2)/log(r) + n)/2;								// start with m = ceil ((log_r(A/2) + n)/2)
			int QPP[3] = {3,2,1};															// pick QPP, gcd(N,Q) = 1, gcd (P1, P2, Q) = 1
			if (get_RST_i (RST_m, m, QPP, A, r, rEXPn, gamma_n_r, eta, N)) {
				m++;
				while (m <= n) {
					mpz_ui_pow_ui (rEXPm, r, m);											// do using binomial theorem ?
					mpz_set_si (tmp_val[0], -1);
					mpz_set_ui (tmp_val[1], r);
					if (h_lift_root (gamma_m_r, tmp_val[0], tmp_val[1], m)) {			// = sqrt(-1) mod r^m
						if (mpz_tstbit(gamma_m_r, 0) == 0) {								// if root found is even
							mpz_neg (gamma_m_r, gamma_m_r);
							mpz_add (gamma_m_r, gamma_m_r, rEXPm);							// to get odd sqrt(-1)
						}																	// gamma_n(r) = sqrt(-1) (mod r^n)
						if (mpz_cmp_ui (gamma_m_r, 1) != 0 && mpz_divisible_p(N, gamma_m_r)) {		// check that gamma_m(r) !| N
							prime = 0;
							break;
						}
					}
					get_next_RST_i (RST_m, RST_m, tmp_val, QPP, r, N);
					mpz_mul (tmp_val[0], RST_m[2], RST_m[2]);								// = (T_m)^2 (if N is prime will be 4 (mod N))
					mpz_set_ui (tmp_val[1], 4);
					if (mpz_congruent_p(tmp_val[0], tmp_val[1], N) && mpz_divisible_p(RST_m[0], N)){	// if (T_m)^2 = 4 (mod N) && N | R_m
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
	mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
	mpz_clear (N);
	mpz_clear (gamma_n_r);
	mpz_clear (gamma_m_r);
	mpz_clear (rEXPn);
	mpz_clear (rEXPm);
	return prime;																			// 1 if prime, 0 if composite, -1 if unknown, -2 if N doesn't exist
}
*/

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
		
Runtime:	log(log(N)) * log(N)^2
	from hensel lifting to find N
*/
int primality_test_7_2_4 (int A, int r, int n, int eta) {
	if (A % 2 != 0 || eta*eta != 1 || A % r == 0) {
		return -1;
	}
	mpz_t rEXPn; mpz_init (rEXPn);
	mpz_ui_pow_ui (rEXPn, r, n);															// r^n
	if (mpz_cmp_ui (rEXPn, A/2) <= 0) {
		mpz_clear (rEXPn);
		return -1;
	}
	mpz_t gamma_n_r; mpz_init (gamma_n_r);
	mpz_t tmp_val[2]; mpz_init (tmp_val[0]); mpz_init (tmp_val[1]);
	mpz_set_si (tmp_val[0], -1);
	mpz_set_ui (tmp_val[1], r);
	if (! h_lift_root (gamma_n_r, tmp_val[0], tmp_val[1], n)) {							// = sqrt(-1) if r^n = 1 (mod 4)
		mpz_clear (rEXPn);
		mpz_clear (gamma_n_r);
		mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
		return -2;
	}
	if (mpz_tstbit(gamma_n_r, 0) == 0) {												// if root found is even
		mpz_neg (gamma_n_r, gamma_n_r);
		mpz_add (gamma_n_r, gamma_n_r, rEXPn);											// to get odd sqrt(-1)
	}																					// gamma_n(r) = sqrt(-1) (mod r^n)
	mpz_t N; mpz_init (N);
	mpz_mul_ui (N, rEXPn, A);															// = Ar^n
	mpz_mul_si (tmp_val[0], gamma_n_r, eta);
	mpz_add (N, N, tmp_val[0]);														// computed N
//************************************************************************************************************************************************************************************
//if (mpz_probab_prime_p (N, 20)) {gmp_printf("n = %d gives a likely prime\n", n);}
//************************************************************************************************************************************************************************************
	int divRes = trial_div(N, 0);
	if (divRes) {													// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
		if (mpz_cmpabs_ui (N, divRes) == 0) {
			mpz_clear(N);
			mpz_clear (rEXPn);
			mpz_clear (gamma_n_r);
			mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
			return 1;
		}
		else {
			mpz_clear(N);
			mpz_clear (rEXPn);
			mpz_clear (gamma_n_r);
			mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
			return 0;
		}
	}
	int QPP[3];
	if (!find_QPP_7_2_4 (QPP, N)) {
		mpz_clear (N);
		mpz_clear (rEXPn);
		mpz_clear (gamma_n_r);
		mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
		return 0;
	}
	mpz_t RST_i[3]; mpz_init (RST_i[0]); mpz_init (RST_i[1]); mpz_init (RST_i[2]);
	mpz_t S_im1; mpz_init (S_im1);
	if (!get_RST_i (RST_i, 0, QPP, A, r, rEXPn, gamma_n_r, eta, N)) {
		mpz_clear (N);
		mpz_clear (rEXPn);
		mpz_clear (gamma_n_r);
		mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
		mpz_clear (RST_i[0]); mpz_clear (RST_i[1]); mpz_clear (RST_i[2]);
		mpz_clear (S_im1);
		return 0;
	}
	int alpha = 0;
	while (alpha++ < n) {														// alpha = 1 ... n
		mpz_set (S_im1, RST_i[1]);
		get_next_RST_i (RST_i, RST_i, tmp_val, QPP, r, N);
		mpz_mul (tmp_val[0], RST_i[2], RST_i[2]);								// = (T_i)^2 (if N is prime will be 4 (mod N))
		mpz_sub_ui (tmp_val[0], tmp_val[0], 4);
		if (mpz_divisible_p(RST_i[0], N) && !mpz_divisible_p(tmp_val[0], N)) {	// based on Theorem 7.4
			mpz_clear (N);
			mpz_clear (rEXPn);
			mpz_clear (gamma_n_r);
			mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
			mpz_clear (RST_i[0]); mpz_clear (RST_i[1]); mpz_clear (RST_i[2]);
			mpz_clear (S_im1);
			return 0;
		}
		if (mpz_divisible_p(tmp_val[0], N)) {									// if N | (T_i)^2 - 4
			if ((mpz_divisible_p(RST_i[1], N) && !mpz_divisible_p(S_im1, N)) || mpz_divisible_p(RST_i[0], N)) {	// if N | S_i && N !| S_i-1 or N | R[i]
				break;
			}
		}
	}
	if (log(A/2)/log(r) < 2*alpha - n) {
		while (alpha++ < n) {
			get_next_RST_i (RST_i, RST_i, tmp_val, QPP, r, N);
		}
		mpz_mul (tmp_val[0], RST_i[2], RST_i[2]);											// = (T_n)^2 (if N is prime will be 4 (mod N))
		mpz_sub_ui (tmp_val[0], tmp_val[0], 4);
		mpz_clear (rEXPn);
		mpz_clear (gamma_n_r);
		mpz_clear (tmp_val[1]);
		mpz_clear (RST_i[0]); mpz_clear (RST_i[2]);
		mpz_clear (S_im1);
		if (mpz_divisible_p(tmp_val[0], N) && mpz_divisible_p(RST_i[1], N)) {				// if N | T_n^2 -4 and S_n
			mpz_clear (N);
			mpz_clear (tmp_val[0]);
			mpz_clear (RST_i[1]); 
			return 1;																		// N is prime
		}
		else {
			mpz_clear (N);
			mpz_clear (tmp_val[0]);
			mpz_clear (RST_i[1]); 
			return 0;
		}
	}
/*
	int i = alpha;
	while (i++ < n) {																// uncomment to potentially get information about prime divisors of N
		get_next_RST_i (RST_i, RST_i, tmp_val, QPP, r, N);
	}
	if (mpz_divisible_p(tmp_val[0], N) && mpz_divisible_p(RST_i[1], N)) {				// if N | T_n^2 -4 and S_n
		gmp_printf ("A prime divisor of N = %d*%d^%d + (%d)gamma_n(r) must satisfy p^4 = 1 (mod %d^%d\n", A, r, n, eta, r, alpha);
	}
*/
	mpz_clear (N);
	mpz_clear (rEXPn);
	mpz_clear (gamma_n_r);
	mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
	mpz_clear (RST_i[0]); mpz_clear (RST_i[1]); mpz_clear (RST_i[2]);
	mpz_clear (S_im1);
	return 0;
}

/*			Theorem 7.5 Primality Test from RWG
Parameters:
	A	an even positive integer not divisible by r
	n	a positive integer
	eta	= +/-1

Returns:
	integer		based on N = A*5^n + eta*gamma_n(5)
				-2 if no such number exists (gamma_n(r) cannot be computed (sqrt(-1) doesn't exist for r = 3 (mod 4)) )
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

Runtime:	log(log(N)) * log(N)^2
	from hensel lifting to find N
*/
int primality_test_7_5 (int A, int n, int eta) {
	int r = 5;
	if (A % 2 != 0 || eta*eta != 1 || A % r == 0) {
		return -1;
	}
	mpz_t N; mpz_init (N);
	mpz_t rEXPn; mpz_init (rEXPn);
	mpz_ui_pow_ui (rEXPn, r, n);															// r^n
	if (mpz_cmp_ui (rEXPn, A/2) <= 0) {
		mpz_clear (N);
		mpz_clear (rEXPn);
		return -1;
	}
	mpz_t gamma_n_r; mpz_init (gamma_n_r);
	mpz_t tmp_val[3]; mpz_init (tmp_val[0]); mpz_init (tmp_val[1]); mpz_init (tmp_val[2]);
	mpz_set_si (tmp_val[0], -1);
	mpz_set_ui (tmp_val[1], r);
	if (! h_lift_root (gamma_n_r, tmp_val[0], tmp_val[1], n)) {							// = sqrt(-1) if r^n = 1 (mod 4)
		mpz_clear (N);
		mpz_clear (rEXPn);
		mpz_clear (gamma_n_r);
		mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]); mpz_clear (tmp_val[2]);
		return -2;
	}
	if (mpz_tstbit(gamma_n_r, 0) == 0) {												// if root found is even
		mpz_neg (gamma_n_r, gamma_n_r);
		mpz_add (gamma_n_r, gamma_n_r, rEXPn);											// to get odd sqrt(-1)
	}																					// gamma_n(r) = sqrt(-1) (mod r^n)
	mpz_mul_ui (N, rEXPn, A);															// = Ar^n
	mpz_mul_si (tmp_val[0], gamma_n_r, eta);
	mpz_add (N, N, tmp_val[0]);														// computed N
	mpz_clear (gamma_n_r);
//************************************************************************************************************************************************************************************
//if (mpz_probab_prime_p (N, 20)) {gmp_printf("n = %d gives a likely prime:\n", n);}
//************************************************************************************************************************************************************************************
	int divRes = trial_div(N, 0);
	if (divRes) {													// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
		if (mpz_cmpabs_ui (N, divRes) == 0) {
			mpz_clear(N);
			mpz_clear (rEXPn);
			mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]); mpz_clear (tmp_val[2]);
			return 1;
		}
		else {
			mpz_clear(N);
			mpz_clear (rEXPn);
			mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]); mpz_clear (tmp_val[2]);
			return 0;
		}
	}
	int QPP[3];
	if (!find_QPP_7_5 (QPP, tmp_val[0], N)) {
//gmp_printf("Couldn't find Q, P1, and P2 for n = %d. Not implemented.\n", n);
		return primality_test_7_2_4 (A, r, n, eta);								// slower than this method but no chance of failure on find_QPP
	}
	int Delta = QPP[1]*QPP[1] - 4*QPP[2];										// Delta = P1^2 - 4P2
	mpz_t S_i; mpz_init (S_i);
	mpz_t T_i; mpz_init (T_i);
	get_ST_0 (S_i, T_i, QPP, rEXPn, tmp_val, N);
	mpz_clear (rEXPn);
	for (int i = 0; i < n - 1; i++) {
		get_next_ST_i (S_i, T_i, Delta, tmp_val, N);								// compute S_(i+1) and T_(i+1)
	}																				// have S_(n-1) and T_(n-1)
	mpz_clear (tmp_val[2]);
	mpz_powm_ui (tmp_val[0], S_i, 2, N);
	mpz_mul_si (tmp_val[0], tmp_val[0], 4 * Delta);
	mpz_sub_ui (tmp_val[0], tmp_val[0], 5);										// 4*Delta*S_(n-1)^2 - 5 (mod N)
	mpz_mul_ui (tmp_val[1], T_i, 2);
	mpz_add_ui (tmp_val[1], tmp_val[1], 1);										// 2*T_(n-1) + 1 (mod N)
	mpz_clear (S_i);
	mpz_clear (T_i);
	if (mpz_divisible_p(tmp_val[0], N) && mpz_divisible_p(tmp_val[1], N)) {		// if two above (commented) values are 0 (mod N)
		mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
		mpz_clear (N);
		return 1;
	}
	mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]); 
	mpz_clear (N);
	return 0;
}













