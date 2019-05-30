#include <gmp.h>
#include <math.h>
#include "utils.h"
#include "RWG2S.h"
#include "RWG7.h"
#include "RWG7S.h"

/*			get K and L sequences from RWG section 4
Parameters:
	m		index of which sequence item is required
	QPP		array of (small) integers selected with properties described in sections 3 and 7 of RWG
	N		arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)
	
Returns:
	K_m		mth term of the K sequence from section 4 of RWG
	L_m		mth term of the L sequence from section 4 of RWG
Runtime:	O(log(N)^2)
	inversion (mod N)
*/
void get_KL_m (mpz_t K_m, mpz_t L_m, mpz_t m, int QPP[3], mpz_t N) {
	mpz_t L_1; mpz_init (L_1);
	mpz_t K_1; mpz_init (K_1);
	mpz_t W_i[4]; mpz_init (W_i[0]); mpz_init (W_i[1]); mpz_init (W_i[2]); mpz_init (W_i[3]);
	mpz_t A; mpz_init (A);
	mpz_t B; mpz_init (B);
	mpz_t C; mpz_init (C);
	mpz_t D; mpz_init (D);
	mpz_set_si (K_1, 2*QPP[0]);
	mpz_invert (K_1, K_1, N);
	mpz_mul_si (L_1, K_1, QPP[1]*QPP[1] - 2*QPP[2] - 4*QPP[0]);						// = L_1 = V_2/2Q
	mpz_mul_si (K_1, K_1, QPP[1]);													// = K_1 = U_2/2Q
	int Delta = QPP[1] * QPP[1] - 4 * QPP[2];
	mpz_set (W_i[0], L_1);															// W_0 = {L_1, K_1, L_2, K_2} where L_2 = L_1^2 + (P1^2 - 4P2)K_1^2 - 2, K_2 = 2K_1L_1
	mpz_set (W_i[1], K_1);
	mpz_mul (W_i[3], L_1, L_1);														// temporary (used to compute L_2)
	mpz_mul (W_i[2], K_1, K_1);
	mpz_mul_si (W_i[2], W_i[2], Delta);
	mpz_add (W_i[2], W_i[2], W_i[3]);
	mpz_sub_ui (W_i[2], W_i[2], 2);													// = L_2
	mpz_mul (W_i[3], K_1, L_1);
	mpz_mul_ui (W_i[3], W_i[3], 2);													// = K_2
	int h = mpz_sizeinbase (m, 2);												// get position of leftmost set bit
	int i = 1;
mpz_mod (W_i[0], W_i[0], N); mpz_mod (W_i[1], W_i[1], N); mpz_mod (W_i[2], W_i[2], N); mpz_mod (W_i[3], W_i[3], N); 
	while (i++ < h) {																	// for each bit of m starting from the left
		mpz_mod (A, W_i[0], N);
		mpz_mod (B, W_i[1], N);
		mpz_mod (C, W_i[2], N);
		mpz_mod (D, W_i[3], N);
		if (!mpz_tstbit (m, h - i)) {												// if bit h-i of m is 0
			mpz_mul (W_i[0], A, A);
			mpz_mul (W_i[1], B, B);													// will use W_i[1] as a temporary value while computation is going on
			mpz_mul_si (W_i[1], W_i[1], Delta);
			mpz_add (W_i[0], W_i[0], W_i[1]);
			mpz_sub_ui (W_i[0], W_i[0], 2);											// = W_i[0]
			mpz_mul (W_i[2], A, C);
			mpz_mul (W_i[1], B, D);
			mpz_mul_si (W_i[1], W_i[1], Delta);
			mpz_add (W_i[2], W_i[2], W_i[1]);
			mpz_sub (W_i[2], W_i[2], L_1);											// = W_i[2]
			mpz_mul (W_i[3], B, C);
			mpz_mul (W_i[1], A, D);
			mpz_add (W_i[3], W_i[3], W_i[1]);
			mpz_sub (W_i[3], W_i[3], K_1);											// = W_i[3]
			mpz_mul (W_i[1], A, B);
			mpz_mul_ui (W_i[1], W_i[1], 2);											// = W_i[1]
		}
		else {																		// if bit h-i of m is 1
			mpz_mul (W_i[0], A, C);
			mpz_mul (W_i[3], B, D);													// will use W_i[3] as a temporary value while computation is going on
			mpz_mul_si (W_i[3], W_i[3], Delta);
			mpz_add (W_i[0], W_i[0], W_i[3]);
			mpz_sub (W_i[0], W_i[0], L_1);											// = W_i[0]
			mpz_mul (W_i[1], B, C);
			mpz_mul (W_i[3], A, D);
			mpz_add (W_i[1], W_i[1], W_i[3]);
			mpz_sub (W_i[1], W_i[1], K_1);											// = W_i[1]
			mpz_mul (W_i[2], C, C);
			mpz_mul (W_i[3], D, D);
			mpz_mul_si (W_i[3], W_i[3], Delta);
			mpz_add (W_i[2], W_i[2], W_i[3]);
			mpz_sub_ui (W_i[2], W_i[2], 2);											// = W_i[2]
			mpz_mul (W_i[3], C, D);
			mpz_mul_ui (W_i[3], W_i[3], 2);											// = W_i[3]
		}
mpz_mod (W_i[0], W_i[0], N); mpz_mod (W_i[1], W_i[1], N); mpz_mod (W_i[2], W_i[2], N); mpz_mod (W_i[3], W_i[3], N); 
	}
	mpz_mod (L_m, W_i[0], N);
	mpz_mod (K_m, W_i[1], N);
	mpz_clear (L_1);
	mpz_clear (K_1);
	mpz_clear (W_i[0]); mpz_clear (W_i[1]); mpz_clear (W_i[2]); mpz_clear (W_i[3]); 
	mpz_clear (A);
	mpz_clear (B);
	mpz_clear (C);
	mpz_clear (D);
}

/*			Get the order 4 lucas like sequences defined in section 3 of RWG
Parameters:
	QPP		array of (small) integers selected with properties described in sections 3 and 7 of RWG
	i		even positive integer sequence index 
	N		arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)

Returns:
	U_i		order 4 lucas like sequence defined in section 3 of RWG
	V_i		order 4 lucas like sequence defined in section 3 of RWG
Runtime:	O(log(N)^2)
	2 inversions (mod N), one in get_uv_i, the other in either invert, or get_KL_m
*/
_Bool get_cUV_i (mpz_t U_i, mpz_t V_i, int QPP[3], mpz_t i, mpz_t N) {
	int Delta = QPP[1] * QPP[1] - 4 * QPP[2];							// Delta = P1^2 - 4 * P2
	int E = (QPP[2] + 4*QPP[0])*(QPP[2] + 4*QPP[0]) - 4 * QPP[0] * QPP[1] * QPP[1];			// E = (P2 + 4Q)^2 - 4 * Q * P1^2
//	int D = E * Delta * Delta * QPP[0] * QPP[0];						// discriminant D = E * Delta^2 * Q^2
	if (Delta == 0) {
		if (!get_uv_i (U_i, V_i, i, QPP[1]/2, QPP[1]*QPP[1]/4 - 4*QPP[0], N)) {
			return 0;
		}
		mpz_mul (U_i, U_i, i);									// U_i = n * u_i(P1/2, Q)
		mpz_mul_ui (V_i, V_i, 2);								// V_i = 2 * v_i(P1/2, Q)
	}
	else {
		mpz_t i_div_2; mpz_init (i_div_2);
		mpz_tdiv_q_2exp (i_div_2, i, 1);
		int R = (int) sqrt(QPP[0]);
		if (E == 0 && Delta != 0) {								// if this is true then Q is a perfect square
			if (mpz_tstbit(i, 0) == 0) {						// if i is even
				mpz_t xcy[3]; mpz_init (xcy[0]); mpz_init (xcy[1]); mpz_init (xcy[2]);
				mpz_t a; mpz_init (a);							// let (x + sqrt(y))^i/2 = a + b sqrt(y)
				mpz_t b; mpz_init (b);
				mpz_set_si (xcy[0], QPP[1] - 2*R);				// where x = P1 - 2R
				mpz_set_ui (xcy[1], 1);
				mpz_set_si (xcy[2], QPP[1]*QPP[1] - 4*QPP[1]*R);// and y = P1^2 - 4*P1*R
				mpzrz_exp (a, b, xcy, i_div_2, N);				// compute (x + sqrt(y))^i/2
				mpz_mul (U_i, b, b);
				mpz_mul (U_i, U_i, xcy[2]);
				mpz_mul_ui (U_i, U_i, 4);
				mpz_invert (xcy[1], xcy[0], N);
				mpz_mul (U_i, U_i, xcy[1]);						// U_i = 4(b^2)y/x
				mpz_mod (U_i, U_i, N);
				mpz_mul (V_i, a, a);
				mpz_mul_2exp (V_i, V_i, 2);						// V_i = 4a^2
				mpz_mod (V_i, V_i, N);
				mpz_clear (xcy[0]); mpz_clear (xcy[1]); mpz_clear (xcy[2]);
				mpz_clear (a);
				mpz_clear (b);
			}
			else {
				gmp_printf("Do need to implement Lehmer's algorithm\n");
				return 0;										//can be fixed but isn't implemented
				// need actual implementation of Lehmer's algorithm here, use Lehmer 3.1 to compute the u_i and v_i values
			}
		}
		else {
			if (mpz_tstbit(i, 0) == 0) {						// if i is even
				mpz_t tmp_val; mpz_init (tmp_val);
				get_KL_m (U_i, V_i, i_div_2, QPP, N);			// get K_(i/2) and L_(i/2)
				mpz_mul_ui (U_i, U_i, 2);
				mpz_mul_ui (V_i, V_i, 2);
				mpz_set_ui (tmp_val, QPP[0]);
				mpz_powm (tmp_val, tmp_val, i_div_2, N);		// = Q^(i/2)
				mpz_mul (U_i, U_i, tmp_val);					// U_i = 2 * Q^(i/2) * K_(i/2)
				mpz_mul (V_i, V_i, tmp_val);					// V_i = 2 * Q^(i/2) * L_(i/2)
				mpz_clear (tmp_val);
			}
			else {
				gmp_printf ("Odd case for Delta != 0\n");		// not sure what to do here
				return 0;										// something has probably gone wrong
			}
		}
		mpz_clear (i_div_2);
	}
	mpz_mod (U_i, U_i, N);
	mpz_mod (V_i, V_i, N);
	return 1;
}

/*			get first term of R, S, T sequences defined in eq 6.4 of RWG
Parameters:
	A			positive even integer
	rEXPn		r^n for small prime r, positive integer n
	gamma_n_r	sqrt(-1) (mod r^n)
	eta			= +/- 1
	QPP			array of (small) integers selected with properties described in sections 3 and 7 of RWG
	N			arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)

Returns:
	R0, S0, T0	as defined in eq 6.4 of RWG
	_Bool		0 if inversion fails (N not prime), else 1
	
Runtime:	O(log(N)^2)
	2 inversions (1 in get_cUV_i)
*/
_Bool get_RST_0 (mpz_t RST[3], int A, mpz_t rEXPn, mpz_t gamma_n_r, int eta, int QPP[3], mpz_t N) {
	mpz_t tmp_val; mpz_init (tmp_val);
	mpz_t i; mpz_init (i);
	mpz_mul (i, N, N);
	mpz_add_ui (i, i, 1);
	mpz_divexact (i, i, rEXPn);								// i = (N^2 + 1) / (r^n)
	if (!get_cUV_i (RST[0], RST[2], QPP, i, N)) {			// inversion failed, N not prime
		return 0;
	}
	mpz_set_ui (tmp_val, QPP[0]);
	mpz_tdiv_q_2exp (i, i, 1);								// = i/2
	mpz_powm (tmp_val, tmp_val, i, N);						// tmp_val = Q^(i/2)
	mpz_mul_ui (tmp_val, tmp_val, 2);						// = 2*Q^(i/2)
	if (!mpz_invert (tmp_val, tmp_val, N)) {
		return 0;
	}
	mpz_mul (RST[0], RST[0], tmp_val);
	mpz_mod (RST[0], RST[0], N);
	mpz_set (RST[1], RST[0]);
	mpz_mul (RST[2], RST[2], tmp_val);
	mpz_mod (RST[2], RST[2], N);
	mpz_clear (i);
	mpz_clear (tmp_val);
	return 1;
}

/*			increment sequence term for R, S, T sequences (RWG eq 6.5-7)
Parameters:
	oldRST		values of R_i, S_i, T_i
	tmp_vals	memory already allocated for this method (will be called many times in a row)
	QPP			array of (small) integers selected with properties described in sections 3 and 7 of RWG
	r			small odd prime
	N			arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)
	
Returns:
	newRST		values of R_(i+1), S_(i+1), T_(i+1)
Runtime:		log(r) * r^2 * M(N)
	due to get_H_k
*/
void get_next_RST_i (mpz_t newRST[3], mpz_t oldRST[3], mpz_t tmp_val[2], int QPP[3], int r, mpz_t N) {
	int delta = QPP[1] * QPP[1] - 4 * QPP[2];				// = P1^2 - 4*P2
	int k = (r-1)/2;
	mpz_mul (tmp_val[0], oldRST[2], oldRST[2]);
	mpz_mul (tmp_val[1], oldRST[1], oldRST[1]);
	mpz_mul_si (tmp_val[1], tmp_val[1], delta);
	mpz_mod (tmp_val[0], tmp_val[0], N);
	mpz_mod (tmp_val[1], tmp_val[1], N);
	get_HI_k (newRST[0], tmp_val[0], tmp_val[0], tmp_val[1], k, N);		// R_i+1 = H_k(T_i^2, delta*S_i^2) (mod N)
	mpz_mul (newRST[2], oldRST[2], tmp_val[0]);							// T_i+1 = T_i * H_k(delta*S_i^2, T_i^2) (mod N)
	mpz_mul (newRST[1], oldRST[1], newRST[0]);							// S_i+1 = S_i * R_i+1
	mpz_mod (newRST[0], newRST[0], N);
	mpz_mod (newRST[1], newRST[1], N);
	mpz_mod (newRST[2], newRST[2], N);
}

/*			get particular term of R, S, T sequences (RWG 6.4 - 6.11)
eq 6.8 - 6.10 of RWG give explicit formulas for these but computing these sequentially is faster (no inversion)

Parameters:
	i			index of sequences
	QPP			array of (small) integers selected with properties described in sections 3 and 7 of RWG
	A			even integer r !| A
	r			small prime
	rEXPn		r^n
	gamma_n_r	sqrt(-1) (mod r^n)
	eta			= +/- 1
	N			arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)
	
Returns:
	rop			array of mpz_t [R_i, S_i, T_i]
	
Runtime:	O(log(N)^2)
	2 inversions in get_RST_0
*/
_Bool get_RST_i (mpz_t rop[3], int i, int QPP[3], int A, int r, mpz_t rEXPn, mpz_t gamma_n_r, int eta, mpz_t N) {
//		************************************		would probably be more efficient to implement eq 6.8 - 6.10 of RWG here			************************************
	mpz_t RST[3]; mpz_init (RST[0]); mpz_init (RST[1]); mpz_init (RST[2]);
	mpz_t tmp_val[2]; mpz_init (tmp_val[0]); mpz_init (tmp_val[1]);		// allocate space for get_next_RST_i method
	if (!get_RST_0 (RST, A, rEXPn, gamma_n_r, eta, QPP, N)) {					// compute R_0, S_0, T_0
		return 0;
	}
	int j = 0;
	while (j++ < i) {
		get_next_RST_i (RST, RST, tmp_val, QPP, r, N);						// tmp_val are space reserved in memory to avoid reallocaiton every iteration
	}
	mpz_mod (rop[0], RST[0], N);
	mpz_mod (rop[1], RST[1], N);
	mpz_mod (rop[2], RST[2], N);
	mpz_clear (RST[0]); mpz_clear (RST[1]); mpz_clear (RST[2]);
	mpz_clear (tmp_val[0]); mpz_clear (tmp_val[1]);
	return 1;
}

/*			get first term of S, T sequences defined in Theorem 7.5 of RWG
Parameters:
	QPP			array of integers as required in Theorem 7.5 of RWG
	rEXPn		5^n for N = A*5^n+eta*gamma_5(n)
	tmp_val	memory already allocated
	N			integer to be shown prime (or not)

Returns:
	S_i			S_0 as defined in Theorem 7.5 of RWG
	T_i			T_0 as defined in Theorem 7.5 of RWG
	_Bool		will be 0 if N is not prime (non-invertable member of Z_N found) else 1

Runtime:	O(log(N)^2)
	inversion
*/
_Bool get_ST_0 (mpz_t S_i, mpz_t T_i, int QPP[3], mpz_t rEXPn, mpz_t tmp_val[2], mpz_t N) {
	mpz_mul (tmp_val[0], N, N);
	mpz_add_ui (tmp_val[0], tmp_val[0], 1);
	mpz_divexact (tmp_val[0], tmp_val[0], rEXPn);								// = (N^2+1) / 5^n
	get_cUV_i (S_i, T_i, QPP, tmp_val[0], N);									// get U_((N^2+1) / 5^n) and V_((N^2+1) / 5^n)
	mpz_divexact_ui (tmp_val[0], tmp_val[0], 2);								// = (N^2+1) / 2*5^n

	mpz_set_ui (tmp_val[1], QPP[0]);
	mpz_powm (tmp_val[0], tmp_val[1], tmp_val[0], N);							// = Q^(N^2+1) / 2*5^n
	mpz_mul_ui (tmp_val[0], tmp_val[0], 2);										// = 2Q^(N^2+1) / 2*5^n
	if (!mpz_invert (tmp_val[0], tmp_val[0], N)) {								// = 1 / 2Q^(N^2+1) / 2*5^n
		return 0;																// not prime if inversion is not possible
	}
	mpz_mul (S_i, S_i, tmp_val[0]);												// got S_0
	mpz_mul (T_i, T_i, tmp_val[0]);												// got T_0
	mpz_mod (S_i, S_i, N);
	mpz_mod (T_i, T_i, N);
	return 1;
}

/*			get first term of S, T sequences defined in Theorem 7.5 of RWG
Parameters:
	Delta		as defined in Theorem 7.5 of RWG
	tmp_val	memory already allocated
	N			integer to be shown prime

Returns:
	S_i			S_i+1 as defined in Theorem 7.5 of RWG
	T_i			T_i+1 as defined in Theorem 7.5 of RWG
Runtime:	O(M(N))
	small exponents (<5) and multiplication (mod N)
*/
void get_next_ST_i (mpz_t S_i, mpz_t T_i, int Delta, mpz_t tmp_val[3], mpz_t N) {
	mpz_powm_ui (tmp_val[0], S_i, 4, N);
	mpz_mul_si (tmp_val[1], tmp_val[0], Delta * Delta);		// 1st term S_i+1 (D^2 * S^4)
	mpz_mul_ui (tmp_val[2], tmp_val[1], 5);					// 1st term T_i+1 (5 * D^2 * S^4)

	mpz_powm_ui (tmp_val[0], S_i, 2, N);
	mpz_mul (tmp_val[0], tmp_val[0], T_i);
	mpz_mul (tmp_val[0], tmp_val[0], T_i);
	mpz_mul_si (tmp_val[0], tmp_val[0], 10 * Delta);		// 2nd term S_i+1 and T_i+1 (10 * D * S^2 * T^2)
	mpz_add (tmp_val[1], tmp_val[1], tmp_val[0]);
	mpz_add (tmp_val[2], tmp_val[2], tmp_val[0]);
	
	mpz_powm_ui (tmp_val[0], T_i, 4, N);					// 3rd term T_i+1 (T^4)
	mpz_add (tmp_val[2], tmp_val[2], tmp_val[0]);
	mpz_mul_ui (tmp_val[0], tmp_val[0], 5);					// 3rd term S_i+1 (5 * T^4)
	mpz_add (tmp_val[1], tmp_val[1], tmp_val[0]);
	
	mpz_powm_ui (tmp_val[0], S_i, 2, N);
	mpz_mul_si (tmp_val[0], tmp_val[0], -5 * Delta);		// 4th term S_i+1 (-5 * D * S^2)
	mpz_add (tmp_val[1], tmp_val[1], tmp_val[0]);
	mpz_mul_si (tmp_val[0], tmp_val[0], 3);					// 4th term T_i+1 (-15 * D * S^2)
	mpz_add (tmp_val[2], tmp_val[2], tmp_val[0]);
	
	mpz_powm_ui (tmp_val[0], T_i, 2, N);
	mpz_mul_si (tmp_val[0], tmp_val[0], -5);				// 5th term T_i+1 (-5 * T^2)
	mpz_add (tmp_val[2], tmp_val[2], tmp_val[0]);
	mpz_mul_si (tmp_val[0], tmp_val[0], 3);					// 5th term S_i+1 (-15 * T^2)
	mpz_add (tmp_val[1], tmp_val[1], tmp_val[0]);
	
	mpz_add_ui (tmp_val[1], tmp_val[1], 5);					// 6th term S_i+1 (5)
	mpz_add_ui (tmp_val[2], tmp_val[2], 5);					// 6th term T_i+1 (5)
	
	mpz_mul (tmp_val[1], tmp_val[1], S_i);
	mpz_mod (S_i, tmp_val[1], N);
	mpz_mul (tmp_val[2], tmp_val[2], T_i);
	mpz_mod (T_i, tmp_val[2], N);
}














