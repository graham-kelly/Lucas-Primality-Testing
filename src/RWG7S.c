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
	long int h = mpz_sizeinbase (m, 2);												// get position of leftmost set bit
	long int i = 1;
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
	}
	mpz_mod (L_m, W_i[0], N);			//should this be W_i[2] and [3] ?
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
*/
_Bool get_cUV_i (mpz_t U_i, mpz_t V_i, int QPP[3], mpz_t i, mpz_t N) {
	long int Delta = QPP[1] * QPP[1] - 4 * QPP[2];							// Delta = P1^2 - 4 * P2
	long int E = (QPP[2] + 4*QPP[0])*(QPP[2] + 4*QPP[0]) - 4 * QPP[0] * QPP[1] * QPP[1];			// E = (P2 + 4Q)^2 - 4 * Q * P1^2
	long int D = E * Delta * Delta * QPP[0] * QPP[0];						// discriminant D = E * Delta^2 * Q^2
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
		long int R = (long int) sqrt(QPP[0]);
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
				gmp_printf("Do need to implement Lehmer's algorithm");
				return 0;										//can be fixed but isn't implemented
				// need actual implementation of Lehmer's algorithm here, use Lehmer 3.1 to compute the u_i and v_i values
			}
		}
		else {
			if (mpz_tstbit(i, 0) == 0) {						// if i is even
				mpz_t constant; mpz_init (constant);
				get_KL_m (U_i, V_i, i_div_2, QPP, N);			// get K_(i/2) and L_(i/2)
				mpz_mul_ui (U_i, U_i, 2);
				mpz_mul_ui (V_i, V_i, 2);
				mpz_set_ui (constant, QPP[0]);
				mpz_powm (constant, constant, i_div_2, N);		// = Q^(i/2)
				mpz_mul (U_i, U_i, constant);					// U_i = 2 * Q^(i/2) * K_(i/2)
				mpz_mul (V_i, V_i, constant);					// V_i = 2 * Q^(i/2) * L_(i/2)
				mpz_clear (constant);
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
	
*/
_Bool get_RST_0 (mpz_t RST[3], mpz_t rEXPn, int QPP[3], mpz_t N) {
	mpz_t constant; mpz_init (constant);
	mpz_t i; mpz_init (i);
	mpz_mul (i, N, N);
	mpz_add_ui (i, i, 1);
	mpz_divexact (i, i, rEXPn);								// i = (N^2 + 1) / (r^n)
	if (!get_cUV_i (RST[0], RST[2], QPP, i, N)) {
		return 0;
	}
	mpz_set_ui (constant, QPP[0]);
	mpz_tdiv_q_2exp (i, i, 1);								// = i/2
	mpz_powm (constant, constant, i, N);					// constant = Q^(i/2)
	mpz_mul_ui (constant, constant, 2);						// = 2*Q^(i/2)
	mpz_invert (constant, constant, N);
	mpz_mul (RST[0], RST[0], constant);
	mpz_mod (RST[0], RST[0], N);
	mpz_set (RST[1], RST[0]);
	mpz_mul (RST[2], RST[2], constant);
	mpz_mod (RST[2], RST[2], N);
	mpz_clear (i);
	mpz_clear (constant);
	return 1;
}

/*			increment sequence term for R, S, T sequences (RWG eq 6.5-7)
Parameters:
	oldRST		values of R_i, S_i, T_i
	constants	memory already allocated for this method (will be called many times in a row)
	QPP			array of (small) integers selected with properties described in sections 3 and 7 of RWG
	r			small odd prime
	N			arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)
	
Returns:
	newRST		values of R_(i+1), S_(i+1), T_(i+1)
*/
void get_next_RST_i (mpz_t newRST[3], mpz_t oldRST[3], mpz_t constants[2], int QPP[3], long int r, mpz_t N) {
	long int delta = QPP[1] * QPP[1] - 4 * QPP[2];				// = P1^2 - 4*P2
	long int k = (r-1)/2;
	mpz_mul (constants[0], oldRST[2], oldRST[2]);
	mpz_mul (constants[1], oldRST[1], oldRST[1]);
	mpz_mul_si (constants[1], constants[1], delta);
	mpz_mod (constants[0], constants[0], N);
	mpz_mod (constants[1], constants[1], N);
	get_H_k (newRST[0], constants[0], constants[1], k, N);		// R_i+1 = H_k(T_i^2, delta*S_i^2) (mod N)
	get_H_k (constants[0], constants[1], constants[0], k, N);
	mpz_mul (newRST[2], oldRST[2], constants[0]);				// T_i+1 = T_i * H_k(delta*S_i^2, T_i^2) (mod N)
	mpz_mul (newRST[1], oldRST[1], newRST[0]);					// S_i+1 = S_i * R_i+1
	mpz_mod (newRST[0], newRST[0], N);
	mpz_mod (newRST[1], newRST[1], N);
	mpz_mod (newRST[2], newRST[2], N);
}


//_Bool new_get_RST_i (mpz_t rop[3], long int i, int QPP[3], long int r, mpz_t rEXPn, mpz_t N);
/*			get particular term of R, S, T sequences (RWG 6.4 - 6.11)
Parameters:
	i			index of sequences
	QPP			array of (small) integers selected with properties described in sections 3 and 7 of RWG
	r			small prime
	rEXPn		r^n
	N			arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)
	
Returns:
	rop			array of mpz_t [R_i, S_i, T_i]
*/
_Bool get_RST_i (mpz_t rop[3], long int i, int QPP[3], long int r, mpz_t rEXPn, mpz_t N) {		//		iterative appraoch
	mpz_t RST[3]; mpz_init (RST[0]); mpz_init (RST[1]); mpz_init (RST[2]);
	mpz_t constants[2]; mpz_init (constants[0]); mpz_init (constants[1]);		// allocate space for get_next_RST_i method
	if (!get_RST_0 (RST, rEXPn, QPP, N)) {										// compute R_0, S_0, T_0
		return 0;
	}
	long int j = 0;
	while (j++ < i) {
		get_next_RST_i (RST, RST, constants, QPP, r, N);						// constants are space reserved in memory to avoid reallocaiton every iteration
	}
	mpz_mod (rop[0], RST[0], N);
	mpz_mod (rop[1], RST[1], N);
	mpz_mod (rop[2], RST[2], N);
	
	

/*
if (new_get_RST_i (RST, i, QPP, r, rEXPn, N)) {
	if (mpz_cmp(RST[0], rop[0]) != 0 || mpz_cmp(RST[1], rop[1]) != 0 || mpz_cmp(RST[2], rop[2]) != 0) {
		gmp_printf ("new R_i = %Zd\nold R_i = %Zd\n", RST[0], rop[0]);
		gmp_printf ("new S_i = %Zd\nold S_i = %Zd\n", RST[1], rop[1]);
		gmp_printf ("new T_i = %Zd\nold T_i = %Zd\n\n", RST[2], rop[2]);
	}
}
*/
	
	
	mpz_clear (RST[0]); mpz_clear (RST[1]); mpz_clear (RST[2]);
	mpz_clear (constants[0]); mpz_clear (constants[1]);
	return 1;
}

/*			get particular term of R, S, T sequences (RWG 6.8 - 6.10)
Parameters:
	i			index of sequences
	QPP			array of (small) integers selected with properties described in sections 3 and 7 of RWG
	r			small prime
	rEXPn		r^n
	N			arithmetic done (mod N), N = Ar^n + eta*gamma_n(r)
	
Returns:
	rop			array of mpz_t [R_i, S_i, T_i]
*/
/*
_Bool new_get_RST_i (mpz_t rop[3], long int i, int QPP[3], long int r, mpz_t rEXPn, mpz_t N) {
gmp_printf("i = %d\nN = %Zd\n", i, N);
	if (i == 0) {
		if (get_RST_0 (rop, rEXPn, QPP, N)) {
			return 1;
		}
		else {
			return 0;
		}
	}
	mpz_t constants[3]; mpz_init (constants[0]); mpz_init (constants[1]); mpz_init (constants[2]);
	mpz_mul (constants[0], N, N);
	mpz_add_ui (constants[0], constants[0], 1);					// N^2 + 1
	mpz_ui_pow_ui (constants[1], r, i);							// r^i
	mpz_divexact (constants[1], rEXPn, constants[1]);			// r^(n-i) == r^n / r^i
	mpz_divexact (constants[0], constants[0], constants[1]);	// N^2 + 1 / r^(n-i)
	if (!get_cUV_i (rop[0], rop[2], QPP, constants[0], N)) {	// rop[0] and rop[2] = U_(N^2 + 1 / r^(n-i)) and V_(N^2 + 1 / r^(n-i)) respectively
		mpz_clear (constants[0]); mpz_clear (constants[1]); mpz_clear (constants[2]);
		return 0;
	}
gmp_printf ("1st term\tU_%Zd = %Zd\n\t\tV_%Zd = %Zd\n", constants[0], constants[1], constants[0], constants[2]);
	mpz_cdiv_q_ui (constants[1], N, 2);							// 1/2
	mpz_mul (constants[2], constants[0], constants[1]);			// = N^2 + 1 / 2r^(n-i)
	mpz_set_si (rop[1], QPP[0]);
	mpz_powm (constants[1], rop[1], constants[2], N);			// Q^(N^2 + 1 / 2r^(n-i))
	mpz_mul_ui (constants[1], constants[1], 2);					// 2Q^(N^2 + 1 / 2r^(n-i))
	if (!mpz_invert (rop[1], constants[1], N)) {				// 1 / 2Q^(N^2 + 1 / 2r^(n-i))
		mpz_clear (constants[0]); mpz_clear (constants[1]); mpz_clear (constants[2]);
		return 0;
	}
//gmp_printf ("\nU_i = %Zd\nV_i = %Zd\n", rop[0], rop[2]);
	mpz_mul (rop[2], rop[1], rop[2]);							// = T_i = V_(N^2 + 1 / r^(n-i)) / (2 (N^2 + 1 / r^(n-i)))
	mpz_mul (rop[1], rop[0], rop[1]);							// = S_i = U_(N^2 + 1 / r^(n-i)) / (2 (N^2 + 1 / r^(n-i)))
//gmp_printf ("subscript = %Zd\n", constants[0]);
	mpz_divexact_ui (constants[0], constants[0], r);			// N^2 + 1 / r^(n-i+1)
//gmp_printf ("subscript = %Zd\n", constants[0]);
	if (!get_cUV_i (constants[1], constants[2], QPP, constants[0], N)) {		// constants[1] = U_(N^2 + 1 / r^(n-i+1)), don't need V_(...) term
		mpz_clear (constants[0]); mpz_clear (constants[1]); mpz_clear (constants[2]);
		return 0;
	}
gmp_printf ("U_%Zd = %Zd\nV_%Zd = %Zd\n", constants[0], constants[1], constants[0], constants[2]);
	if (!mpz_invert (constants[1], constants[1], N)) {
		mpz_clear (constants[0]); mpz_clear (constants[1]); mpz_clear (constants[2]);
		return 0;
	}
gmp_printf ("1/U_i = %Zd\n", constants[1]);
	mpz_mul (rop[0], rop[0], constants[1]);			// = U_(N^2 + 1 / r^(n-i)) / U_(N^2 + 1 / r^(n-i+1))
	mpz_divexact_ui (constants[0], constants[0], 2);			// N^2 + 1 / 2r^(n-i+1)
	mpz_mul_ui (constants[0], constants[0], r-1);				// (r-1) (N^2 + 1 / 2r^(n-i+1))
	mpz_set_si (constants[2], QPP[0]);
	mpz_powm (constants[0], constants[2], constants[0], N);		// = Q^ (r-1) (N^2 + 1) / (2r^(n-i+1))
gmp_printf ("Q term = %Zd\n", constants[0]);
	mpz_mul (constants[0], constants[0], constants[1]);			// U_(N^2 + 1 / r^(n-i+1)) * Q^ (r-1) (N^2 + 1) / (2r^(n-i+1))
	if (!mpz_invert (constants[0], constants[0], N)) {			// 1 / (line above)
		mpz_clear (constants[0]); mpz_clear (constants[1]); mpz_clear (constants[2]);
		return 0;
	}
	mpz_mul (rop[0], rop[0], constants[0]);						// = R_i
	mpz_mod (rop[0], rop[0], N);
	mpz_mod (rop[1], rop[1], N);
	mpz_mod (rop[2], rop[2], N);
	mpz_clear (constants[0]); mpz_clear (constants[1]); mpz_clear (constants[2]);
	return 1;
}

*/




