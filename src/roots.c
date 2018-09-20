#include <gmp.h>
#include "utils.h"

/*			Verrify roots (mod p)
checks that root^2 = k (mod p) for integers root, k and p
*/
_Bool verrify_sqrt (mpz_t root, mpz_t k, mpz_t p) {
	mpz_t constant; mpz_init (constant);
	mpz_mul (constant, root, root);
	if (mpz_congruent_p (constant, k, p)) {
		mpz_clear (constant);
		return 1;
	}
	else {
		mpz_clear (constant);
		return 0;
	}
}

/*			Cipolla's algorithm for finding square roots (mod p)
Parameters:
	n	the integer to find the square root of
	p	the integer to find the square root modulo (if p is not prime result != sqrt(n) (mod p) always, also this may be indicated by false return value)

Return:
	result		= sqrt(n) (mod p)
	bool		if p is not prime then this method may fail in this case return 0 (may also return 1 and wrong result if p is not prime)
	
Efficiency:	Cipolla's if S(S-1)>8m+20 where		S=max(s such that 2^s devides p-1) and m = number of binary digits of (p+1)/2 according to Theorem 3.3 on page 433 of Square Roots Modulo P
Latin American Symposium on Theoretical Informatics
LATIN 2002: LATIN 2002: Theoretical Informatics pp 430-434

Primality test from RWG's Theorem 2.10
	=>	A < 2^n and S = n
	=>	use Cipolla's if	n^2 - 9n - 8logA > 21
							n^2 > 17n + 21
							n > 18
*/
_Bool c_get_sqrt (mpz_t result, mpz_t n, mpz_t p) {
	mpz_t mod; mpz_init(mod);
	mpz_t a; mpz_init(a);
	mpz_t i; mpz_init(i);
	mpz_t base[3]; mpz_init(base[0]); mpz_init(base[1]); mpz_init(base[2]);
	mpz_t exp; mpz_init(exp);
	mpz_t check; mpz_init(check);
	mpz_set_ui (a, 2);												// a will found such that a^2 - n is a quadratic non-residue (mod p)
	mpz_mul (i, a, a);
	mpz_sub (i, i, n);												// i = a^2-n
	while (mpz_jacobi(i, p) != -1 && mpz_cmp(a, p) != 0) {			// find a such that jacobi(a, p) = -1
		mpz_add_ui (a, a, 1);
		mpz_mul (i, a, a);
		mpz_sub (i, i, n);											// a = i^2-n		
	}
	mpz_mul (mod, p, p);											// arithmetic in exponentiation is mod p^2
	mpz_cdiv_q_ui (exp, p, 2);										// exp = (p+1)/2, (note that p is odd)
	mpz_set (base[0], a);
	mpz_set_ui (base[1], 1);
	mpz_mul (base[2], a, a);
	mpz_sub (base[2], base[2], n);									// base = a + 1 * sqrt(a^2-n)
	mpz_mod (base[2], base[2], mod);
	mpzrz_exp (result, check, base, exp, mod);						// binary exponentiation for numbers x + y * sqrt(z), integer part is sqrt(n) (mod p)
	mpz_mod (check, check, p);										// non-integer part of eq should be 0 (mod p)
	if (mpz_cmp_ui (check, 0) != 0) {
		return 0;
	}
	mpz_clear (mod);
	mpz_clear (a);
	mpz_clear (i);	
	mpz_clear (base[0]); mpz_clear (base[1]); mpz_clear (base[2]);
	mpz_clear (exp);
	mpz_clear (check);
	return 1;
}

/*			Tonelli–Shanks algorithm for finding square roots (mod p)
Parameters:
	a	the integer to find the square root of
	q	the integer to find the square root modulo (if q is not prime result != sqrt(n) (mod q) always, also this may be indicated by false return value)

Return:
	result		= sqrt(n) (mod q)
	bool		if q is not prime then this method may fail in this case return 0 (may also return 1 and wrong result if q is not prime)

Will usually be faster than Cipolla's algorithm unless p-1 is highly divisible by 2 (occurs in test based on Theorem 2.8)
This method was addapted from The history of the Theory of Numbers volume 1 p215 by Leonard Eugene Dickson
*/
_Bool ts_get_sqrt (mpz_t result, mpz_t a, mpz_t q) {
	mpz_t constants[3]; mpz_init(constants[0]); mpz_init(constants[1]); mpz_init(constants[2]);
	mpz_t g; mpz_init(g);
	mpz_t e; mpz_init(e);
	mpz_t t; mpz_init(t);
	mpz_set_ui (g, 2);
	while (mpz_jacobi(g, q) != -1 && mpz_cmp(g,q) < 0) {	//find g such that jacobi(g, q) = -1
		mpz_add_ui (g, g, 1);								//try next g;
	}														//got g now
	if (mpz_cmp(g, q) == 0) {
		return 0;
	}
	mpz_sub_ui (constants[0], q, 1);						// = q-1
	long int S = mpz_scan1 (constants[0], 0);				// find first non-zero bit (from right to left) in p-1 and store its bit-number in S
	mpz_set_ui (t, 1);
	mpz_mul_2exp (t, t, S);									// = 2^S
	mpz_divexact (t, constants[0], t);						// t = (q-1) / (2^S)
	mpz_set_ui (e, 0);										// initialize e
	long int i = 2;
	mpz_invert (constants[0], g, q);						// constant[0] = g^-1
	mpz_sub_ui (constants[1], q, 1);
	mpz_tdiv_q_2exp (constants[1], constants[1], i);		// constant[1] = (q-1) / 2^i
	while (i <= S) {
		mpz_powm (constants[2], constants[0], e, q);		// = g^-e
		mpz_mul (constants[2], a, constants[2]);			// = ag^-e
		mpz_powm (constants[2], constants[2], constants[1], q);		// = (ag^-e)^(q-1)/2^i
		if (mpz_cmp_ui (constants[2], 1) != 0) {			//if (ag^-e)^(q-1)/2^i != 1
			mpz_set_ui (constants[2], 1);
			mpz_mul_2exp (constants[2], constants[2], i-1);
			mpz_add (e, e, constants[2]);					// e = 2^(i-1) + e
		}
		i++;
		mpz_tdiv_q_2exp (constants[1], constants[1], 1);	// = (q-1)/(2^i), (each iteration left shifts constants[1] by 1)
	}
	mpz_powm (constants[0], constants[0], e, q);			// = g^-e
	mpz_mul (constants[0], a, constants[0]);				// = h = ag^-e
	mpz_tdiv_q_2exp (constants[1], e, 1);					// = e/2
	mpz_powm (constants[1], g, constants[1], q);			// = g^(e/2)
	mpz_add_ui (constants[2], t, 1);						// = t+1
	mpz_tdiv_q_2exp (constants[2], constants[2], 1);		// = (t+1) / 2
	mpz_powm (result, constants[0], constants[2], q);		// = h^((t+1)/2)
	mpz_mul (result, result, constants[1]);					// = (g^(e/2))h^((t+1)/2)	
	mpz_mod (result, result, q);
	mpz_clear(constants[0]); mpz_clear(constants[1]); mpz_clear(constants[2]);
	mpz_clear(g);
	mpz_clear(e);
	mpz_clear(t);
	return 1;
}

/*			Get sqrt(a) (mod p)
Parameters:
	a	the integer to find the square root of
	p	the integer to find the square root modulo (if p is not prime result != sqrt(n) (mod p) always, also this may be indicated by false return value)

Return:
	result		= smaller sqrt(n) (mod p)
	bool		if p is not prime then this method may fail in this case return 0 otherwise the sqrt is correct and 1 is returned

Uses special methods for when p is 3 (mod 4) or 5 (mod 8)
If p is neither of these then either Cipolla's or Tonelli–Shanks method are used to find the sqrt

*/
_Bool get_sqrt (mpz_t result, mpz_t a, mpz_t p) {
	if (mpz_jacobi(a, p) != 1) {						// if a has no sqrt (mod p) then nothing can be done
		return 0;
	}
	mpz_t constant; mpz_init(constant);
	mpz_t k; mpz_init(k);
	mpz_tdiv_qr_ui (k, constant, p, 4);
	if (mpz_cmp_ui (constant, 3) == 0) {				// if p == 3 (mod 4) then let p = 4k + 3
		mpz_add_ui (constant, k, 1);
		mpz_powm (result, a, constant, p);				// root = a^((k+1)/4) (mod p)
	}
	else {
		mpz_tdiv_qr_ui (k, constant, p, 8);
		if (mpz_cmp_ui (constant, 5) == 0) {			// if p == 5 (mod 8)
			mpz_sub_ui (constant, p, 1);
			mpz_tdiv_q_ui (constant, constant, 4);
			mpz_powm (constant, a , constant, p);		// constant = a^((p-1)/4)
			if (mpz_cmp_ui (constant, 1) == 0) {		// if a^((p-1)/4) = 1
				mpz_add_ui(constant, k, 1);
				mpz_powm (result, a, constant, p);		// root = a^(k+1)
			}
			else {										// if a^((p-1)/4) = -1
				mpz_add_ui(constant, k, 1);
				mpz_powm (result, a, constant, p);		// a^(k+1)
				mpz_mul_ui (constant, k, 2);
				mpz_add_ui (constant, constant, 1);
				mpz_set_ui (k, 2);						// (don't need k anymore, using it as constant)
				mpz_powm (constant, k, constant, p);	// 2^(2k+1)
				mpz_mul (result, result, constant);		// root = 2^(2k+1)*a^(k+1)
			}
		}
		else {
			mpz_set_ui (result, 0);
		}
	}
	if (! verrify_sqrt(result, a, p)) {
		mpz_sub_ui (constant, p, 1);
		int fst_set_bit = mpz_scan1 (constant, 0);				// S = max(s such that 2^s|p-1) 
		mpz_add_ui (constant, p, 1);
		mpz_tdiv_q_2exp (constant, constant, 1);
		int num_bin_dig = mpz_sizeinbase(constant, 2) - 1;		// m = number of binary digits of (p+1)/2
		if (fst_set_bit*(fst_set_bit-1) > 8*num_bin_dig+20) {	// Use Cipolla's method if S(S-1)>8m+20
			if (! c_get_sqrt (result, a, p)) {					// where S=max(s such that 2^s devides p-1) and m = number of binary digits of (p+1)/2
				return 0;
			}
			if (! verrify_sqrt(result, a, p)) {
				return 0;
			}
		}
		else {													// Otherwise use Tonelli–Shanks
			if (! ts_get_sqrt (result, a, p)) {
				return 0;
			}
			if (! verrify_sqrt(result, a, p)) {
				return 0;
			}
		}
	}
	mpz_clear(k);
	mpz_mod (result, result, p);
	mpz_fdiv_q_ui (constant, p, 2);						//to get smaller root
	if (mpz_cmp (constant, result) < 0) {				//if p/2 < result
		mpz_neg (result, result);						//use other root (-origional root)
		mpz_mod (result, result, p);					//put root returned in range (0,p-1)
	}
	mpz_clear(constant);
	return 1;
}

/*			Get Sqrt(a) (mod p^pow)
Parameters:
	a	the integer to find the square root of
	p	the integer to find the square root modulo (if p is not prime result != sqrt(n) (mod p) always, also this may be indicated by false return value)
	power	the positive integer to hensel lift the sqrt to

Return:
	result		= sqrt(n) (mod p^pow)
	bool		if p is not prime then this method may fail in this case return 0 otherwise the sqrt is correct and 1 is returned

Uses hensel lifting and the get_sqrt method above to find the sqrt(a) (mod p^pow)
*/
_Bool h_lift_root (mpz_t rop, mpz_t a, mpz_t p, long int power) {
	mpz_t root; mpz_init (root);
	if (!get_sqrt(root, a, p)) {				// get sqrt(a) (mod p)
		mpz_clear (root);
		return 0;
	}
	int i = 1;
	mpz_t prev_res; mpz_init (prev_res);
	mpz_t constant; mpz_init (constant);
	mpz_t mod; mpz_init (mod);
	mpz_set (mod, p);
	mpz_set (prev_res, root);
	while (i < power) {
		mpz_mul (mod, mod, mod);					// current power = p^i
		mpz_invert (constant, root, mod);
		mpz_mul (root, root, root);
		mpz_neg (root, root);
		mpz_add (root, root, a);					// = a - x_i^2
		mpz_mul (root, root, constant);				// = (a - x_i^2)/x_i
		mpz_cdiv_q_2exp (constant, mod, 1);			// = 1/2
		mpz_mul (root, root, constant);				// = (a - x_i^2)/2x_i
		i <<= 1;									// i *= 2
		mpz_mod (root, root, mod);
		mpz_add (root, root, prev_res);				// + prevoius root
		mpz_set (prev_res, root);
	}
	mpz_pow_ui (mod, p, power);						// = p^power
	mpz_mod (rop, root, mod);						// = sqrt(a) (mod p^power)
	mpz_clear (root);
	mpz_clear (prev_res);
	mpz_clear (constant);
	mpz_clear (mod);
	return 1;
}















