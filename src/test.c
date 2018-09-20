#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <gmp.h>

_Bool h_lift_root (mpz_t rop, mpz_t a, mpz_t p, long int power);
_Bool verrify_root (mpz_t root, mpz_t k, mpz_t p);
_Bool c_get_root (mpz_t result, mpz_t n, mpz_t p);
_Bool ts_get_root (mpz_t result, mpz_t a, mpz_t q);
_Bool get_root (mpz_t result, mpz_t a, mpz_t p);
void mpzrz_mul (mpz_t z[], mpz_t x[], mpz_t y[], mpz_t p);
void mpzrz_sqr (mpz_t z[], mpz_t x[], mpz_t p);
void mpzrz_exp (mpz_t int_result, mpz_t root_result, mpz_t base[3], mpz_t exp, mpz_t mod);

int main(int argc, char *argv[]) {
	int r = 13;								// must be 1 (mod 4) !!!
	int set_n[] = {1, 30, 32, 33, 40, 50, 65, 66, 83, 93};
	int i = 0;
	mpz_t constants[2]; mpz_init (constants[0]); mpz_init (constants[1]);
	mpz_t gamma_n_r; mpz_init (gamma_n_r);
	mpz_t rEXPn; mpz_init (rEXPn);	
	
	mpz_set_si (constants[0], -1);
	mpz_set_ui (constants[1], r);
	
	while (i < 10) {
		h_lift_root (gamma_n_r, constants[0], constants[1], set_n[i]);							// = sqrt(-1) if r^n = 1 (mod 4)
		if (mpz_tstbit(gamma_n_r, 0) == 0) {												// if root found is even
			mpz_neg (gamma_n_r, gamma_n_r);
			mpz_add (gamma_n_r, gamma_n_r, rEXPn);											// to get odd sqrt(-1)
		}																					// gamma_n(r) = sqrt(-1) (mod r^n)
		gmp_printf ("gamma_%d(%d) = %Zd\n", set_n[i++], r, gamma_n_r);
	}
	mpz_clear (constants[0]); mpz_clear (constants[1]);
	mpz_clear (gamma_n_r);
	mpz_clear (rEXPn);
	
	return 0;
}

_Bool h_lift_root (mpz_t rop, mpz_t a, mpz_t p, long int power) {
	mpz_t root; mpz_init (root);
	if (!get_root(root, a, p)) {
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
		mpz_add (root, root, a);				// = a - x_i^2
		mpz_mul (root, root, constant);			// = (a - x_i^2)/x_i
		mpz_cdiv_q_2exp (constant, mod, 1);			// = 1/2
		mpz_mul (root, root, constant);			// = (a - x_i^2)/2x_i
		i <<= 1;
		mpz_mod (root, root, mod);
		mpz_add (root, root, prev_res);			// + prevoius root
		mpz_set (prev_res, root);
	}
	mpz_pow_ui (mod, p, power);						// = p^power
	mpz_mod (rop, root, mod);
	mpz_clear (root);
	mpz_clear (prev_res);
	mpz_clear (constant);
	mpz_clear (mod);
	return 1;
}

_Bool verrify_root (mpz_t root, mpz_t k, mpz_t p) {
	mpz_t constant; mpz_init (constant);
	mpz_mul (constant, root, root);
	mpz_sub (constant, constant, k);
	mpz_mod (constant, constant, p);
	if (mpz_cmp_ui (constant, 0) == 0) {
		mpz_clear (constant);
		return 1;
	}
	else {
		mpz_clear (constant);
		return 0;
	}
}

//Cipolla's algorithm for finding square roots mod p (better for k>18 where p = A * 2^k +/- 1 and A < 2^k)
_Bool c_get_root (mpz_t result, mpz_t n, mpz_t p) {
	mpz_t mod; mpz_init(mod);
	mpz_t a; mpz_init(a);
	mpz_t i; mpz_init(i);
	mpz_t base[3]; mpz_init(base[0]); mpz_init(base[1]); mpz_init(base[2]);
	mpz_t exp; mpz_init(exp);
	mpz_t check; mpz_init(check);
	mpz_set_ui (a, 2);
	mpz_mul (i, a, a);
	mpz_sub (i, i, n);												// i = a^2-n
	while (mpz_jacobi(i, p) != -1 && mpz_cmp(a, p) != 0) {			// find a such that jacobi(a, p) = -1
		mpz_add_ui (a, a, 1);
		mpz_mul (i, a, a);
		mpz_sub (i, i, n);											// a = i^2-n		
	}
	mpz_mul (mod, p, p);											// arithmetic in exponentiation is mod p^2
	mpz_add_ui (exp, p, 1);
	mpz_tdiv_q_ui (exp, exp, 2);									// exp = (p+1)/2
	mpz_set (base[0], a);
	mpz_set_ui (base[1], 1);
	mpz_mul (base[2], a, a);
	mpz_sub (base[2], base[2], n);									// base = a+root(a^2-n)
	mpz_mod (base[2], base[2], mod);
	mpzrz_exp (result, check, base, exp, mod);						// binary exponentiation for numbers x+y*root(z)
	mpz_mod (check, check, p);
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

_Bool ts_get_root (mpz_t result, mpz_t a, mpz_t q) {
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
	mpz_set_ui(e, 0);
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
		mpz_tdiv_q_2exp (constants[1], constants[1], 1);	// = (q-1)/(2^i)
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

//Get the square root (mod p) of a and store in result
_Bool get_root (mpz_t result, mpz_t a, mpz_t p) {
	if (mpz_jacobi(a, p) != 1) {
		return 0;
	}
	mpz_t constant; mpz_init(constant);
	mpz_t k; mpz_init(k);
	mpz_tdiv_qr_ui (k, constant, p, 4);
	if (mpz_cmp_ui (constant, 3) == 0) {				//if p == 3 (mod 4) then let p = 4k + 3
		mpz_add_ui (constant, k, 1);
		mpz_powm (result, a, constant, p);				//root = a^((k+1)/4) mod p
	}
	else {
		mpz_tdiv_qr_ui (k, constant, p, 8);
		if (mpz_cmp_ui (constant, 5) == 0) {			//if p == 5 (mod 8)
			mpz_sub_ui (constant, p, 1);
			mpz_tdiv_q_ui (constant, constant, 4);
			mpz_powm (constant, a , constant, p);		//constant = a^((p-1)/4)
			if (mpz_cmp_ui (constant, 1) == 0) {		//if a^((p-1)/4) = 1
				mpz_add_ui(constant, k, 1);
				mpz_powm (result, a, constant, p);		//root = a^(k+1)
			}
			else {										//if a^((p-1)/4) = -1
				mpz_add_ui(constant, k, 1);
				mpz_powm (result, a, constant, p);		//a^(k+1)
				mpz_mul_ui (constant, k, 2);
				mpz_add_ui (constant, constant, 1);
				mpz_set_ui (k, 2);						//(don't need k anymore, using it as constant)
				mpz_powm (constant, k, constant, p);	//2^(2k+1)
				mpz_mul (result, result, constant);		//root = 2^(2k+1)*a^(k+1)
			}
		}
		else {
			mpz_set_ui (result, 0);
		}
	}
	if (! verrify_root(result, a, p)) {
		mpz_sub_ui (constant, p, 1);
		int fst_set_bit = mpz_scan1 (constant, 0);				// S = max(s such that 2^s|p-1) 
		mpz_add_ui (constant, p, 1);
		mpz_tdiv_q_2exp (constant, constant, 1);
		int num_bin_dig = mpz_sizeinbase(constant, 2) - 1;		// m = number of binary digits of (p+1)/2
		if (fst_set_bit*(fst_set_bit-1) > 8*num_bin_dig+20) {	// Use Cipolla's method if S(S-1)>8m+20
			if (! c_get_root (result, a, p)) {					// where S=max(s such that 2^s devides p-1) and m = number of binary digits of (p+1)/2
				return 0;
			}
			if (! verrify_root(result, a, p)) {
				return 0;
			}
		}
		else {													// Otherwise use Tonelli–Shanks
			if (! ts_get_root (result, a, p)) {
				return 0;
			}
			if (! verrify_root(result, a, p)) {
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

void mpzrz_mul (mpz_t z[], mpz_t x[], mpz_t y[], mpz_t p) {
	mpz_t constant1; mpz_init(constant1);
	mpz_t constant2; mpz_init(constant2);
	mpz_t result; mpz_init(result);
	mpz_mul (constant1, x[0], y[0]);
	mpz_mul (constant2, x[1], y[1]);
	mpz_mul (constant2, constant2, x[2]);
	mpz_add (result, constant1, constant2);
	mpz_mul (constant1, x[0], y[1]);
	mpz_mul (constant2, x[1], y[0]);
	
	mpz_set (z[0], result);
	mpz_add (z[1], constant1, constant2);
	mpz_set (z[2], x[2]);
	
	mpz_clear(constant1);
	mpz_clear(constant2);
	mpz_clear(result);
	
	mpz_mod (z[0], z[0], p);
	mpz_mod (z[1], z[1], p);
}

//square for numbers a + b * root(c)
void mpzrz_sqr (mpz_t z[], mpz_t x[], mpz_t p) {
	mpz_t constant; mpz_init(constant);
	mpz_t result; mpz_init(result);
	
	mpz_mul (constant, x[0], x[0]);
	mpz_mul (result, x[1], x[1]);
	mpz_mul (result, result, x[2]);
	mpz_add (result, constant, result);
	
	mpz_mul (constant, x[0], x[1]);
	
	mpz_set (z[0], result);
	mpz_mul_ui (z[1], constant, 2);
	mpz_set (z[2], x[2]);
	
	mpz_clear(constant);
	mpz_clear(result);
	
	mpz_mod (z[0], z[0], p);
	mpz_mod (z[1], z[1], p);
}

//binary exponentiation algorithm used in Cipolla's algorithm and u_i/v_i (for numbers a + b * root(c))
void mpzrz_exp (mpz_t int_result, mpz_t root_result, mpz_t base[3], mpz_t exp, mpz_t mod) {
	mpz_t cexp; mpz_init(cexp);
	mpz_t term[3]; mpz_init(term[0]); mpz_init(term[1]); mpz_init(term[2]);
	mpz_t product[3]; mpz_init(product[0]); mpz_init(product[1]); mpz_init(product[2]);
	mpz_set (cexp, exp);									//need a copy of the exponent to work with
	mpz_set (term[0], base[0]);								//initialize term
	mpz_set (term[1], base[1]);
	mpz_set (term[2], base[2]);
	while (! mpz_tstbit (cexp, 0)) {						//before first set bit of exponent
		mpzrz_sqr (term, term, mod);						//square the term to multiply by
		mpz_tdiv_q_2exp (cexp, cexp, 1);					//right shift exponent by 1
	}														//now at first set bit of exponent => can initialize product
	mpz_set (product[0], term[0]);							//initialize product
	mpz_set (product[1], term[1]);
	mpz_set (product[2], term[2]);
	mpz_tdiv_q_2exp (cexp, cexp, 1);						//right shift exponent by 1
	while (mpz_cmp_ui (cexp, 0) != 0) {						//while exponent has more bits
		mpzrz_sqr (term, term, mod);						//square the term to multiply by
		if (mpz_tstbit (cexp, 0)) {							//if rightmost bit is 1
			mpzrz_mul (product, product, term, mod);		//do multiplication
		}
		mpz_tdiv_q_2exp (cexp, cexp, 1);					//right shift exponent by 1
	}
	mpz_set (int_result, product[0]);
	mpz_set (root_result, product[1]);
	mpz_clear(cexp);
	mpz_clear(term[0]); mpz_clear(term[1]); mpz_clear(term[2]);
	mpz_clear(product[0]); mpz_clear(product[1]); mpz_clear(product[2]);
}

