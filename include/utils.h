#ifndef MPZRZ_MUL_H_
#define MPZRZ_MUL_H_

void mpzrz_mul (mpz_t z[], mpz_t x[], mpz_t y[], mpz_t p);

#endif		//MPZRZ_MUL_H_

#ifndef MPZRZ_SQR_H_
#define MPZRZ_SQR_H_

void mpzrz_sqr (mpz_t z[], mpz_t x[], mpz_t p);

#endif		//MPZRZ_SQR_H_

#ifndef MPZRZ_EXP_H_
#define MPZRZ_EXP_H_

void mpzrz_exp (mpz_t result, mpz_t check, mpz_t base[3], mpz_t exp, mpz_t mod);

#endif		//MPZRZ_EXP_H_

#ifndef TRIAL_DIV_H_
#define TRIAL_DIV_H_

int trial_div (mpz_t N, int ndiv);

#endif		//TRIAL_DIV_H_

#ifndef TRIAL_DIV_SMALL_H_
#define TRIAL_DIV_SMALL_H_

int trail_div_small(int A, int r, int n, int y);

#endif		//TRIAL_DIV_SMALL_H_

#ifndef GCD_H_
#define GCD_H_

int gcd(int x, int y);

#endif		//GCD_H_

#ifndef CMP_H_
#define CMP_H_

int cmp(const void * a, const void * b);

#endif		//CMP_H_

#ifndef REMOVE_DUP_ARRAY_H_
#define REMOVE_DUP_ARRAY_H_

int remove_dup_array(int * arr, int init_size);

#endif		//REMOVE_DUP_ARRAY_H_

#ifndef GET_NOT_TO_RUN_2_8_10_H_
#define GET_NOT_TO_RUN_2_8_10_H_

int get_not_to_run_2_8_10 (int * not_to_run, int A, int r, int n, int y);

#endif		//GET_NOT_TO_RUN_2_8_10_H_

#ifndef GET_N_TO_RUN_2_8_10_H_
#define GET_N_TO_RUN_2_8_10_H_

int get_n_to_run_2_8_10 (int * to_run, int A, int r, int n, int nf, int y);

#endif		//GET_N_TO_RUN_2_8_10_H_

