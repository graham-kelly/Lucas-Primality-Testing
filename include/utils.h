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

#ifndef GCD_H_
#define GCD_H_

int gcd(int x, int y);

#endif		//GCD_H_

