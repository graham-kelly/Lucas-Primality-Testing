#ifndef PRIMALITYRWG_H_
#define PRIMALITYRWG_H_

#ifndef PRIMALITY_TEST_2_8_H_
#define PRIMALITY_TEST_2_8_H_

int primality_test_2_8 (mpz_t A, int r, int n, short y);

#endif		//PRIMALITY_TEST_2_8_H_

#ifndef PRIMALITY_TEST_2_10_H_
#define PRIMALITY_TEST_2_10_H_

int primality_test_2_10 (mpz_t A, int n, short y);

#endif		//PRIMALITY_TEST_2_10_H_

#ifndef GET_U_I_H_
#define GET_U_I_H_

_Bool get_u_i (mpz_t result, mpz_t i, int p, int d, mpz_t N);

#endif		//GET_U_I_H_

#ifndef GET_UV_I_H_
#define GET_UV_I_H_

_Bool get_uv_i (mpz_t u_i, mpz_t v_i, mpz_t i, int p, int d, mpz_t N);

#endif		//GET_UV_I_H_

#ifndef GET_V_I_H_
#define GET_V_I_H_

_Bool get_v_i (mpz_t result, mpz_t i, int p, int d, mpz_t N);

#endif		//GET_V_I_H_

#ifndef PRIMALITY_TEST_7_2_4_H_
#define PRIMALITY_TEST_7_2_4_H_

int primality_test_7_2_4 (int A, int r, int n, int eta);

#endif		//PRIMALITY_TEST_7_2_4_H_

#ifndef PRIMALITY_TEST_7_5_H_
#define PRIMALITY_TEST_7_5_H_

int primality_test_7_5 (int A, int n, int eta);

#endif		//PRIMALITY_TEST_7_5_H_

#ifndef GET_CUV_I_H_
#define GET_CUV_I_H_

_Bool get_cUV_i (mpz_t U_n, mpz_t V_n, int QPP[3], mpz_t i, mpz_t N);

#endif		//GET_CUV_I_H_

#ifndef TRIAL_DIV_H_
#define TRIAL_DIV_H_

int trial_div (mpz_t N, int ndiv);

#endif		//TRIAL_DIV_H_


#endif		//PRIMALITYRWG_H_

