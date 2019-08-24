#ifndef GET_H_K_H_
#define GET_H_K_H_

void get_HI_k (mpz_t rop1, mpz_t rop2, int k, mpz_t XY_array[], mpz_t tmp_val[6], mpz_t N);

#endif		//GET_H_K_H_

#ifndef PRECOMPUTE_XY_EXP_H_
#define PRECOMPUTE_XY_EXP_H_

void get_XY_exp (mpz_t XY_array[], mpz_t X, mpz_t Y, int k, mpz_t N);

#endif		//PRECOMPUTE_XY_EXP_H_

/*
#ifndef PRIMALITY_TEST_7_1_H_
#define PRIMALITY_TEST_7_1_H_

int primality_test_7_1 (int A, int r, int n, int eta);

#endif		//PRIMALITY_TEST_7_1_H_
*/

#ifndef PRIMALITY_TEST_7_H_
#define PRIMALITY_TEST_7_H_

int primality_test_7 (mpz_t A, int r, int n, int eta);

#endif		//PRIMALITY_TEST_7_H_

#ifndef PRIMALITY_TEST_7_2_4_H_
#define PRIMALITY_TEST_7_2_4_H_

int primality_test_7_2_4 (mpz_t A, int r, int n, int eta);

#endif		//PRIMALITY_TEST_7_2_4_H_

#ifndef PRIMALITY_TEST_7_5_H_
#define PRIMALITY_TEST_7_5_H_

int primality_test_7_5 (mpz_t A, int n, int eta);

#endif		//PRIMALITY_TEST_7_5_H_

#ifndef FIND_QPP_7_2_4_H_
#define FIND_QPP_7_2_4_H_

_Bool find_QPP_7_2_4 (int QPP[3], mpz_t N);

#endif		//FIND_QPP_7_2_4_H_

#ifndef TRY_BEST_QPP_7_2_4_H_
#define TRY_BEST_QPP_7_2_4_H_

_Bool try_best_QPP_7_2_4(int QPP[3], mpz_t N);

#endif		//TRY_BEST_QPP_7_2_4_H_

#ifndef FIND_QPP_7_5_H_
#define FIND_QPP_7_5_H_

_Bool find_QPP_7_5 (int QPP[3], mpz_t constant, mpz_t N);

#endif		//FIND_QPP_7_5_H_

