#ifndef GET_CUV_I_H_
#define GET_CUV_I_H_

_Bool get_cUV_i (mpz_t U_n, mpz_t V_n, int QPP[3], mpz_t i, mpz_t N);

#endif		//GET_CUV_I_H_

#ifndef GET_KL_M_H_
#define GET_KL_M_H_

void get_KL_m (mpz_t K_m, mpz_t L_m, mpz_t m, int QPP[3], mpz_t N);

#endif		//GET_KL_M_H_

#ifndef GET_RST_0_H_
#define GET_RST_0_H_

_Bool get_RST_0 (mpz_t RST[3], mpz_t rEXPn, int QPP[3], mpz_t N);

#endif		//GET_RST_0_H_

#ifndef GET_NEXT_RST_I_H_
#define GET_NEXT_RST_I_H_

void get_next_RST_i (mpz_t newRST[3], mpz_t oldRST[3], mpz_t tmp_val[6], int QPP[3], int r, mpz_t XY_array[], mpz_t N);

#endif		//GET_NEXT_RST_I_H_

#ifndef GET_RST_I_H_
#define GET_RST_I_H_

_Bool get_RST_i (mpz_t rop[3], int i, int QPP[3], int r, mpz_t rEXPn, mpz_t XY_array[], mpz_t N);

#endif		//GET_RST_I_H_

#ifndef GET_ST_0_H_
#define GET_ST_0_H_

_Bool get_ST_0(mpz_t S_i, mpz_t T_i, int QPP[3], mpz_t rEXPn, mpz_t tmp_val[2], mpz_t N);

#endif		//GET_ST_0_H_

#ifndef GET_NEXT_ST_I_H_
#define GET_NEXT_ST_I_H_


void get_next_ST_i (mpz_t S_i, mpz_t T_i, int Delta, mpz_t tmp_val[3], mpz_t N);

#endif		//GET_NEXT_ST_I_H_






