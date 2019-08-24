#ifndef OEIS_TESTING_H_
#define OEIS_TESTING_H_

void OEIS_testing (char * fileStr, int test_num, mpz_t A, int r, int n, int nf, int y_eta);

#endif		//OEIS_TESTING_H_

#ifndef SIDH_TESTING_H_
#define SIDH_TESTING_H_

void SIDH_testing (char * fileStr, int bit_length, int r);

#endif		//SIDH_TESTING_H_

#ifndef GET_NOT_TO_RUN_2_8_10_H_
#define GET_NOT_TO_RUN_2_8_10_H_

int get_not_to_run_2_8_10 (int * not_to_run, mpz_t A, int r, int n, int y);

#endif		//GET_NOT_TO_RUN_2_8_10_H_

#ifndef GET_N_TO_RUN_2_8_10_H_
#define GET_N_TO_RUN_2_8_10_H_

int get_n_to_run_2_8_10 (int * to_run, mpz_t A, int r, int n, int nf, int y);

#endif		//GET_N_TO_RUN_2_8_10_H_

#ifndef ADD_PRIME_FILE_OEIS_H_
#define ADD_PRIME_FILE_OEIS_H_

void add_prime_file_OEIS(int exp, char * fileStr, _Bool prime, time_t time);

#endif		//ADD_PRIME_FILE_OEIS_H_

#ifndef ADD_PRIME_FILE_SIDH_H_
#define ADD_PRIME_FILE_SIDH_H_

void add_prime_file_SIDH(int f, int x, int y, char * fileStr);

#endif		//ADD_PRIME_FILE_SIDH_H_

#ifndef GET_FILE_NAME_OEIS_H_
#define GET_FILE_NAME_OEIS_H_

void get_file_name_OEIS(char * fileStr, char *argv[]);

#endif		//GET_FILE_NAME_OEIS_H_

#ifndef GET_FILE_NAME2_OEIS_H_
#define GET_FILE_NAME2_OEIS_H_

void get_file_name2_OEIS(char * fileStr, int test, mpz_t A, int r, int y_eta);

#endif		//GET_FILE_NAME2_OEIS_H_

#ifndef GET_FILE_NAME_SIDH_H_
#define GET_FILE_NAME_SIDH_H_

void get_file_name_SIDH(char * fileStr, char * argv[]);

#endif		//GET_FILE_NAME_SIDH_H_

#ifndef PRINT_REQ_ARGS_H_
#define PRINT_REQ_ARGS_H_

void print_req_args();

#endif		//PRINT_REQ_ARGS_H_

