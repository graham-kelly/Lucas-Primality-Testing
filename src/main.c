#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "RWG2.h"
#include "RWG7.h"
#include "ptestio.h"
#include "utils.h"

void OEIS_testing (char * fileStr, int test_num, int A, int r, int n, int nf, int y_eta);
void SIDH_testing (char * fileStr, int test_num, int A, int r, int n, int nf, int y_eta);

//	Run primality tests in RWG paper based on Lucas sequences
int main(int argc, char *argv[]) {
	if (argc == 7) {
		char fileStr[100];
		get_file_str(fileStr, argv);
		OEIS_testing (fileStr, atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		remove_dup_file(fileStr);
		return 0;
	}
	printf("\n\tUse command line arguments: 2/7, A, r, n, nf, gamma/eta\n\n\tThese correspond to section 2/7\n\tA a small positive integer with gcd(A,r) = 1\n\tr a small prime\n\tn > nf positive integers\n\tgamma/eta = +/-1\n\n\n");
	return 1;
}

int LARGEST_IN_SMALL_PRIMES = 15485863;
/*
Parameters:
	fileStr			file to output results to
	test_num		2 => do tests from section 2 of RWG, 7 => tests from section 7
	A, r, y_eta		integers as described in RWG Theorems 2.8, 2.10, 7.2, 7.4, 7.5
	n, nf			test for primes using above constants and exponent n, n+1, n+2, ... , nf-2, nf-1, nf
*/
void OEIS_testing (char * fileStr, int test_num, int A, int r, int n, int nf, int y_eta) {
	if (test_num == 2) {															// do the tests from section 2
		mpz_t N; mpz_init (N);
		while (n < (log(LARGEST_IN_SMALL_PRIMES - y_eta)-log(A))/log(r)) {				// if the potential prime is <= this big it will be in smallPrimes.txt
			mpz_ui_pow_ui (N, r, n);
			mpz_mul_ui (N, N, A);
			if (y_eta == 1) {
				mpz_add_ui (N, N, 1);
			}
			else {
				mpz_sub_ui (N, N, 1);
			}						// N = A*r^n+y
			int divRes = trial_div(N, 0);
			if (divRes) {			// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
				if (mpz_cmpabs_ui (N, divRes) == 0) {
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file(n, fileStr);			// add to file
				}
			}
			printf("%d,", n);
			fflush(stdout);
			n++;
		}
		mpz_clear(N);
		int *to_run = (int *) calloc((nf - n + 1), sizeof(int));
		int size_to_run = get_n_to_run_2_8_10 (to_run, A, r, n, nf, y_eta);
		if (r == 2) {
			for (int i = 0; i < size_to_run; i++) {
				if (primality_test_2_10 (A, to_run[i], y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", to_run[i]);
					add_prime_file(to_run[i], fileStr);			// add to file
				}
				printf("%d,", to_run[i]);
				fflush(stdout);
			}
		}
		else {
			for (int i = 0; i < size_to_run; i++) {
				if (primality_test_2_8 (A, r, to_run[i], y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", to_run[i]);
					add_prime_file(to_run[i], fileStr);			// add to file
				}
				printf("%d,", to_run[i]);
				fflush(stdout);
			}
		}
	}
	else {																// use the tests from section 7
		if (r == 2 || (r%4 != 1)) {
			printf ("Section 7 primality tests do not apply to r = 2 or r = 1 (mod 4)\n");
			return;
		}
		if (r == 5) {
			while (n <= nf) {
				if (primality_test_7_5 (A, n, y_eta) == 1) {			// specific to r = 5
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file(n, fileStr);			// add to file						}
				}
				else {
					printf("%d,", n);
				}
				fflush(stdout);
				n++;
			}
		}
		else {			// use primality testing from both Theorem 7.2 and 7.4 as they require computing the same sequences are have the same requirements
			while (n <= nf) {
				if (primality_test_7_2_4 (A, r, n, y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file(n, fileStr);			// add to file
				}
				else {
					printf("%d,", n);
				}
				fflush(stdout);
				n++;
			}
		}
	}
	printf("\n");
	return;
}

void SIDH_testing (char * fileStr, int test_num, int A, int r, int n, int nf, int y_eta) {
	return;
}




















