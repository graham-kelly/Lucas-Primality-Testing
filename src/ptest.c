#include <gmp.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "primalityRWG.h"
#include "ptest.h"

int main(int argc, char *argv[]) {
	char fileStr[100];
	if (argc == 7) {
		get_file_name_OEIS (fileStr, argv);
		mpz_t A; mpz_init (A);
		mpz_set_ui (A, atoi(argv[2]));
		OEIS_testing (fileStr, atoi(argv[1]), A, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));		// fileStr, test_num, A, r, n, nf, y
		mpz_clear (A);
		return 0;
	}
	if (argc == 3) {
		get_file_name_SIDH (fileStr, argv);
		printf("%s\n",fileStr + 7);
		SIDH_testing(fileStr, atoi(argv[1]), atoi(argv[2]));
		return 0;
	}
	if (argc == 9) {
		mpz_t A; mpz_init (A);
		mpz_set_ui (A, atoi(argv[2]));
		mpz_ui_pow_ui (A, atoi(argv[2]), atoi(argv[3]));
		if (atoi(argv[4]) >= 0) {
			mpz_add_ui (A, A, atoi(argv[4]));
		}
		else {
			mpz_sub_ui (A, A, abs(atoi(argv[4])));
		}
		get_file_name2_OEIS (fileStr, atoi(argv[1]), A, atoi(argv[5]), atoi(argv[8]));
		OEIS_testing (fileStr, atoi(argv[1]), A, atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));		// fileStr, test_num, A, r, n, nf, y
		mpz_clear (A);
		return 0;
	}
	print_req_args();
	return 1;
}

int LARGEST_IN_SMALL_PRIMES = 15485863;
/*
Parameters:
	fileStr			file to output results to
	test_num		2 => do tests from section 2 of RWG, 7 => tests from section 7
	A, r, y_eta		integers as described in RWG Theorems 2.8, 2.10, 7.2, 7.4, 7.5
	n, nf			test for primes using above constants and exponent n, n+1, n+2, ... , nf-2, nf-1, nf
Returns:
	prints results of primality test to file specified in fileStr
Runtime:	O(log(log(N)) * log(N)^2)						for tests from section 7
			O(min(r * log(N)^2, log(log(N)) * log(N)^2))	for tests from section 2
*/
void OEIS_testing (char * fileStr, int test_num, mpz_t A, int r, int n, int nf, int y_eta) {
	time_t ti;
	if (test_num == 2) {															// do the tests from section 2
		int *to_run = (int *) calloc((nf - n + 1), sizeof(int));
		int size_to_run = get_n_to_run_2_8_10 (to_run, A, r, n, nf, y_eta);
		for (int i = 0; i < size_to_run; i++) {
			ti = time(NULL);
			if (primality_test_2 (A, r, to_run[i], y_eta) == 1) {
				add_prime_file_OEIS(to_run[i], fileStr, 1, time(NULL) - ti);			// add to file
			}
			else {
				add_prime_file_OEIS(to_run[i], fileStr, 0, time(NULL) - ti);			// add to file
			}
		}
	}
	else {																// use the tests from section 7
		if (r%4 != 1) {
			printf ("Section 7 primality tests only apply to r == 1 (mod 4)\n");
			return;
		}
		for (; n <= nf; n++) {
			ti = time(NULL);
			if (primality_test_7 (A, r, n, y_eta) == 1) {
				add_prime_file_OEIS(n, fileStr, 1, time(NULL) - ti);			// add to file
			}
			else {
				add_prime_file_OEIS(n, fileStr, 0, time(NULL) - ti);			// add to file
			}
		}
	}
	return;
}

/*
Parameters:
	fileStr			file to output results to
	bit_length		2 => do tests from section 2 of RWG, 7 => tests from section 7
	r				integer as described in RWG Theorems 2.8, 2.10
Returns:
	prints results of primality test to file specified in fileStr
Runtime:	O(bit_length * min(r * log(N)^2, log(log(N)) * log(N)^2))
*/
void SIDH_testing (char * fileStr, int bit_length, int r) {
	mpz_t A; mpz_init (A);
	double log2r = log2(r);
	int x = 0;
	int y;
	if (bit_length - 40 > 0) {
		x = (bit_length / 2) - 20;
	}
	while (x < bit_length / 2 + 20) {
		y = 0;
		if (bit_length - x - 20 > 0) {
			y = (int) ceil((bit_length - x - 20) / log2r);
		}
		int loopControl = (int) floor ((bit_length - x) / log2r);
		while (x >= y * log2r && y < loopControl) {
			int innerLoopControl = bit_length - x - ((int) ceil(y * log2r));
			for (int f = 1; f < innerLoopControl; f+=2) {
				if (f % r == 0) {
					f+=2;
				}
				mpz_ui_pow_ui (A, r, y);
				mpz_mul_ui (A, A, f);
				if (primality_test_2 (A, 2, x, -1) == 1) {
					add_prime_file_SIDH(f, x, y, fileStr);				// add to file
				}
				printf("%d/%d/%d,", f, x, y);
			}
			y++;
		}
		while (y < loopControl) {
			int innerLoopControl = bit_length - x - ((int) ceil(y * log2r));
			for (int f = 1; f < innerLoopControl; f+=2) {
				if (f % r == 0) {
					f+=2;
				}
				mpz_set_ui (A, f);
				mpz_mul_2exp (A, A, x);
				if (primality_test_2 (A, r, y, -1) == 1) {
					add_prime_file_SIDH(f, x, y, fileStr);				// add to file
				}
				printf("%d/%d/%d,", f, x, y);
			}
			y++;
		}
		x++;
	}
	mpz_clear (A);
	printf("\n");
	return;
}

/*	Get information about which groups of n will always give composite numbers for N=Ar^n+y
Parameters:
	A,r,n,y		are positive integers as given in RWG section 2 primality tests (theorems 2.8 / 2.10)
Returns:
	not_to_run			contains triples of p, x, z such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
	size_not_to_run		
Runtime:	O(p_n * M(A * r^(p_n))) 
	where p_n is the nth prime number

Note that initial size of not_to_run should be 2*n
*/
int get_not_to_run_2_8_10 (int * not_to_run, mpz_t A, int r, int n, int y) {
	extern _Bool *p_loaded_primes;
	extern int *p_small_primes;
	if (!*p_loaded_primes) {
		load_small_primes();
	}
	int size_not_to_run = 0;
	mpz_t test_N; mpz_init (test_N);
	mpz_t rEXPj; mpz_init (rEXPj);
	_Bool found_x;
	_Bool found_z;
	for (int i = 0; i < n; i++) {
		if (*p_small_primes + i >= n) {
			break;
		}
		mpz_tdiv_r_ui (test_N, A, *p_small_primes + i);
		if (y == 1) {
			mpz_add_ui (test_N, test_N, 1);
		}
		else {
			mpz_sub_ui (test_N, test_N, 1);
		}						// test_N = A+y
		mpz_set_ui (rEXPj, 1);
		found_x = 0;
		found_z = 0;
		for (int j = 1; j < *p_small_primes + i; j++) {					// find x, z such that A*r^x+z = 0 (mod p) and r^z = 1 (mod p)
			if (y == 1) {
				mpz_sub_ui (test_N, test_N, 1);
				mpz_mul_ui (test_N, test_N, r);
				mpz_add_ui (test_N, test_N, 1);
			}
			else {
				mpz_add_ui (test_N, test_N, 1);
				mpz_mul_ui (test_N, test_N, r);
				mpz_sub_ui (test_N, test_N, 1);
			}						// test_N = A*(r^j) + y
			if (!found_x && mpz_divisible_ui_p (test_N, *p_small_primes + i)) {
				*(not_to_run + 2 * size_not_to_run) = j;
				found_x = 1;
			}
			mpz_mul_ui (rEXPj, rEXPj, r);
			if (!found_z && mpz_congruent_ui_p(rEXPj, 1, *p_small_primes + i)) {
				*(not_to_run + 2 * size_not_to_run + 1) = j;
				found_z = 1;
			}
			if (found_x && found_z) {					// not_to_run[size] = (x, z) such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
				size_not_to_run++;
				break;
			}
		}
	}
	mpz_clear (test_N);
	mpz_clear (rEXPj);
	not_to_run = (int *) realloc (not_to_run, sizeof(int) * 2 * size_not_to_run);
	return size_not_to_run;
}

double SIEVING_POW = 0.8

int get_n_to_run_2_8_10 (int * to_run, mpz_t A, int r, int n, int nf, int y) {
	int *not_to_run = (int *) malloc(sizeof(int) * 2 * nf);						// get p, x, z such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
	int size_not_to_run = get_not_to_run_2_8_10 (not_to_run, A, r, (int) pow(nf, SIEVING_POW), y);	// no n = x+kz for any integer k should be tested (all divisble by p)
	int k = 0;
	int exp;
	for (int i = 0; i < size_not_to_run; i++) {		// for each small prime tested
		exp = (int) ceil ((n - *(not_to_run + 2 * i)) / (double) *(not_to_run + 2 * i + 1)) * (*(not_to_run + 2 * i + 1)) + *(not_to_run + 2 * i) - n;		// get initial exponent's index (x - n + z*ceil( (n-x)/z) )
		while (exp < nf - n + 1) {
			*(to_run + exp) = 1;					// this N with exponent (exp+n) will be divisible by some small prime
			exp += *(not_to_run + 2 * i + 1);		// add z (next iteration remove new exponent from the list if still in range)
		}
	}
	int size_to_run = 0;
	for (int i = 0; i < nf - n + 1; i++) {
		if (!*(to_run + i)) {						// has not been marked as divisible by some x + kz above
			*(to_run + size_to_run) = i + n;		// this integer should be run
			size_to_run++;
		}
	}
	to_run = (int *) realloc (to_run, sizeof(int) * size_to_run);
	return size_to_run;
}

//	Adds prime's exponent to correct .txt file
void add_prime_file_OEIS (int exp, char * fileStr, _Bool prime, time_t time) {
	FILE * f = fopen (fileStr, "a");
	if (f == NULL) {
		printf ("Error opening file:\n\t%s\nn = %d, prime = %d, time = %d\n", fileStr, exp, prime, time);
		return;
	}
	if (prime) {
		fprintf (f, "p%d:%d\n", exp, time);
	}
	else {
		fprintf (f, "%d:%d\n", exp, time);
	}
	fclose(f);
	return;
}

//	Adds prime's exponent to correct .txt file
void add_prime_file_SIDH (int f, int x, int y, char * fileStr) {
	FILE * file = fopen (fileStr, "a");
	if (file == NULL) {
		printf ("Error opening file:\n\t%s\nf = %d, x = %d, y = %d is prime\n", fileStr, f, x, y);
		return;
	}
	fprintf (file, "%d, %d, %d\n", f, x, y);
	fclose(file);
	return;
}

//	Converts the parameters passed from command line to a correct fileStr to add newly found primes to
void get_file_name_OEIS(char * fileStr, char *args[]) {
	strcpy(fileStr,"../res/");
	strcat(fileStr,args[2]);					// = A value used
	strcat(fileStr,"x");
	strcat(fileStr,args[3]);					// = r value used
	if (atoi(args[6]) == 1) {					// if gamma / eta is 1
		strcat(fileStr,"^n+");
	}
	else {
		strcat(fileStr,"^n-");
	}
	if (atoi(args[1]) == 2) {					// if doing section 2 tests
		strcat(fileStr,"1.txt");
	}
	else {
		strcat(fileStr,"y_n.txt");
	}
	return;
}

//	Converts the parameters to a correct fileStr to add newly found primes to (with mpz_t A)
void get_file_name2_OEIS(char * fileStr, int test, mpz_t A, int r, int y_eta) {
	strcpy(fileStr,"../res/");
	mpz_get_str (fileStr + 7, 10, A);			// A
	strcat(fileStr,"x");
	sprintf(fileStr + strlen(fileStr), "%d", r);				// r value used
	if (y_eta == 1) {					// if gamma / eta is 1
		strcat(fileStr,"^n+");
	}
	else {
		strcat(fileStr,"^n-");
	}
	if (test == 2) {					// if doing section 2 tests
		strcat(fileStr,"1.txt");
	}
	else {
		strcat(fileStr,"y_n.txt");
	}
	return;
}

void get_file_name_SIDH(char * fileStr, char * args[]) {
	strcpy(fileStr,"../res/SIDH-");
	strcat(fileStr,args[1]);					// = bitlength
	strcat(fileStr,"-bit:f(2^n)(");
	strcat(fileStr,args[2]);					// = r used (small prime)
	strcat(fileStr,"^m)-1");
	return;
}

void print_req_args() {
	printf("\nFor OEIS testing use one of the command line arguments:\n");
	printf("\t2/7, A, r, n, nf, gamma/eta\n");
	printf("\t2/7, A_base, A_exp, A_offset, r, n, nf, gamma/eta\n");
	printf("\nSection 2 or 7 of RWG\n");
	printf("A a small positive integer with gcd(A,r) = 1\n");
	printf("r a small prime\n");
	printf("n < nf positive integers\n");
	printf("gamma/eta = +/-1\n\n\n");
	printf("For SIDH testing use command line arguments:\n");
	printf("\tBitLength, r\n");
	printf("Bitlength = total bitlength for prime numbers to be found\n");
	printf("r a small prime (!=2)\n\n\n");
	return;
}

