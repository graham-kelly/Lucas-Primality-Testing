#include <gmp.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "primalityRWG.h"
#include "ptest.h"

int main(int argc, char *argv[]) {
	char fileStr[100];
	if (argc == 7) {
		get_file_name_OEIS (fileStr, argv);
		printf("%s\n",fileStr + 7);
		OEIS_testing (fileStr, atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		remove_dup_file(fileStr);
		return 0;
	}
	if (argc == 3) {
		get_file_name_SIDH (fileStr, argv);
		printf("%s\n",fileStr + 7);
		SIDH_testing(fileStr, atoi(argv[1]), atoi(argv[2]));
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
			int divRes = trial_div(N, 0);			// trial division by first million primes (can change 0 to any number between 1 and 1,000,000 to only try dividing by that many primes)
			if (divRes) {
				if (mpz_cmpabs_ui (N, divRes) == 0) {
					printf("\n\nPrime for exponent = %d\n\n", n);
					add_prime_file_OEIS(n, fileStr);			// add to file
				}
			}
			printf("%d,", n);
			fflush(stdout);
			n++;
		}
		mpz_clear(N);
		mpz_t Az; mpz_init (Az);
		mpz_set_ui (Az, A);
		int *to_run = (int *) calloc((nf - n + 1), sizeof(int));
		int size_to_run = get_n_to_run_2_8_10 (to_run, Az, r, n, nf, y_eta);
		if (r == 2) {
			for (int i = 0; i < size_to_run; i++) {
				if (primality_test_2_10 (Az, to_run[i], y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", to_run[i]);
					add_prime_file_OEIS(to_run[i], fileStr);			// add to file
				}
				printf("%d,", to_run[i]);
				fflush(stdout);
			}
		}
		else {
			for (int i = 0; i < size_to_run; i++) {
				if (primality_test_2_8 (Az, r, to_run[i], y_eta) == 1) {
					printf("\n\nPrime for exponent = %d\n\n", to_run[i]);
					add_prime_file_OEIS(to_run[i], fileStr);			// add to file
				}
				printf("%d,", to_run[i]);
				fflush(stdout);
			}
		}
		mpz_clear (Az);
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
					add_prime_file_OEIS(n, fileStr);			// add to file						}
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
					add_prime_file_OEIS(n, fileStr);			// add to file
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
				if (primality_test_2_10 (A, x, -1) == 1) {
					printf("\n\nPrime for values: f = %d, x = %d, y = %d\n\n", f, x, y);
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
				if (primality_test_2_8 (A, r, y, -1) == 1) {
					printf("\n\nPrime for values: f = %d, x = %d, y = %d\n\n", f, x, y);
					add_prime_file_SIDH(f, x, y, fileStr);				// add to file
				}
				printf("%d/%d/%d,", f, x, y);
			}
			y++;
		}
		x++;
		fflush(stdout);													// ensure to inform the user of progress
	}
	// select A = f * 2^x, r^n such that r^n has about x bits (>x => use 2.8 else use 2.10)
	// can use small f, want 
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
	FILE *f = fopen ("primessmall.txt", "r");
	if (f == NULL) {printf("Error opening file.\n"); return -1;}
	int size_not_to_run = 0;
	int p;
	mpz_t test_N; mpz_init (test_N);
	mpz_t rEXPj; mpz_init (rEXPj);
	_Bool found_x;
	_Bool found_z;
	for (int i = 0; i < n; i++) {
		fscanf (f, "%d", &p);														// get prime from primesSmall
		if (p >= n) {
			break;
		}
		mpz_tdiv_r_ui (test_N, A, p);
		if (y == 1) {
			mpz_add_ui (test_N, test_N, 1);
		}
		else {
			mpz_sub_ui (test_N, test_N, 1);
		}						// test_N = A+y
		mpz_set_ui (rEXPj, 1);
		found_x = 0;
		found_z = 0;
		for (int j = 1; j < p; j++) {					// find x, z such that A*r^x+z = 0 (mod p) and r^z = 1 (mod p)
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
			if (!found_x && mpz_divisible_ui_p (test_N, p)) {
				*(not_to_run + 2 * size_not_to_run) = j;
				found_x = 1;
			}
			mpz_mul_ui (rEXPj, rEXPj, r);
			if (!found_z && mpz_congruent_ui_p(rEXPj, 1, p)) {
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
	fclose(f);
	not_to_run = (int *) realloc (not_to_run, sizeof(int) * 2 * size_not_to_run);
	return size_not_to_run;
}

int get_n_to_run_2_8_10 (int * to_run, mpz_t A, int r, int n, int nf, int y) {
	int *not_to_run = (int *) malloc(sizeof(int) * 2 * nf);						// get p, x, z such that A*r^x+y = 0 (mod p) and r^z = 1 (mod p)
	int size_not_to_run = get_not_to_run_2_8_10 (not_to_run, A, r, nf, y);	// no n = x+kz for any integer k should be tested (all divisble by p)
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
void add_prime_file_OEIS(int exp, char * fileStr) {
	FILE * f = fopen (fileStr, "a");
	if (f == NULL) {
		printf ("Error oepning file:\n\t%s\nn = %d is prime\n", fileStr, exp);
		return;
	}
	fprintf (f, "%d\n", exp);
	fclose(f);
	return;
}

//	Adds prime's exponent to correct .txt file
void add_prime_file_SIDH(int f, int x, int y, char * fileStr) {
	FILE * file = fopen (fileStr, "a");
	if (file == NULL) {
		printf ("Error oepning file:\n\t%s\nf = %d, x = %d, y = %d is prime\n", fileStr, f, x, y);
		return;
	}
	fprintf (file, "%d, %d, %d\n", f, x, y);
	fclose(file);
	return;
}

//	Comparison of two integers for use with qsort (equivalent to <)
int cmp (const void * a, const void * b) {
	return (*(int *)a - *(int *)b);
}

//	Removes any duplicated prime's exponents from the specified .txt file
void remove_dup_file(char * fileStr) {
	int size = 100;
	int nPrimes = 0;
	int * primes = (int *) calloc(size, sizeof(int));
	FILE * f = fopen (fileStr, "r");
	if (f == NULL) {printf ("Error opening file, no duplicates removed.\n"); return;}
	fscanf (f, "%d", primes + nPrimes);
	while (!feof(f)) {																// get exponents for all the prime numbers found
		nPrimes++;																	// count number found
		fscanf (f, "%d", primes + nPrimes);
		if (nPrimes == size) {														// increase allocation for array of primes' exponents if necessary
			size *= 2;
			primes = (int *) realloc (primes, size * sizeof(int));					// should implement error checking as well
		}
	}
	qsort(primes, nPrimes, sizeof(int), cmp);										// sort primes found
	int current = *primes;
	int new_size = 1;
	for (int i = 1; i < nPrimes; i++) {												// for each number
		if (current != *(primes + i)) {												// if not yet seen
			current = *(primes + i);												// update current
			*(primes + new_size) = current;											// update new list of primes
			new_size++;																// update size of new list of primes
		}
	}
	primes = (int *) realloc (primes, sizeof(int) * new_size);						// reallocate only the necessary memory
	nPrimes = new_size;
	f = freopen (fileStr, "w", f);													// overwrite the old file with the new data
	if (f == NULL) {printf ("Error reopening file, no duplicates removed.\n"); return;}
	for (int i = 0; i < nPrimes; i++) {
		fprintf (f, "%d\n", *(primes + i));
	}
	fclose(f);
	return;
}

//	Converts the parameters passed from command line to a correct fileStr to add newly found primes to
void get_file_name_OEIS(char * fileStr, char *args[]) {
	strcpy(fileStr,"../res/");
	strcat(fileStr,args[2]);					// = A value used
	strcat(fileStr," x ");
	strcat(fileStr,args[3]);					// = r value used
	if (atoi(args[6]) == 1) {					// if gamma / eta is 1
		strcat(fileStr,"^n + ");
	}
	else {
		strcat(fileStr,"^n - ");
	}
	if (atoi(args[1]) == 2) {					// if doing section 2 tests
		strcat(fileStr,"1.txt");
	}
	else {
		strcat(fileStr,"y_n.txt");
	}
	return;
}

void get_file_name_SIDH(char * fileStr, char * args[]) {
	strcpy(fileStr,"../res/SIDH ");
	strcat(fileStr,args[1]);					// = bitlength
	strcat(fileStr,"-bit: f x 2^n x ");
	strcat(fileStr,args[2]);					// = r used (small prime)
	strcat(fileStr,"^m - 1");
	return;
}

void print_req_args() {
	printf("\nFor OEIS testing use command line arguments:\n");
	printf("\t2/7, A, r, n, nf, gamma/eta\n");
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

