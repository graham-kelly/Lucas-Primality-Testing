#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "RWG2.h"
#include "RWG7.h"

/*
This program is a primality test for numbers of the form Ar^n+gamma or Ar^n+eta(gamma_n(r))
It implements the Primality tests from sections 2 and 7 of Roetteger, Williams, and Guy's paper "Some primality tests that eluded Lucas"
Specifically from Theorems 2.8, 2.10, 7.2, 7,4 and 7.5
Runtime for the algorithm is:

_____________ (based on input)

Find integers n such that Ar^n+gamma or Ar^n+eta(gamma_n(r)) is prime in primes.txt after running the program
*/
int main(int argc, char *argv[]) {
	if (argc == 7) {
		long int A = atoi(argv[2]);
		long int r = atoi(argv[3]);				// a small prime (r = 2 is not available for the test from section 7)
		long int n = atoi(argv[4]);				// first value of n to be checked
		long int nf = atoi(argv[5]);			// last value of n to be checked (all intermediate values will also be checked)
		short y_eta = atoi(argv[6]);			// = +/- 1
		FILE *f = fopen ("primes.txt", "w");
		if (f == NULL) {printf ("Error oepning file.\n"); return 1;}
		if (atoi(argv[1]) == 2) {											// do the tests from section 2
			fprintf (f, "Some primes of the type: %d * %d ^ n + (%d) have n =\n", A, r, y_eta);
			if (r == 2) {
				while (n <= nf) {
					if (primality_test_2_10 (A, n, y_eta) == 1) {
						fprintf (f, "%d\n", n);
					}
					n++;
				}
			}
			else {
				while (n <= nf) {
					if (primality_test_2_8 (A, r, n, y_eta) == 1) {
						fprintf (f, "%d\n", n);
					}
					n++;
				}
			}
		}
		else {																// use the tests from section 7
			fprintf (f, "Some primes of the type: %d * %d ^ n + (%d) * gamma_r(n) have n =\n", A, r, y_eta);
			if (r == 2 || (r%4 != 1)) {
				printf ("Section 7 primality tests do not apply to r = 2 or r != 1 (mod 4)\n");
			}
			else {
				if (r == 5 && 0) {										// ******************	REMOVE 0 WHEN 7.5 IS COMPLETE	******************
					while (n <= nf) {
						if (primality_test_7_5 (A, n, y_eta)) {			// specific to r = 5
							fprintf (f, "%d\n", n);
						}
						n++;
					}
				}
				else {			// use primality testing from both Theorem 7.2 and 7.4 as they require computing the same sequences are have the same requirements
					while (n <= nf) {
						if (primality_test_7_2_4 (A, r, n, y_eta)) {
							fprintf (f, "%d\n", n);
						}
						n++;
					}
				}
			}
		}
		fclose(f);
	}
	else {
		printf("\n\tUse command line arguments: 2/7, A, r, n, nf, gamma/eta\n\n\tThese correspond to section 2/7\n\tA a small positive integer with gcd(A,r) = 1\n\tr a small prime\n\tn > nf positive integers\n\tgamma/eta = +/-1\n\n\n");
	}
	return 0;
}

