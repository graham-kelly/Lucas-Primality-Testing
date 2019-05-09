#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "RWG2.h"
#include "RWG7.h"
#include "ptestio.h"

#define A atoi(argv[2])
#define r atoi(argv[3])
#define nf atoi(argv[5])
#define y_eta atoi(argv[6])

//	Run primality tests in RWG paper based on Lucas sequences
int main(int argc, char *argv[]) {
	if (argc == 7) {
		long int n = atoi(argv[4]);				// first value of n to be checked
		char fileStr[100];
		get_file_str(fileStr, argv);
		if (atoi(argv[1]) == 2) {											// do the tests from section 2
			if (r == 2) {
				while (n <= nf) {
					if (primality_test_2_10 (A, n, y_eta) == 1) {
						printf("\n\nPrime for exponent = %d\n\n", n);
						add_prime_file(n, fileStr);			// add to file
					}
					else {
						printf("%d\t", n);
					}
					n++;
				}
			}
			else {
				while (n <= nf) {
					if (primality_test_2_8 (A, r, n, y_eta) == 1) {
						printf("\n\nPrime for exponent = %d\n\n", n);
						add_prime_file(n, fileStr);			// add to file
					}
					else {
						printf("%d\t", n);
					}
					n++;
				}
			}
		}
		else {																// use the tests from section 7
			if (r == 2 || (r%4 != 1)) {
				printf ("Section 7 primality tests do not apply to r = 2 or r = 1 (mod 4)\n");
			}
			else {
				if (r == 5) {
					while (n <= nf) {
						if (primality_test_7_5 (A, n, y_eta) == 1) {			// specific to r = 5
							printf("\n\nPrime for exponent = %d\n\n", n);
							add_prime_file(n, fileStr);			// add to file						}
						}
						else {
							printf("%d\t", n);
						}
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
							printf("%d\t", n);
						}
						n++;
					}
				}
			}
		}
		printf("\n");
		remove_dup_file(fileStr);
	}
	else {
		printf("\n\tUse command line arguments: 2/7, A, r, n, nf, gamma/eta\n\n\tThese correspond to section 2/7\n\tA a small positive integer with gcd(A,r) = 1\n\tr a small prime\n\tn > nf positive integers\n\tgamma/eta = +/-1\n\n\n");
	}
	return 0;
}

