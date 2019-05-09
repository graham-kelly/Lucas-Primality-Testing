#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//	Adds prime's exponent to correct .txt file
void add_prime_file(int exp, char * fileStr) {
	FILE * f = fopen (fileStr, "a");
	if (f == NULL) {
		printf ("Error oepning file:\n\t%s\nn = %d is prime", fileStr, exp);
		return;
	}
	fprintf (f, "%d\n", exp);
	fclose(f);
	return;
}

//	Used by qsort in removeDupFile to order files (equivalent to <)
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
	while (!feof(f)) {
		nPrimes++;
		fscanf (f, "%d", primes + nPrimes);
		if (nPrimes == size) {
			size *= 2;
			primes = (int *) realloc (primes, size * sizeof(int));					// should implement error checking as well
		}
	}
	
	qsort(primes, nPrimes, sizeof(int), cmp);
	
	int newNPrimes = 0;
	int current = -1;
	for (int i = 0; i < nPrimes; i++) {
		if (current != *(primes + i)) {
			current = *(primes + i);
			*(primes + newNPrimes) = current;
			newNPrimes++;
		}
	}
	
	f = freopen (fileStr, "w", f);
	if (f == NULL) {printf ("Error reopening file, no duplicates removed.\n"); return;}
	for (int i = 0; i < newNPrimes; i++) {
		fprintf (f, "%d\n", *(primes + i));
	}
	fclose(f);
	return;
}

//	Converts the parameters passed from command line to a correct fileStr to add newly found primes to
void get_file_str(char * fileStr, char *argv[]) {
	strcat(fileStr,"../res/");
	strcat(fileStr,argv[2]);					// = A value used
	strcat(fileStr," x ");
	strcat(fileStr,argv[3]);					// = r value used
	if (atoi(argv[6]) == 1) {					// if gamma / eta is 1
		strcat(fileStr,"^n + ");
	}
	else {
		strcat(fileStr,"^n - ");
	}
	if (atoi(argv[1]) == 2) {					// if doing section 2 tests
		strcat(fileStr,"1.txt");
	}
	else {
		strcat(fileStr,"y_n.txt");
	}
	return;
}

