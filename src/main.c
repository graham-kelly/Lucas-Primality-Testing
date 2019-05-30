#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "ptestio.h"
#include "utils.h"



//	Run primality tests in RWG paper based on Lucas sequences
int main(int argc, char *argv[]) {
	char fileStr[100];
	if (argc == 7) {
		get_file_name_OEIS (fileStr, argv);
		OEIS_testing (fileStr, atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		remove_dup_file(fileStr);
		return 0;
	}
	if (argc == 3) {
		get_file_name_SIDH (fileStr, argv);
		SIDH_testing(fileStr, atoi(argv[1]), atoi(argv[2]));
		return 0;
	}
	print_req_args();
	return 1;
}



















