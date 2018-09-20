These files implement the primality tests from (RWG) Roetteger, Williams and Guy's paper 'Some Primality Tests that Eluded Lucas'

To compile the program, add the src, include and a new empty directory 'obj' to a primality testing file and run 'make' from a terminal window (while in src directory)

To test primality of numbers of the form
N = A * r^i + gamma						for section 2 where n <= i <= nf
N = A * r^i + eta * gamma_n(r)			for section 7 where n <= i <= nf, gamma_n(r) = sqrt(-1) (mod r^n)

Run the program with the following command line arguments (you will be prompted to do so if you forget)
2/7			section number of RWG which primality test is done from
A			a small positive integer with gcd(A,r) = 1
r			r a small prime (r != 2 and r != 3 (mod 4) for section 7 test)				* note that this is required for gamma_n(r) and N to exist
n			positive integer
nf			positive integer (nf >= n)
gamma/eta	+/-1


n such that N is prime will be printed to the file primes.txt
possible prime numbers will have n printed to the terminal window plus some additional constraints found
