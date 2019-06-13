These files implement the primality tests from (RWG) Roetteger, Williams and Guy's paper 'Some Primality Tests that Eluded Lucas'

To compile the program, add the src, inc, obj, and res folders (last two may be empty) to a primality testing file and run 'make' from a terminal window (while in src directory)
This uses the gcc compiler and makes use of the gnu mp library (https://gmplib.org/manual/) and has been run on the University of Calgary's Linux servers

To test primality of numbers of the form:
N = A * r^i + gamma						for section 2 where n <= i <= nf
N = A * r^i + eta * gamma_n(r)			for section 7 where n <= i <= nf, gamma_n(r) = sqrt(-1) (mod r^n)

Run the program with the following command line arguments (you will be prompted to do so if you forget)
2/7			section number of RWG which primality test is done from
A			a small positive integer with gcd(A,r) = 1
r			r a small prime (r != 2 and r != 3 (mod 4) for section 7 test)				* note that this is required for gamma_n(r) and N to exist
n			positive integer
nf			positive integer (nf >= n)
gamma/eta	+/-1


integers such that N is prime will be printed to an appropreately named file (A * r^n +/- 1.txt or A * r^n +/- gamma_r.txt) in the res folder
possible prime numbers will have n printed to the terminal window plus some additional constraints found

Alternatively, this program can find prime numbers that may be useful for the SIDH cryptosystem using the same primality testing algorithms.
In this case run

./main bit_length r

to find N such that
N = f * 2^x * r^y - 1 is prime
for integers f, x, y where |r^y| ~= x bits and |N| ~< bit_length bits

Such integers will be printed to an appropreately named file (SIDH-bit_length bit: f*2^x*r^y-1.txt)
The numbers in the file in each line will be f,x,y for a prime N