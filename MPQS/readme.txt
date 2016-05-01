Homework #3: The Quadratic Sieve Algorithm
Language/Version: Python 2.7.11
Author: Steven Dao, Tien Ho, and Brandon Hostetter
Date last modified: 25 April 2016

This program factors a number using the quadratic sieve algorithm.
Given that N is the number we would like to factor, R is the square root of N,
this program first uses the function (R + x)^2 - N to sieve through values of x
in a chosen interval and find smooth numbers that factor complete over a factor
base. We build a matrix whose rows contain the exponents of the primes in the
factor base for each of the smooth numbers. Then we find a combination of rows
that would yield a perfect square by performing the Gauss-Jordan elimination.
We use these rows to compute x and y in the relation x^2 == y^2 (mod N).
Taking the gcd(x - y, N) or the gcd(x + y, N) would give a factor of N.

We tested our program with the number 239812014798221, which gives us two
factors 15485863 and 15485867.

COMMAND:
To run the program, type the command:

python MPQS.py


OUTPUT:
Running the program will log output to the console. We put the output into a
text file for easy reference.
