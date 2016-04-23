Homework #4: The Index Calculus Algorithm
Language/Version: Python 2.7.11
Author: Tien Ho
Date last modified: 23 April 2016

This program implements the index calculus algorithm to solve discrete
logarithms. The algorithm consists of three stages. The first stage searches
for a set of linearly independent smooth relations between the factor base and
the powers of the generator of the finite field. The second stage solves the
linear system of smooth relations using Gauss-Jordan elimination to find
the logarithms of the factor base. The third stage computes the discrete
logarithm of a desired field element.

We test the program using the prime field 14087, the primitive root 5, and the
desired field element 5872. In short, our program seeks to solve for the
discrete problem 5^x == 5872 (mod 14087).

COMMAND:
To run the program, type the command:

python index_calculus_discrete_log.py


OUTPUT:
Running the program will log output to the console. The output lists the three
stages of the index calculus algorithm.

The first stage shows the search for
a list of linearly independent smooth relations between the factor base and the
powers of the generator of the prime field. A linear system that is dependent
is marked as 'LINEARLY DEPENDENT'. The first stage terminates when a set of
linearly independent smooth relations is found.

The second stage solves for the logarithms of the factor base by performing
the Gauss-Jordan elimination on the linearly independent smooth relations.
Since our log smooth relations are the congruences modulo (prime - 1) which is
a ring, we perform the Gauss-Jordan elimination independently for each of the
factors of prime - 1. We then obtain the final solution using the Chinese
Remainder Theorem.

Once the logarithms of the factor base is found, we compute the discrete log
of the desired field element by searching for a power of the generator such
that its product with the desired field element is smooth. Then we apply
the logarithms of the factor base to obtain the final result.

For our example in which to solve for x such that 5^x == 5872 (mod 14087),
the answer to the discrete log of 5872 is 11608.
