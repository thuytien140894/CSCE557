Homework #2: The Simplified DES Algorithm
Language/Version: Python 2.7.11
Author: Tien Ho
Date last modified: 03 March 2016

This program creates a static class DES that has functions to perform
four rounds of a Feistel System on a 12-bit input bitstring using a 9-bit key.
The class also tests for weak keys for this four-round DES-typed algorithm.

COMMAND:
To run the program in the same directory that contains the ciphertext file,
type the command:

python DES.py

RESULT:
Running the program will log output to a file. The first part of the log file
is the four-round encryption on the input bitstring of 001111100111 and the key
010011001. Below are the ciphertexts for each round:

E_1(001111100111, 01001100) = 100111100101
E_2(100111100101, 10011001) = 100101010011
E_3(100101010011, 00110010) = 010011100101
E_4(010011100101, 01100101) = 100101010101

The second part of the log file is the four-round decryption on the ciphertext
that we obtained from the four-round encryption above. Below are the outputs
for each round:

D_4(100101010101, 01100101) = 010011100101
D_3(010011100101, 00110010) = 100101010011
D_2(100101010011, 10011001) = 100111100101
D_1(100111100101, 01001100) = 001111100111

The output of the four-round decryption returns the original plaintext.

The third part of the log file confirms that there are no weak keys for the
four-round simplified DES algorithm. However, if we modify the algorithm such
that the left and right halves of the ciphertext from the four-round simplified
DES encryption are swapped, then there are two weak keys for this algorithms.
They are 000000000 and 111111111.
