import collections
import sys

"""
Program to decrypt a ciphertext encrypted using the Vigenere cipher.

This program contains functions to find the key length, find the key used to
encrypted the plaintext, and decrypt the cipher_text using the determined key.
The main method reads in the file that contains the cipher_text and invokes
all these methods to output the dycrypted text.

Author/copyright: Tien Ho. All rights reserved.
Date last modified: 04 February 2015
"""

# Global constants
_FREQUENCIES_OF_LETTERS = [0.082, 0.015, 0.028, 0.043, 0.127, 0.022,
                           0.020, 0.061, 0.070, 0.002, 0.008, 0.040,
                           0.024, 0.067, 0.075, 0.019, 0.001, 0.060,
                           0.063, 0.091, 0.028, 0.010, 0.023, 0.001,
                           0.020, 0.001]
_ROTATING_TEMPLATE = _FREQUENCIES_OF_LETTERS * 2
_HIGHEST_POSSIBLE_KEY_LENGTH = 10
_OFFSET = 97
_NUMBER_OF_ALPHABET_LETTERS = 26

def perform_dot_product (a, b):
	"""
	Function 'perform_dot_product' computes the dot product of two vectors.

	Parameters:
		a, b: Two lists of the same lengths which represent two vectors

	Output:
		A sum of all the products of the two corresponding elements in the
	    two lists

	Returns:
		The dot product of the input lists.
	"""
	sum = 0.0
	for x, y in zip(a, b):
        sum += (x*y)

	return sum

def find_key_length(text):
	"""
	Function 'find_key_length' determines the length of the key used to encrypt
	the plaintext.

	This function compares the ciphertext itself with a n-position right shift
	of the ciphertext and records the number of coincidences (where the letters
	between the two texts agree) for each shift of n positions. The shift of n
	positions with the highest number of coincidences is the best guess for the
	length of key.

	At first, I tried to determine all the possible shifts up to the entire
	length of the ciphertext and obtained either 20 or 35 as the possible key
	lengths because they have the highest number of coincidences, which is 28.
	However, finding the key using these key lengths made me realize that they
	just resulted in the key with repeated 5-letter keys.  Therefore, I capped
	the number of shifts that this function should run up to only 10.

	Parameter:
		text: A string that is the ciphertext

	Ouput:
		The possible key length

	Returns:
		A non-negative integer that is the key length
	"""
	list_of_displaced_text = []
	number_of_coincidences = []
	key_length = 0

	# create a list of cipher_texts displaced up to 10 positions to the right
	# of the original cipher_text
	for i in range(1, _HIGHEST_POSSIBLE_KEY_LENGTH):
		space = ' ' * i
		new_text = space + text
		list_of_displaced_text.append(new_text)

	# record the number of coincidences between each displaced cipher_text and
	# the cipher_text itself
	for i in range(len(list_of_displaced_text)):
		count = 0
		# iterate through both the texts at the same time to access their
		# corresponding characters.  The iterator stops when the shortest
		# string (in this case the original text) exhausts.
		for j, k in zip(list_of_displaced_text[i], text):
			if j == k:
				count += 1

		number_of_coincidences.append(count)

	# Each index n in number_of_coincidences corresponds to an n + 1 shift
	key_length = number_of_coincidences.index(max(number_of_coincidences)) + 1

	return key_length

def find_key(text, key_length):
	"""
	Function 'find_key' determines the key used for the encryption of the
	plaintext.

	This function uses the second method described in the textbook.  Basically,
	using the possible key length that we find, in order to find the nth
	element of the key, we first determine the occurrences of the letters at
	the positions in the ciphertext where this nth element is supposed to line
	up for ecryption.  The letters at these positions represent a random sample
	of English letters all shifted by the same amount.  In order to determine
	how much of this shift, we compute the dot product between the determined
	vector of letter occurrences and the m-position cyclic right shift vector
	containing the frequencies of English letters (found in the textbook on
	page 17) for n = 0, 2, ..., 25.  The m shift that results in the highest
	dot product with the vector of letter occurrences at the nth element of
	the key suggests the best guess for the value of that element.  We repeat
	this algorithm to find all the key elements.

	Parameters:
		text: A string that is the ciphertext
		key_length: A positive integer that is the key length

	Output:
		The key

	Returns:
		A list that contains the numbers of shifts at each element of the key
	"""
	key = []
	alphabet = collections.OrderedDict([('a', 0), ('b', 0), ('c', 0), ('d', 0),
                                        ('e', 0), ('f', 0), ('g', 0), ('h', 0),
                                        ('i', 0), ('j', 0), ('k', 0), ('l', 0),
                                        ('m', 0), ('n', 0), ('o', 0), ('p', 0),
                                        ('q', 0), ('r', 0), ('s', 0), ('t', 0),
                                        ('u', 0), ('v', 0), ('w', 0), ('x', 0),
                                        ('y', 0), ('z', 0)])

	for j in range(key_length):
		# reset the dot_products and the alphabet frequency
		dot_products = []
		alphabet = alphabet.fromkeys(alphabet, 0)

		# count the frequency of letters at positions in the cipher_text
		# that line up with the jth element of the key
		number_of_letters_counted = 0
		for i in range(j, len(text), key_length):
			number_of_letters_counted += 1
			alphabet[text[i]] += 1

		vector_of_occurrences = ([x / float(number_of_letters_counted)
                                  for x in alphabet.values()])

		for i in range(_NUMBER_OF_ALPHABET_LETTERS):
			# in order to perform cyclic right shift of the vector of English
			# letter frequencies for i positions, we use rotating template
			# that contains this frequency vectors twice and take the
			# appropriate sublist.
			start_position = len(_FREQUENCIES_OF_LETTERS) - i
			end_position = len(_ROTATING_TEMPLATE) - i
			shiftedvector_of_occurrences = (_ROTATING_TEMPLATE[start_position:
                                                               end_position])

			dot_product = perform_dot_product(vector_of_occurrences,
                                              shiftedvector_of_occurrences)
			dot_products.append(dot_product)

		print dot_products
		key.append(dot_products.index(max(dot_products)))

	return key

def convert_key(key):
	"""
	Function 'convert_key' converts the key digits into text.

	This function iterates over all the key digits found in the method
	'find_key', converts them to the appropriate Unicode of the characters
	they represent.

	Parameter:
		key: A list containing the key digits

	Ouput:
		The key phrase used for encryption

	Returns:
		A string that is the key phrase
	"""
	key_letters = [chr(n + _OFFSET) for n in key]
	key_text = ''.join(key_letters)

	return key_text

def decrypt(text, key):
	"""
	Function 'decrypt' decrypts the ciphertext using the found key.

	This function loops through all the letters in the ciphertext, converts
	them to an integer that corresponds to their position in the alphabet by
	subtracting 97 from their Unicode point), subtracting this number by the
	number of shifts caused by the key element lined up with the letter at its
	position, takes the modulus of this difference (in case of negativity)
	to get the alphabet position of the decrypted letter, adds back 97 to get
	the correct Unicode of the letter.

	Parameters:
		text: A string that is the ciphertext
		key: A list that contains the elements of the key

	Output:
		The decrypted text

	Returns:
		A string that is the decrypted text
	"""
	decrypted_letters = []
	decrypted_text = ''
	key_length = len(key)
	for i in range(len(text)):
		letter_number = ord(text[i]) - _OFFSET
		key_letter_number = key[i % key_length]
		decrypted_letter = chr((letter_number - key_letter_number) %
		                        _NUMBER_OF_ALPHABET_LETTERS + _OFFSET)
		decrypted_letters.append(decrypted_letter)

	decrypted_text = ''.join(decrypted_letters)

	return decrypted_text

def print_output(cipher_text, key_text, decrypted_text):
	"""
	Function 'print_output' outputs the ciphertext, the key used for
	encryption, and the decrypted text.

	Parameters:
		cipher_text: A string that is the ciphertext
        key_text: A string that is the key
		decrypted_text: A string that is the decrypted text

	Output:
		The ciphertext
		The key
		The decrypted text

	Returns:
		None
	"""
	print 'Ciphertext: \n' + cipher_text + '\n'
	print 'Key Used for Encryption: \n' + key_text + '\n'
	print 'Decrypted Text: \n' + decrypted_text + '\n'

def main():
	"""
	Function 'main' to perform invoke all the function necessary to decrypt
	the ciphertext file read from the command line.

	Parameters:
		None

	Output:
		The ciphertext
		The key
		The decrypted text

	Returns:
		None
	"""
	with open(sys.argv[1], 'r') as file:
		cipher_text = file.read()
		cipher_text = cipher_text.replace('\n', '')

	key_length = find_key_length(cipher_text);
	key = find_key(cipher_text, key_length)
	key_text = convert_key(key)
	decrypted_text = decrypt(cipher_text, key)

	print_output(cipher_text, key_text, decrypted_text)

	return

# make sure that the main() method is only executed after the module
# gets imported
if __name__ == '__main__':
	# make sure the ciphertext file is included in the command line
	if len(sys.argv) != 2:
		print 'Missing ciphertext'
		sys.exit()

	main()
