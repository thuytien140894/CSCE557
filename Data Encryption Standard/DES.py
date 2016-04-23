"""
Program to implement the simplied DES-type algorithm for encryption
and decryption.

This program creates a static class DES that has functions to perform
four rounds of a Feistel System on a 12-bit input bitstring using a 9-bit key.
The class also tests for weak keys for this four-round DES-typed algorithm.

Author/copyright: Tien Ho. All rights reserved.
Date last modified: 03 March 2016
"""
class DES(object):
    """
    Class 'DES' has class methods to encrypt and decrypt using the simplified
    DES algorithm. The class also tests for weak keys in this algorithm.
    """
    S_1 = [['101', '010', '001', '110', '011', '100', '111', '000'],
           ['001', '100', '110', '010', '000', '111', '101', '011']]
    S_2 = [['100', '000', '110', '101', '111', '001', '011', '010'],
           ['101', '011', '000', '111', '110', '010', '001', '100']]
    log_file = open('logfile', 'w')
    input_for_next_round = ''
    rotating_key_template = ''
    key_position = 0
    round_number = 0
    key_length = 0
    bitstring_length = 0

    @classmethod
    def initialize(cls, bitstring, key):
        """
        Function 'initialize' initializes the class variables in order to
        perform an encryption of a different bitstring and key.

        Parameters:
            bitstring: a 12-bit string of 1s and 0s
            key: a 9-bit string of 1s and 0s

        Returns:
            None
        """
        cls.input_for_next_round = bitstring
        cls.rotating_key_template = key * 2
        cls.key_position = 0
        cls.round_number = 0
        cls.key_length = len(key)
        cls.bitstring_length = len(bitstring)

    @classmethod
    def encrypt(cls, number_of_rounds, print_output):
        """
        Function 'encrypt' performs an n-round encryption on an input bitstring
        using the DES-typed algorithm.

        Assuming that the input bitstring is represented as L(i)R(i) and the
        output ciphertext is represented as L(r+1)R(i+1). This function first
        assigns R(i) to L(i+1). To determine R(i+1), we XOR L(i) with
        f(R(i), K(i+1)). This f function first expands R(i) to an 8-bit
        bitstring, which is then XORed with the key used for that round. Given
        the 9-bit key for the whole encryption, K(i+1) used
        specifically for round (i+1) is determined by using 8 bits of the key,
        starting with the (i+1) bit and wrapping around the key if necessarily.
        The XOR output from the expanded R(i) and K(i+1) is then fed into the
        s boxes to compute its corresponding 6-bit output.

        Depending on the number of rounds that we would like to encrypt, the
        output of each round is then fed to the next round. This
        recursive function terminates when we have performed all the rounds.

        Parameters:
            number_of_rounds: an integer indicating the number of rounds to
                              perform encryption
            print_output: a boolean indicating whether we would like to print
                          out the steps and output of each encryption round.

        Returns:
            a 12-bit string representing the ciphertext
        """
        if (number_of_rounds == 0):
            return cls.input_for_next_round

        half_length = cls.bitstring_length / 2
        previous_input = cls.input_for_next_round
        previous_left_bits = cls.input_for_next_round[0:half_length]
        previous_right_bits = cls.input_for_next_round[half_length:]
        key_position = cls.key_position % cls.key_length

        key = cls.rotating_key_template[key_position: key_position +
                                        (cls.key_length - 1)]
        left_bits = previous_right_bits
        function_output = cls.compute_function(previous_right_bits, key)
        right_bits = cls.perform_xor(previous_left_bits, function_output)

        cls.input_for_next_round = left_bits + right_bits
        cls.key_position += 1
        number_of_rounds -= 1
        cls.round_number += 1

        if (print_output):
            cls.print_encryption(cls.round_number, previous_input, key,
                                 function_output, cls.input_for_next_round)

        return cls.encrypt(number_of_rounds, print_output)

    @classmethod
    def decrypt(cls, output, number_of_rounds):
        """
        Function 'decrypt' performs an n-round decryption on a ciphertext
        using the DES-typed algorithm.

        Assuming that the ciphertext is represented as L(i)R(i) and the
        input ciphertext of the previous round is represented as L(r-1)R(i-1).
        This function first assigns L(i) to R(i-1). To determine L(i-1), we XOR
        R(i) with f(L(i), K(i)). This f function first expands L(i) to an 8-bit
        bitstring, which is then XORed with the key used for that round. Given
        the 9-bit key for the whole encryption, K(i) used
        specifically for ith round is determined by using 8 bits of the key,
        starting with the ith bit and wrapping around the key if necessarily.
        The XOR output from the expanded L(i) and K(i) is then fed into the
        s boxes to compute its corresponding 6-bit output.

        Depending on the number of rounds that we would like to decrypt, the
        output of each round is then fed to the next round. This
        recursive function terminates when we have performed all the rounds.

        Parameters:
            output: a 12-bit string
            number_of_rounds: an integer indicating the number of rounds to
                              perform decryption

        Returns:
            a 12-bit string that is the plaintext
        """
        if (number_of_rounds == 0):
            return output

        half_length = cls.bitstring_length / 2
        left_bits = output[0:half_length]
        right_bits = output[half_length:]
        key_position = number_of_rounds - 1
        key = cls.rotating_key_template[key_position: key_position +
                                        (cls.key_length - 1)]
        function_output = cls.compute_function(left_bits, key)
        previous_left_bits = cls.perform_xor(right_bits, function_output)
        previous_right_bits = left_bits

        previous_output = previous_left_bits + previous_right_bits

        cls.print_decryption(number_of_rounds, output, key, function_output,
                             previous_output)

        number_of_rounds -= 1

        return cls.decrypt(previous_output, number_of_rounds)

    @classmethod
    def expand(cls, bitstring):
        """
        Function 'expand' expands a 6-bit bitstring into an 8-bit bitstring.

        The new expanded bitstring has its third and fifth bit set to the fourth
        bit of the input bitstring and its fourth and sixth bit set to the
        third bit of the input bitstring. All the other bits are kept the same.

        Parameters:
            bitstring: a 6-bit bitstring

        Returns:
            an 8-bit bitstring
        """
        newBitstring = [bitstring[0:2], bitstring[3], bitstring[2],
                        bitstring[3], bitstring[2], bitstring[4:]]
        return ''.join(newBitstring)

    @classmethod
    def compute_function(cls, right_bits, key):
        """
        Function 'compute_function' computes f(R(i-1), K(i)) used to determine
        the second half of the ciphertext.

        This function first expands the right half of the input bitstring,
        and XORs it with the key. This XOR output is then fed to the s boxes
        to find the corresponding 6-bit bitstring.

        Parameters:
            right_bits: the 6-bit right half of an input string
            key: an 8-bit bitstring

        Returns:
            a 6-bit bitstring
        """
        expandedBitstring = cls.expand(right_bits)

        s_input = cls.perform_xor(expandedBitstring, key)
        s_1_input = s_input[0:4]
        s_2_input = s_input[4:]
        s_1_output = cls.S_1[int(s_1_input[0])][int(s_1_input[1:], 2)]
        s_2_output = cls.S_2[int(s_2_input[0])][int(s_2_input[1:], 2)]

        output = s_1_output + s_2_output

        return output

    @classmethod
    def perform_xor(cls, a, b):
        """
        Function 'perform_xor' XORs two input bitstrings of the same length.

        Parameters:
            a, b: two bitstrings

        Returns:
            a bitstring
        """
        result = []
        for x, y in zip(a, b):
            result.append(str(int(x) ^ int(y)))

        return ''.join(result)

    @classmethod
    def encrypt_four_rounds(cls, bitstring, key):
        """
        Function 'encrypt_four_rounds' performs a four-round DES encryption.

        Parameters:
            bitstring: a 12-bit bitstring
            key: a 9-bit bitstrings

        Returns:
            a 12-bit bitstring that is the ciphertext
        """
        cls.initialize(bitstring, key)

        return cls.encrypt(4, False)

    @classmethod
    def encrypt_four_rounds_and_swap(cls, bitstring, key):
        """
        Function 'encrypt_four_rounds_and_swap' swaps the left and right halves
        after the four Feistel encryption rounds.

        Paramters:
            bitstring: a 12-bit bitstring
            key: a 9-bit bitstring

        Returns:
            a 12-bit bitstring that is the ciphertext
        """
        output = cls.encrypt_four_rounds(bitstring, key)

        return cls.swap(output)

    @classmethod
    def swap(cls, bitstring):
        """
        Function 'swap' swaps the left half and the right half of a bitstring.

        Parameters:
            bitstring: a 12-bit bitstring

        Returns:
            a 12-bit bitstring whose left and right halves are swapped
        """
        half_length = len(bitstring) / 2
        swapped_bitstring = bitstring[half_length:] + bitstring[0:half_length]

        return swapped_bitstring

    @classmethod
    def generate_all_possible_bitstrings(cls, number_of_bits):
        """
        Function 'generate_all_possible_bitstrings' generates all possible
        bitstrings given its length in bits.

        Parameters:
            number_of_bits: an integer representing the number of bits
                            to generate all possible bitstrings

        Returns:
            a list of all possible n-bit bitstrings
        """
        possible_bitstrings = ['0', '1']

        for i in range(number_of_bits - 1):
            temp = []
            for j in range(len(possible_bitstrings)):
                temp.append(possible_bitstrings[j] + '0')
                temp.append(possible_bitstrings[j] + '1')

            possible_bitstrings = []
            possible_bitstrings += temp

        return possible_bitstrings

    @classmethod
    def generate_all_possible_input_bitstrings(cls):
        possible_input_bitstrings = ['0', '1']

        for i in range(11):
            temp = []
            for j in range(len(possible_input_bitstrings)):
                temp.append(possible_input_bitstrings[j] + '0')
                temp.append(possible_input_bitstrings[j] + '1')

            possible_input_bitstrings = []
            possible_input_bitstrings += temp

        return possible_input_bitstrings

    @classmethod
    def test_for_weak_key(cls, output_is_swapped):
        """
        Function 'test_for_weak_key' checks for weak keys for the simplified
        DES-type algorithm implemented in two different ways. The first one
        is just a four-round Feistel encryption. The second one is similar to
        the first one, but the left and right halves are swapped after
        the four Feistel rounds.

        The function first generates all possible 9-bit keys and 12-bit
        input bitstrings. The function then doubles encrypts each key with an
        input bitstring. A key is considered weak when its double encryption
        with every possible bitstrings returns the bitstring itself. In other
        word, as long as there is some bitstring M such that E(K)(E(K)(M))
        does not give M back, the key K is not weak.

        Parameters:
            output_is_swapped: a boolean indicating whether the ciphertext
                               after the four-round Feistel encryption is
                               swapped. This boolean tells the function which
                               of the two implementations aforesaid should
                               be tested for weak keys.

        Returns:
            None
        """
        weak_keys = []
        all_possible_keys = cls.generate_all_possible_bitstrings(cls.key_length)
        all_possible_input_bitstrings = (
                cls.generate_all_possible_bitstrings(cls.bitstring_length))

        for i in range(len(all_possible_keys)):
            key = all_possible_keys[i]
            for j in range(len(all_possible_input_bitstrings)):
                bitstring = all_possible_input_bitstrings[j]
                if (not output_is_swapped):
                    output = cls.encrypt_four_rounds(
                        cls.encrypt_four_rounds(bitstring, key), key)
                else:
                    output = cls.encrypt_four_rounds_and_swap(
                        cls.encrypt_four_rounds_and_swap(bitstring, key), key)

                if (output != bitstring):
                    break

                if (j == len(all_possible_input_bitstrings) - 1):
                    weak_keys.append(key)

        if (len(weak_keys) == 0):
            DES.log_file.write('--> There are no weak keys. \n\n')
        else:
            DES.log_file.write('--> The weak keys are: ')
            for i in range(len(weak_keys)):
                DES.log_file.write(weak_keys[i] + ' ')

        return

    @classmethod
    def print_encryption(cls, round, input, key, function_output, output):
        """
        Function 'print_encryption' prints the steps and output of each
        ecryption round to a log file.

        Parameters:
            round: an integer indicating the round number
            input: the 12-bit input bitstream for each round
            key: the 8-bit key used for each round
            function_output: the 6-bit bitstring resulted from f(R(i-1), K(i))
            output: the 12-bit output bitstring after each round

        Returns:
            None
        """
        previous_left = input[0:6]
        previous_right = input[6:]
        output_left = output[0:6]
        output_right = output[6:]

        cls.log_file.write('ROUND %s \n\n' % (str(round)))
        cls.log_file.write('Input Bitstring: %s \n\n' % (input))
        cls.log_file.write('L_%s = %s          '
                           % (str(round - 1), previous_left))
        cls.log_file.write('R_%s = %s \n\n' % (str(round - 1), previous_right))
        cls.log_file.write('L_%s = R_%s \n'
                           % (str(round), str(round - 1)))
        cls.log_file.write('    = %s \n\n' % (output_left))
        cls.log_file.write('R_%s = L_%s ^ f(R_%s, K_%s) \n'
                           % (str(round), str(round - 1), str(round - 1),
                           str(round)))
        cls.log_file.write('    = %s ^ f(%s,%s) \n'
                           % (previous_left, previous_right, key))
        cls.log_file.write('    = %s ^ %s \n'
                           % (previous_left, function_output))
        cls.log_file.write('    = %s \n\n' % (output_right))

        cls.log_file.write('Output Bitstring: %s \n\n' % (output))
        divider = '-' * 80
        cls.log_file.write(divider)
        cls.log_file.write('\n\n')

        return

    @classmethod
    def print_decryption(cls, round, input, key, function_output, output):
        """
        Function 'print_decryption' prints the steps and output of each
        decryption round to a log file.

        Parameters:
            round: an integer indicating the round number
            input: the 12-bit input bitstream for each round
            key: the 8-bit key used for each round
            function_output: the 6-bit bitstring resulted from f(L(i), K(i))
            output: the 12-bit output bitstring after each round

        Returns:
            None
        """
        left = input[0:6]
        right = input[6:]
        previous_left = output[0:6]
        previous_right = output[6:]

        cls.log_file.write('ROUND %s \n\n' % (str(round)))
        cls.log_file.write('Output Bitstring: %s \n\n' % (input))
        cls.log_file.write('L_%s = %s          ' % (str(round), left))
        cls.log_file.write('R_%s = %s \n\n' % (str(round), right))
        cls.log_file.write('L_%s = R_%s ^ f(L_%s, K_%s) \n'
                           % (str(round - 1), str(round), str(round),
                           str(round)))
        cls.log_file.write('    = %s ^ f(%s,%s) \n' % (right, left, key))
        cls.log_file.write('    = %s ^ %s \n' % (right, function_output))
        cls.log_file.write('    = %s \n\n' % (previous_left))
        cls.log_file.write('R_%s = L_%s \n' % (str(round - 1), str(round)))
        cls.log_file.write('    = %s \n\n' % (previous_right))
        cls.log_file.write('Input Bitstring: %s \n\n' % (output))

        divider = '-' * 80
        cls.log_file.write(divider)
        cls.log_file.write('\n\n')

        return

def main():
    """
    Function 'main' uses the DES class to encrypt and decrypt a sample input
    bitstring and a key for four Feistel rounds. We also test for weak keys,
    first for the usual four-round encryption, and secondly for the four-round
    encryption in which the left and right halves of the output are swapped.

    Parameters:
        None

    Returns:
        None
    """
    divider = '*' * 80

    DES.initialize('001111100111', '010011001')
    DES.log_file.write(divider + '\n\n')
    DES.log_file.write('FOUR-ROUND ENCRYPTION \n\n')
    DES.log_file.write(divider + '\n\n')
    output = DES.encrypt(4, True)

    DES.log_file.write(divider + '\n\n')
    DES.log_file.write('FOUR-ROUND DECRYPTION \n\n')
    DES.log_file.write(divider + '\n\n')
    input_bitstring = DES.decrypt(output, 4)

    DES.log_file.write(divider + '\n\n')
    DES.log_file.write('TESTING FOR WEAK KEYS FOR THE FOUR-ROUND '
                       'SIMPLIFIED DES DECRYPTION \n\n')
    DES.log_file.write(divider + '\n\n')
    DES.test_for_weak_key(False)

    DES.log_file.write(divider)
    DES.log_file.write('\n\n')
    DES.log_file.write('TESTING FOR WEAK KEYS FOR THE SWAPPED FOUR-ROUND '
                       'SIMPLIFIED DES ENCRYPTION \n\n')
    DES.log_file.write(divider + '\n\n')
    DES.test_for_weak_key(True)

    return

if __name__ == '__main__':
    main()
