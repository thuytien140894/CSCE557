class DES(object):
    s_1 = [['101', '010', '001', '110', '011', '100', '111', '000'],
           ['001', '100', '110', '010', '000', '111', '101', '011']]
    s_2 = [['100', '000', '110', '101', '111', '001', '011', '010'],
           ['101', '011', '000', '111', '110', '010', '001', '100']]
    log_file = open('logfile', 'w')

    def __init__ (self, bitstring, key):
        self.plaintext = bitstring
        self.key = key
        self.input_for_next_round = bitstring
        self.rotating_key_template = key * 2
        self.key_position = 0
        self.round = 0

    def encrypt(self, number_of_rounds, print_output):
        if (number_of_rounds == 0):
            return self.input_for_next_round

        previous_input = self.input_for_next_round
        previous_left_bits = self.input_for_next_round[0:6]
        previous_right_bits = self.input_for_next_round[6:]
        key_position = self.key_position % 9

        key = self.rotating_key_template[key_position: key_position + 8]
        left_bits = previous_right_bits
        function_output = self.computeFunction(previous_right_bits, key)
        right_bits = self.perform_xor(previous_left_bits, function_output)

        self.input_for_next_round = left_bits + right_bits
        self.key_position += 1
        number_of_rounds -= 1
        self.round += 1

        if (print_output):
            self.print_encryption(self.round, previous_input, key,
                                  function_output, self.input_for_next_round)

        return self.encrypt(number_of_rounds, print_output)

    def decrypt(self, output, number_of_rounds):
        if (number_of_rounds == 0):
            return output

        left_bits = output[0:6]
        right_bits = output[6:]
        key_position = number_of_rounds - 1
        key = self.rotating_key_template[key_position: key_position + 8]
        function_output = self.computeFunction(left_bits, key)
        previous_left_bits = self.perform_xor(right_bits, function_output)
        previous_right_bits = left_bits

        previous_output = previous_left_bits + previous_right_bits

        self.print_decryption(number_of_rounds, output, key, function_output,
                              previous_output)

        number_of_rounds -= 1

        return self.decrypt(previous_output, number_of_rounds)

    def expand(self, bitstring):
        newBitstring = [bitstring[0:2], bitstring[3], bitstring[2],
                        bitstring[3], bitstring[2], bitstring[4:]]
        return ''.join(newBitstring)

    def computeFunction(self, right_bits, key):
        expandedBitstring = self.expand(right_bits)

        s_input = self.perform_xor(expandedBitstring, key)
        s_1_input = s_input[0:4]
        s_2_input = s_input[4:]
        s_1_output = DES.s_1[int(s_1_input[0])][int(s_1_input[1:], 2)]
        s_2_output = DES.s_2[int(s_2_input[0])][int(s_2_input[1:], 2)]

        output = s_1_output + s_2_output

        return output

    def perform_xor(self, a, b):
        result = []
        for x, y in zip(a, b):
            result.append(str(int(x) ^ int(y)))

        return ''.join(result)

    def encrypt_four_rounds(self, bitstring, key):
        self.input_for_next_round = bitstring
        self.rotating_key_template = key * 2
        self.key_position = 0

        return self.encrypt(4, False)

    def encrypt_four_rounds_and_swap(self, bitstring, key):
        output = self.encrypt_four_rounds(bitstring, key)

        return self.swap(output)

    def swap(self, bitstring):
        swapped_bitstring = bitstring[6:] + bitstring[0:6]

        return swapped_bitstring

    def generate_all_possible_keys(self):
        possible_keys = ['0', '1']

        for i in range(8):
            temp = []
            for j in range(len(possible_keys)):
                temp.append(possible_keys[j] + '0')
                temp.append(possible_keys[j] + '1')

            possible_keys = []
            possible_keys += temp

        return possible_keys

    def generate_all_possible_input_bitstrings(self):
        possible_input_bitstrings = ['0', '1']

        for i in range(11):
            temp = []
            for j in range(len(possible_input_bitstrings)):
                temp.append(possible_input_bitstrings[j] + '0')
                temp.append(possible_input_bitstrings[j] + '1')

            possible_input_bitstrings = []
            possible_input_bitstrings += temp

        return possible_input_bitstrings

    def test_for_weak_key(self, output_is_swapped):
        weak_keys = []
        all_possible_keys = self.generate_all_possible_keys()
        all_possible_input_bitstrings = (
                self.generate_all_possible_input_bitstrings())

        for i in range(len(all_possible_keys)):
            key = all_possible_keys[i]
            for j in range(len(all_possible_input_bitstrings)):
                bitstring = all_possible_input_bitstrings[j]
                if (not output_is_swapped):
                    output = self.encrypt_four_rounds(
                        self.encrypt_four_rounds(bitstring, key), key)
                else:
                    output = self.encrypt_four_rounds_and_swap(
                        self.encrypt_four_rounds_and_swap(bitstring, key), key)

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

    def print_encryption(self, round, input, key, function_output, output):
        previous_left = input[0:6]
        previous_right = input[6:]
        output_left = output[0:6]
        output_right = output[6:]

        DES.log_file.write('Round %s \n\n' % (str(round)))
        DES.log_file.write('Input Bitstring: %s \n\n' % (input))
        DES.log_file.write('L_%s = %s \n' % (str(round - 1), previous_left))
        DES.log_file.write('R_%s = %s \n\n' % (str(round - 1), previous_right))
        DES.log_file.write('L_%s = R_%s = %s \n' %
                            (str(round), str(round - 1), output_left))
        DES.log_file.write('R_%s = L_%s ^ f(R_%s,K_%s) = %s ^ f(%s,%s)'
                           ' = %s ^ %s = %s \n\n'
                           % (str(round), str(round - 1), str(round - 1),
                           str(round), previous_left, previous_right, key,
                           previous_left, function_output, output_right))

        DES.log_file.write('Output Bitstring: %s \n\n' % (output))
        divider = '*' * 80
        DES.log_file.write(divider)
        DES.log_file.write('\n\n')

        return

    def print_decryption(self, round, input, key, function_output, output):
        left = input[0:6]
        right = input[6:]
        previous_left = output[0:6]
        previous_right = output[6:]

        DES.log_file.write('Round %s \n\n' % (str(round)))
        DES.log_file.write('Output Bitstring: %s \n\n' % (input))
        DES.log_file.write('L_%s = %s \n' % (str(round), left))
        DES.log_file.write('R_%s = %s \n\n' % (str(round), right))
        DES.log_file.write('L_%s = R_%s ^ f(L_%s,K_%s) = %s ^ f(%s,%s) '
                           '= %s ^ %s = %s \n'
                           % (str(round - 1), str(round), str(round),
                           str(round), right, left, key, right,
                           function_output, previous_left))
        DES.log_file.write('R_%s = %s \n\n' % (str(round - 1), previous_right))
        DES.log_file.write('Input Bitstring: %s \n\n' % (output))

        divider = '*' * 80
        DES.log_file.write(divider)
        DES.log_file.write('\n\n')

        return

my_DES = DES('001111100111', '010011001')
DES.log_file.write('FOUR-ROUND ENCRYPTION \n\n')
output = my_DES.encrypt(4, True)
DES.log_file.write('FOUR-ROUND DECRYPTION \n\n')
input_bitstring = my_DES.decrypt(output, 4)
DES.log_file.write('Testing for weak keys for the four-round '
                   'simplified DES encryption \n')
my_DES.test_for_weak_key(False)
divider = '*' * 80
DES.log_file.write(divider)
DES.log_file.write('\n\n')
DES.log_file.write('Testing for weak keys for the swapped four-round '
                   'simplified DES encryption \n')
my_DES.test_for_weak_key(True)
