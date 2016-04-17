import math

def swap_rows(matrix, row_1, row_2):
    row = matrix[row_1]
    matrix[row_1] = matrix[row_2]
    matrix[row_2] = row
    return matrix

def reduce_matrix(matrix):
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol-1)

    for i in range(npivot):
        divisor = matrix[i][i]
        if divisor != 0:
            for j in range(i, ncol):
                matrix[i][j] = matrix[i][j] / divisor

    return matrix

def finite_field_reduce_matrix(matrix, prime):
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol-1)

    for i in range(npivot):
        ratio = multiplicative_inverse(matrix[i][i], prime)
        if ratio != 0:
            for j in range(i, ncol):
                matrix[i][j] = matrix[i][j]*ratio % prime

    return matrix

def sieve_of_Eratosthenes(B):
    B += 1
    upper_limit = int(math.ceil(math.sqrt(B)))
    primes = range(2, B)
    is_prime = [True] * (B - 2)

    for i in range(2, upper_limit):
        if is_prime[i - 2] == True:
            for j in range(i**2, B, i):
                if is_prime[j - 2] == True:
                    is_prime[j - 2] = False
                    primes.remove(j)

    return primes

def gcd(a, b):
    """
    Euclidean Algorithm

    """
    while b != 0:
        temp = b
        b = a % b
        a = temp
    return a

def multiplicative_inverse(n, p):
    if gcd(n, p) == 1:
        b0 = p
        x0 = 0
        x1 = 1

        if (p == 1):
            return 1

        while (n > 1):
            q = n / p
            rem = n % p
            n = p
            p = rem
            temp = x1 - q * x0
            x1 = x0
            x0 = temp

        if x1 < 0:
            x1 += b0

        return x1

def chinese_remainder(set_of_modulus, set_of_congruences):
    product = 1
    sum = 0

    for m in set_of_modulus:
        product *= m

    for i in range(len(set_of_modulus)):
        n = product / set_of_modulus[i]
        inverse = multiplicative_inverse(n, set_of_modulus[i])
        sum += set_of_congruences[i] * inverse * n

    return sum % product

def gauss_jordan_elimination(matrix):
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol-1)

    for k in range(npivot):
        # swap rows to get a pivot of 1
        for j in range(k, nrow):
            if matrix[j][k] == 1:
                matrix = swap_rows(matrix, k, j)
                break

        # swap rows to get a nonzero pivot
        if matrix[k][k] == 0:
            for j in range(k+1,nrow):
                if matrix[j][k] != 0:
                    matrix = swap_rows(matrix,k,j)
                    break

        if matrix[k][k] != 0:
            for i in range(nrow):
                if i != k and matrix[i][k] != 0:
                    ratio = float(matrix[i][k]) / matrix[k][k]
                    for j in range(k, ncol):
                        matrix[i][j] = matrix[i][j] - ratio*matrix[k][j]

        print_matrix(matrix)

    # make the pivot entries all 1s
    return reduce_matrix(matrix)

def finite_field_gauss_elimination(matrix, prime):
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol-1)

    for k in range(npivot):
        # swap rows to get a pivot of 1
        for j in range(k, nrow):
            if matrix[j][k] == 1:
                matrix = swap_rows(matrix, k, j)
                break

        # swap rows to get a nonzero pivot
        if matrix[k][k] == 0:
            for j in range(k+1,nrow):
                if matrix[j][k] != 0:
                    matrix = swap_rows(matrix,k,j)
                    break

        if matrix[k][k] != 0:
            for i in range(nrow):
                if i != k and matrix[i][k] != 0:
                    ratio = (matrix[i][k] *
                        multiplicative_inverse(matrix[k][k], prime) % prime)
                    for j in range(k, ncol):
                        matrix[i][j] = matrix[i][j] - ratio*matrix[k][j]
                        matrix[i][j] = matrix[i][j] % prime

    return matrix

def solve_for_logarithms(log_relations, prime):
    log_relations = finite_field_gauss_elimination(log_relations, prime)
    log_relations = finite_field_reduce_matrix(log_relations, prime)

    return log_relations

def compute_logarithmic_relations(factor_base, prime, generator):
    log_relations = []
    limit = len(factor_base) + 1
    next_exponent = 1

    # collects initial smooth relations
    log_relations = find_relations(factor_base, prime, generator,
                                   log_relations, limit, next_exponent)

    # keeps collecting more relations until m linearly independent relations
    # are found
    while not is_linearly_independent(log_relations, prime):
        limit += 1
        number_of_relations = len(log_relations)
        number_of_primes = len(factor_base)
        next_exponent = log_relations[number_of_relations - 1][number_of_primes]
        next_exponent += 1
        log_relations = find_relations(factor_base, prime, generator,
                                       log_relations, limit, next_exponent)

    return log_relations

def find_relations(factor_base, prime, generator, log_relations, limit, start):
    for i in range(start, prime):
        relation = [0] * len(factor_base)
        sieving_value = (generator ** i) % prime
        for j in range(len(factor_base)):
            exponents = 0
            while sieving_value % factor_base[j] == 0:
                sieving_value /= factor_base[j]
                exponents += 1
            relation[j] = exponents

        if sieving_value == 1:
            relation.append(i)
            log_relations.append(relation)

        if len(log_relations) == limit:
            break

    return log_relations

def is_linearly_independent(log_relations, prime):
    clone = copy_matrix(log_relations)
    factors = [2, 7043]

    for f in factors:
        matrix = finite_field_gauss_elimination(clone, f)
        ncol = len(matrix[0])

        # if one of the pivots is zero, the columns are not linearly independent.
        for i in range(ncol - 1):
            if matrix[i][i] == 0:
                return False

    return True

def copy_matrix(matrix):
    clone = []
    for row in matrix:
        _row = []
        for value in row:
            _row.append(value)
        clone.append(_row)

    return clone

def find_logarithmic_solutions(log_relations, factor_base):
    solutions = []
    # we assume that the factorization is known
    factors = [2, 7043]

    for f in factors:
        clone = copy_matrix(log_relations)
        clone = solve_for_logarithms(clone, f)

        _solution = []
        for i in range(len(factor_base)):
            _solution.append(clone[i][len(factor_base)])

        solutions.append(_solution)

    return solutions

# compute final solutions using Chinese Remainder Theorem
def combine_solutions(solutions):
    factors = [2, 7043]
    solution = []
    for i in range(len(solutions[0])):
        congruences = []
        for j in range(len(solutions)):
            congruences.append(solutions[j][i])

        s = chinese_remainder(factors, congruences)
        solution.append(s)

    return solution

def compute_individual_log():

    return

def compute_discrete_log(prime, generator, output):
    factor_base = sieve_of_Eratosthenes(15)
    print 'Factor base', factor_base
    log_relations = compute_logarithmic_relations(factor_base, prime, generator)

    # for our purpose, we assume that the factorization of p - 1 is known.
    # We solve the linear system of these smooth relations separately for each
    # prime factor and then use the Chinese Remainder Theorem to obtain the
    # final solution.

    # Our prime p is 14087 --> p - 1 = 14086 = 2 * 7043
    solutions = find_logarithmic_solutions(log_relations, factor_base)
    final_solution = combine_solutions(solutions)

    print final_solution

    # sample = [[1, 1, 0, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0, 0],
    # [1, 1, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0],
    # [1, 0, 1, 1, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0],
    # [0, 1, 0, 0, 1, 0, 0]]

    return

def print_matrix(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print '\n'.join(table) + '\n'
    return

def main():
    generator = 5
    prime = 14087
    output = 20

    # matrix = [[1, 1, 0, 0, 0, 0, 2], [1, 0, 0, 1, 1, 0, 2],
    #             [1, 1, 0, 1, 0, 0, 2],[0, 1, 0, 0, 1, 0, 2],
    #             [1, 0, 1, 1, 0, 1, 2], [0, 0, 0, 1, 0, 0, 0],
    #             [0, 1, 0, 0, 1, 0, 2]]

    compute_discrete_log(prime, generator, output)
    # x = chinese_remainder([5, 7, 11], [3, 5, 7])
    # print x
    return

if __name__ == '__main__':
    main()
