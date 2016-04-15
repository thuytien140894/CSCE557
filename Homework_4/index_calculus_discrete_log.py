import math

def swap_rows(matrix, row_1, row_2):
    row = matrix[row_1]
    matrix[row_1] = matrix[row_2]
    matrix[row_2] = row
    return matrix

def reduce_matrix(matrix):
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol)

    for i in range(npivot):
        divisor = matrix[i][i]
        if divisor != 0:
            for j in range(i, ncol):
                matrix[i][j] = matrix[i][j] / divisor

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

def gauss_jordan_elimination(matrix):
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol)

    for k in range(npivot):
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

    # make the pivot entries all 1s
    return reduce_matrix(matrix)

def solve_for_logarithms():
    return

def find_logarithmic_relations(factor_base, prime, generator):
    sieving_values = [(generator ** i) % prime for i in range(prime)]
    log_relations = []
    limit = len(factor_base) + 4

    for i in range(len(sieving_values)):
        relation = [0] * len(factor_base)
        for j in range(len(factor_base)):
            exponents = 0
            while sieving_values[i] % factor_base[j] == 0:
                sieving_values[i] = sieving_values[i] / factor_base[j]
                exponents += 1
            relation[j] = exponents

        if sieving_values[i] == 1:
            relation.append(i)
            log_relations.append(relation)

        if len(log_relations) == limit:
            break

    return log_relations

def is_linearly_independent(new_relation, existing_relations):
    if len(existing_relations) == 0:
        return True


    return

def compute_discrete_log(prime, generator, output):
    factor_base = sieve_of_Eratosthenes(10)
    log_relations = find_logarithmic_relations(factor_base, prime, generator)

    print log_relations

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
    prime = 503
    output = 20

    # x = compute_discrete_log(prime, generator, output)

    # l = sieve_of_Eratosthenes(1000)
    # print l
    #
    # matrix = [[1, 1, 0, 0, 0, 0, 2], [1, 0, 0, 1, 1, 0, 2],
    #             [1, 1, 0, 1, 0, 0, 2],[0, 1, 0, 0, 1, 0, 2],
    #             [1, 0, 1, 1, 0, 1, 2], [0, 0, 0, 1, 0, 0, 0],
    #             [0, 1, 0, 0, 1, 0, 2]]
    #
    # matrix = gauss_jordan_elimination(matrix)
    # print_matrix(matrix)

    compute_discrete_log(prime, generator, output)

    return

if __name__ == '__main__':
    main()
