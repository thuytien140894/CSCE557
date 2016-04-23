import math

"""
Program to compute discrete logarithms using the index calculus algorithm.

This program performs three stages of the index calculus algorithm. The first
stage consists of searching for a set of linearly independent smooth relations
for the logarithms of the factor base. Each relation results from factoring
a power of the generator of the finite field completely over the factor base.
The second stage solves the linear systems of the logarithms using Gauss
elimination. The third stage computes the logarithm of a desired field element. 

Author/copyright: Tien Ho. All rights reserved.
Date last modified: 23 April 2016
"""

def print_matrix(matrix, prime):
    """
    Function to print a formatted matrix
    """
    for row in matrix:
        print '[{0}]'.format(''.join(['{0:^7}'.format(str(row[i]))
                                     for i in range(len(row))])), \
              '(mod {})'.format(prime)

    return

def print_matrix_with_headers(matrix, generator, factor_base, prime):
    """
    Function to print a matrix with the headers for the rows and columns
    """
    clone = copy_matrix(matrix)
    print '{0:>6}'.format(''), ''.join(['{0:>7}'.format('[{0:}]'.format(str(f)))
                                        for f in factor_base])

    for row in clone:
        print '{0:>7}'.format('{}^{}'.format(generator, row[len(row)-1])), \
              '[{0}]'.format(''.join(['{0:^7}'.format(str(row[i]))
                                     for i in range(len(row)-1)])), \
              '(mod {})'.format(prime)

    return

def print_augmented_matrix(matrix, generator, factor_base, prime):
    """
    Function to print an augmented matrix used to perform Gauss elimination
    """
    clone = copy_matrix(matrix)
    print '{0:>6}'.format(''), \
          ''.join(['{0:>7}'.format('L({0:})'.format(str(f)))
                  for f in factor_base]),\
          '{0:>6}'.format('exp')

    for row in clone:
        print '{0:>7}'.format('{}^{}'.format(generator, row[len(row)-1])), \
              '[{0}]'.format(''.join(['{0:^7}'.format(str(row[i]))
                                     for i in range(len(row))])), \
              '(mod {})'.format(prime - 1)

    return

def print_log_solution(solution, factor_base):
    """
    Function to print the solution to the linear system of the logarithms of
    the factor base
    """
    for a, b in zip(factor_base, solution):
        print 'L({}) = {}'.format(a, b)
    print

    return

def print_solution(factor_base, log_coefficients, solution, b, g, s, p, result):
    """
    Function to print the steps in the last stage to find the logarithm of the
    field element
    """
    print 'BRUTE FORCE SEARCH FINDS THAT'
    print '{0:>18}'.format(''), \
          ''.join(['{0:>7}'.format('[{0:}]'.format(str(f)))
                   for f in factor_base])

    clone = log_coefficients[:]
    print '{0:2}'.format(''), \
          '{0:>7}'.format('{} * {}^{} == '.format(b, g, s)), \
          '[{0}]'.format(''.join(['{0:^7}'.format(str(clone[i]))
                                 for i in range(len(clone))])), \
          '(mod {})'.format(p)
    print
    print '{0:2}'.format(''), \
          '{0:>7}'.format('{} * {}^{} == '.format(b, g, s)), \
          ''.join(['{1}^{0} * '.format(i, j)
                   for i, j in zip(log_coefficients, factor_base)]), 1

    print 'L({} * {}^{}) == '.format(b, g, s), \
          ''.join(['{0}*L({1}) + '.format(i, j)
                   for i, j in zip(log_coefficients, factor_base)]), 0, \
          '(mod {})'.format(p - 1)

    print '{0:1}'.format(''), 'L({}) + {} == '.format(b, s), \
          ''.join(['{0}*L({1}) + '.format(i, j)
                   for i, j in zip(log_coefficients, factor_base)]), 0, \
          '(mod {})'.format(p - 1)

    print '{0:7}'.format(''), 'L({}) == '.format(b), \
          ''.join(['{0}*L({1}) + '.format(i, j)
                   for i, j in zip(log_coefficients, factor_base)]), \
          '(-{})'.format(s), '(mod {})'.format(p - 1)

    print '{0:7}'.format(''), 'L({}) == '.format(b), \
          ''.join(['{0}*{1} + '.format(i, j)
                   for i, j in zip(log_coefficients, solution)]), \
          '(-{})'.format(s), '(mod {})'.format(p - 1)

    print '{0:7}'.format(''), 'L({}) == '.format(b), result, \
          '(mod {})'.format(p - 1)

    return

def swap_rows(matrix, row_1, row_2):
    """
    Function to swap the two rows when performing Gauss elimination

    Parameters:
        matrix: a system of linear relations
        row_1, row_2: two rows to be swapped

    Returns:
        a matrix with its two rows swapped
    """
    row = matrix[row_1]
    matrix[row_1] = matrix[row_2]
    matrix[row_2] = row

    return matrix

def copy_matrix(matrix):
    """
    Function to clone a matrix

    Parameters:
        matrix: the original matrix

    Returns:
        a clone of the input matrix
    """
    clone = []
    for row in matrix:
        _row = []
        for value in row:
            _row.append(value)
        clone.append(_row)

    return clone

def finite_field_reduce_matrix(matrix, prime):
    """
    Function to reduce a matrix so that all the pivots become 1

    Parameters:
        matrix: a system of linear relations
        prime: a prime number that is the finite field

    Returns:
        a reduced matrix
    """
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
    """
    Function to find all the primes less than or equal to a certain number

    Parameters:
        B: the given limit to sieve for primes

    Returns:
        a list of primes
    """
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
    Function to find the greatest common divisor for two numbers using
    the Euclidean algorithm

    Parameters:
        a, b: two integers

    Returns:
        the gcd of a and b
    """
    while b != 0:
        temp = b
        b = a % b
        a = temp

    return a

def multiplicative_inverse(n, p):
    """
    Function to find the multiplicative inverse of n in the field p given that
    they are relatively prime

    Parameters:
        n: an element of field p
        p: the field

    Returns:
        an element m of field p such that n*m == 1 (mod p)
    """
    assert(gcd(n, p) == 1)

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
    """
    Function to solve a system of congruences using the Chinese Remainder
    Theorem

    Parameters:
        set_of_modulus: a set of modulus that are relatively prime
        set_of_congruences: a set of congruences

    Returns:
        a solution x that satisifes all the congruences
    """
    product = 1
    sum = 0

    for m in set_of_modulus:
        product *= m

    for i in range(len(set_of_modulus)):
        n = product / set_of_modulus[i]
        inverse = multiplicative_inverse(n, set_of_modulus[i])
        sum += set_of_congruences[i] * inverse * n

    return sum % product

def finite_field_gauss_elimination(matrix, prime):
    """
    Function to solve a system of congruences modulo a prime by performing
    Gauss-Jordan elimination over a finite field

    Parameters:
        matrix: an augmented matrix that consists of a system of linear
                relations and the corresponding solutions
        prime: the prime modulus

    Returns:
        a matrix that has the entries below and above the pivots set to 0
    """
    ncol = len(matrix[0])
    nrow = len(matrix)
    npivot = min(nrow,ncol-1)
    prow = 0
    pcol = prow
    pivot = 0

    while prow < npivot:
        # swap rows to get a nonzero pivot
        if matrix[prow][pcol] == 0:
            for j in range(prow+1,nrow):
                if matrix[j][pcol] != 0:
                    matrix = swap_rows(matrix,prow,j)
                    break

        while matrix[prow][pcol] == 0 and pcol < (ncol - 1):
            pcol += 1
            for j in range(prow+1,nrow):
                if matrix[j][pcol] != 0:
                    matrix = swap_rows(matrix,prow,j)
                    break

        if matrix[prow][pcol] != 0:
            for i in range(nrow):
                if i != prow and matrix[i][pcol] != 0:
                    ratio = (matrix[i][pcol] *
                        multiplicative_inverse(matrix[prow][pcol], prime)
                                               % prime)
                    for j in range(ncol):
                        matrix[i][j] = matrix[i][j] - ratio*matrix[prow][j]
                        matrix[i][j] = matrix[i][j] % prime
        else:
            break

        prow += 1

        if pcol == ncol - 1:
            break

        pcol += 1

    return matrix

def is_smooth(n, factor_base):
    """
    Function to check if a number can be factored completely over a factor
    base

    Parameters:
        n: an integer to check for smoothness
        factor_base: a set of prime factors

    Returns:
        a list of exponents of the factors if a number is smooth
        else, an empty array
    """
    factorization = []
    for p in factor_base:
        exponent = 0
        while n % p == 0:
            n /= p
            exponent += 1
        factorization.append(exponent)

    if n != 1:
        return []

    return factorization

def is_linearly_independent(log_relations, prime):
    """
    Function to check if the rows in a system of linear relations are linearly
    independent so that they yield a unique solution. The function first
    performs Gauss-Jordan elimination on the matrix of logarithmic relations
    for each of the factors of (prime - 1). We then check if all the pivots in
    all the output matrices are nonzero on the diagonal.

    Note: Our logarithmic relations are the congruence classes mod (prime - 1).
    However, since our prime is large, it is certainly odd and so (prime - 1)
    is even. This will cause trouble when performing Gauss elimination because
    some elements have no multiplicative inverses. As a result, we need to
    perform Gauss elimination independently for each of the prime factors of
    (prime - 1). Since (prime - 1) is even, one of the factors is 2. For our
    purpose, we choose a prime such that (prime - 1) has factors of 2 and
    another large prime. We can use quadratic sieve to find the complete
    factorization of our prime.

    Parameters:
        log_relations: a system of linear relations
        prime: the modulus

    Returns:
        True if the pivots are all nonzero on the diagonal, False otherwise
    """
    clone = copy_matrix(log_relations)
    factors = [2, (prime - 1) / 2]

    for f in factors:
        matrix = finite_field_gauss_elimination(clone, f)
        ncol = len(matrix[0])

        # if one of the pivots is zero, the columns are not linearly independent.
        for i in range(ncol - 1):
            if matrix[i][i] == 0:
                return False

    return True

def compute_logarithmic_relations(factor_base, prime, generator):
    """
    Function to search for a set of linearly independent relations between
    the factor base and the powers of the generator. With n primes in the
    factor base, we first collect n + 1 relations to guarantee there could be
    linear independency between the relations. If the relations are linearly
    dependent, we keep collecting more relations until linear independency is
    reached. In other words, there would be a unique solution to these
    relations.

    Args:
        factor_base: a list of primes used for smoothness
        prime: the modulus
        generator: the primitive root of the prime field

    Returns:
        a list of linearly independent smooth relations
    """
    log_relations = []
    limit = len(factor_base) + 1
    next_exponent = 1

    # collects initial smooth relations
    log_relations = find_relations(factor_base, prime, generator,
                                   log_relations, limit, next_exponent)

    # keeps collecting more relations until m linearly independent relations
    # are found
    while not is_linearly_independent(log_relations, prime):
        print '--> LINEARLY DEPENDENT'
        print
        limit += 1
        number_of_relations = len(log_relations)
        number_of_primes = len(factor_base)
        next_exponent = log_relations[number_of_relations - 1][number_of_primes]
        next_exponent += 1
        log_relations = find_relations(factor_base, prime, generator,
                                       log_relations, limit, next_exponent)

    print '--> LINEARLY INDEPENDENT'
    print

    return log_relations

def find_relations(factor_base, prime, generator, log_relations, limit, start):
    """
    Function to search for smooth relations between the factor base and the
    powers of the generator.

    Args:
        factor_base: a list of primes used for smoothness
        prime: the modulus
        generator: the primitive root of the prime field
        log_relations: a list of smooth relations
        limit: a number specifying the number of smooth relations to be found
        start: the next power of the generator to be considered since
               the last found relation

    Returns:
        a list of smooth relations
    """
    for i in range(start, prime):
        relation = []
        sieving_value = (generator ** i) % prime
        relation = is_smooth(sieving_value, factor_base)

        # if the value is smooth
        if not not relation:
            relation.append(i)
            log_relations.append(relation)

        if len(log_relations) == limit:
            break

    print_matrix_with_headers(log_relations, generator, factor_base, prime)
    print

    return log_relations

def solve_for_logarithms(log_relations, prime):
    """
    Function to solve the linear system of linearly independent smooth
    relations using the Gauss-Jordan elimination method and to reduce
    all the pivots to 1.

    Args:
        log_relations: a list of linearly independent smooth relations
        prime: the modulus

    Returns:
        a reduced echelon form of the logarithmic relations
    """
    log_relations = finite_field_gauss_elimination(log_relations, prime)
    log_relations = finite_field_reduce_matrix(log_relations, prime)

    return log_relations

def combine_solutions(solutions, prime):
    """
    Function to find the final solutions from a set of solutions obtained from
    performing Gauss elimination for each of the factors of (prime - 1). We use
    the Chinese Remainder Theorem to combine all the congruences.

    Args:
        solutions: a set of solutions from the factors of (prime - 1)
        prime: the modulus

    Returns:
        the final solution for the linear system of smooth relations over the
        ring (prime - 1)
    """
    factors = [2, (prime - 1) / 2]
    solution = []
    for i in range(len(solutions[0])):
        congruences = []
        for j in range(len(solutions)):
            congruences.append(solutions[j][i])

        s = chinese_remainder(factors, congruences)
        solution.append(s)

    return solution

def find_logarithmic_solution(log_relations, factor_base, prime):
    """
    Function to solve for linear system of smooth relations over the ring
    (prime - 1).

    Args:
        log_relations: a set of linearly independent smooth relations
        factor_base: a list of primes used for smoothness
        prime: the modulus

    Returns:
        the final solution for the linear system of smooth relations
    """
    solutions = []
    factors = [2, (prime - 1) / 2]
    print 'FACTORS OF P - 1:', factors
    print

    for f in factors:
        clone = copy_matrix(log_relations)
        clone = solve_for_logarithms(clone, f)

        _solution = []
        for i in range(len(factor_base)):
            _solution.append(clone[i][len(factor_base)])

        solutions.append(_solution)

        print 'GAUSS ELIMINATION OVER FIELD', f
        print_matrix(clone, f)

        print
        print 'SOLUTION = ', _solution
        print

    print 'SOLUTIONS = ', solutions
    print
    print 'APPLYING THE CHINESE REMAINDER THEOREM YIELDS:'

    combined_solution = combine_solutions(solutions, prime)

    print_log_solution(combined_solution, factor_base)

    return combined_solution

def compute_individual_log(log_solution, prime, generator, output, factor_base):
    """
    Function to solve for the discrete log problem g^x == B (mod p) given
    g as the primitive root of the odd prime field p and B is one of the field
    elements. This function uses the logarithms of the factor base that we
    solved to compute x.

    Args:
        log_solution: the logarithms of the factor base
        prime: the prime field
        generator: the primitive root of the prime field
        output:
        factor_base: a list of primes used for smoothness

    Returns:
        x such that g^x == B (mod p)
    """
    log_coefficients = []
    remainder = 0
    result = 0

    for s in range(1, prime - 2):
        y = (output * generator ** s) % prime
        log_coefficients = is_smooth(y, factor_base)
        if not not log_coefficients:
            remainder = s
            break

    for a, b in zip(log_coefficients, log_solution):
        result += a*b

    result = (result - remainder) % (prime - 1)

    print_solution(factor_base, log_coefficients, log_solution, output,
                   generator, remainder, prime, result)

    return result

def compute_discrete_log(prime, generator, output, B):
    """
    Function to compute the discrete log problem g^x == B (mod p) given
    g as the primitive root of the odd prime field p and B is one of the field
    elements. This function consists of three stages. The first stage consists
    of searching for a set of linearly independent smooth relations between the
    factor base and the powers of the generator of the prime field. The second
    stage finds the logarithms of the factor base using the system of smooth
    relations found in the first stage. The third stage computes the discrete
    log of the desired field element.

    Args:
        prime: the prime field
        generator: the primitive root of the prime field
        output: the desired field element of which to find the discrete log
        B: the smoothness bound to find the factor base

    Returns:
        the solution to the discrete log problem
    """
    factor_base = sieve_of_Eratosthenes(B)
    print 'FACTOR BASE: ', factor_base
    print

    # Stage 1
    print '{:*<95}'.format('')
    print '{0:^80}'.format('STEP 1: FIND LOGARITHMIC SMOOTH RELATIONS')
    print
    log_relations = compute_logarithmic_relations(factor_base, prime, generator)

    print 'SMOOTH LOGARITHMIC RELATIONS'
    print_matrix_with_headers(log_relations, generator, factor_base, prime)
    print
    print '{:*<95}'.format('')
    print '{0:^80}'.format('STEP 2: SOLVE FOR THE LOGARITHMS')
    print
    print 'AUGMENTED MATRIX FOR GAUSS ELIMINATION OVER RING', (prime -  1)
    print_augmented_matrix(log_relations, generator, factor_base, prime)
    print

    # Stage 2
    solution = find_logarithmic_solution(log_relations, factor_base, prime)

    # Stage 3
    print '{:*<95}'.format('')
    print '{0:^80}'.format('STEP 3: COMPUTE THE INDIVIDUAL LOGARITHM')
    print
    result = compute_individual_log(solution, prime,
                                    generator, output, factor_base)

    return result

def main():
    """
    Function to specify the discrete log problem and solve it.

    Args:
        none

    Returns:
        none
    """
    generator = 5
    prime = 14087
    output = 5872
    B = 15

    print 'GENERATOR: ', generator
    print 'FIELD: ', prime
    print 'DISCRETE LOG PROBLEM: {}^X == {} (mod {})'.format(generator, output,
                                                             prime)
    print

    result = compute_discrete_log(prime, generator, output, B)

    print '\nX SUCH THAT {}^X == {} (mod {}) is'.format(generator,
           output, prime), result
    print

    return

if __name__ == '__main__':
    main()
