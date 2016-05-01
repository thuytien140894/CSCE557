"""
Program to factor a number using the quadratic sieve algorithm.

Given that N is the number we would like to factor, R is the square root of N,
this program first uses the function (R + x)^2 - N to sieve through values of x
in a chosen interval and find smooth numbers that factor complete over a factor
base. We build a matrix whose rows contain the exponents of the primes in the
factor base for each of the smooth numbers. Then we find a combination of rows
that would yield a perfect square by performing the Gauss-Jordan elimination.
We use these rows to compute x and y in the relation x^2 == y^2 (mod N).
Taking the gcd(x - y, N) or the gcd(x + y, N) would give a factor of N.

Note: If we cannot find smooth numbers, we need to readjust the sieving
intervals as well as the B bound to find the factor base. Sieving is successful
if a large number of smooth numbers are found and that they are not all
linearly independent. If they are linearly independent, we cannot find any
free variables and as a result cannot find a combination of smooth relations
that yield a perfect square. We check for these two exception cases in our code.
These cases can be avoided using the multiple polynomial quadratic sieve
algorithm in which more than one function is used to perform the sieve.

Authors: Steven Dao, Tien Ho, and Brandon Hostetter
Date: 28 April 2016
"""

import math, random, sys

# Constants
# NUM_TO_FAC = 1692966159481
# NUM_TO_FAC = 5298130281889
# NUM_TO_FAC = 39726901372343
NUM_TO_FAC = 239812014798221
PRIMES_BELOW = 1000
SIEVING_LENGTH = 100
NUM_TO_SET = 1
HALF_SIEVING_LENGTH = int(SIEVING_LENGTH / 2)
SQUARE_ROOT = int(math.sqrt(NUM_TO_FAC))

def sieve_of_eratosthenes():
    """
    Find all primes less than or equal to 'PRIMES_BELOW' using the sieve of
    eratosthenes algorithm.

    Returns:
        list_of_primes: A list of all primes that exist from 0 to 'max_val'.
    """

    max_val = PRIMES_BELOW + 1
    int_range = range(2, int(math.ceil(math.sqrt(max_val))))
    prime_index_list = [True] * (max_val - 2)
    list_of_primes = []

    for i in int_range:
        if prime_index_list[i - 2] == True:
            if i * i > max_val:
                # We won't eliminate anymore values, break
                break
            for j in range(i * i, max_val, i) :
                # Start at i^2, remove every value that is a multiple of i
                prime_index_list[j - 2] = False

    for i in range(len(prime_index_list)):
        if prime_index_list[i] == True:
            list_of_primes.append(i + 2)

    return list_of_primes

def is_quadratic_residue(p):
    """
    Check to see if 'NUM_TO_FAC' is a quadratic residue modulo p. That is,
    if there exists an integer x such that x^2 is congruent to 'NUM_TO_FAC'
    (mod p), then 'NUM_TO_FAC' is a quadratic residue modulo p.

    Args:
        p: A prime number.

    Returns:
        bool: True if 'NUM_TO_FAC' is a quadratic residue modulo p. False if
            'NUM_TO_FAC' is a quadratic non-residue modulo p.
    """

    if (p != 2):
        legendre_symbol = NUM_TO_FAC ** ((p - 1) / 2) % p
        if legendre_symbol == p - 1:
            return False
    return True

def find_factor_base():
    """
    Find all primes less than 'PRIMES_BELOW' and keep a list of all primes
    for which 'NUM_TO_FAC' is a quadratic residue. This is used to eliminate
    primes that are not useful when sieving.

    Returns:
        factor_base: The list of all primes returned by our
            'sieve_of_eratosthenes' method that satisify the condition that
            'NUM_TO_FAC' is a quadratic residue mod p.
    """

    primes = sieve_of_eratosthenes()
    factor_base = []

    for p in primes:
        if is_quadratic_residue(p):
            factor_base.append(p)

    return factor_base

def compute_sieving_values():
    """
    Compute the values to be sieved over by running a loop over
    'SIEVING_LENGTH' using the formula (r-x)^2 - f where r is 'SQUARE_ROOT',
    x in the index in our loop, and f is 'NUM_TO_FAC'. We loop from 0 to
    'SIEVING_LENGTH', however, we want the values of x to range from
    [-'SIEVING_LENGTH' / 2, 'SIEVING_LENGTH' / 2]. To get the correct values
    for x, we subtract a constant 'HALF_SIEVING_LENGTH'.

    Returns:
        sieving_values: A list of integers to be used for sieving.
    """

    print 'ROOT IS	', SQUARE_ROOT
    sieving_values = []

    for i in range(SIEVING_LENGTH):
        x = i - HALF_SIEVING_LENGTH
        f_x = (SQUARE_ROOT + x) ** 2 - NUM_TO_FAC
        sieving_values.append(f_x)
        print '{0:8d} {1:8d} {2:8d}'.format(x, i, f_x)

    return sieving_values

def tonelli_shanks(n, p):
    """
    Solve for a congruence X^2 == N (mod p) where N is a quadratic residue
    mof p and p is a mod prime. We use this method to find the starting points
    in our sieving intervals whose value of the form (R + x)^2 would
    definitely divide by the prime in consideration.

    Args:
        n: An integer that is a quadratic residue of p
        p: An odd prime

    Returns:
        r: one of the two solutions for X^2 == N (mod p)
    """

    # if p is an even prime, then the congruence X^2 == N (mod p) has one
    # solution of 1
    if p % 2 == 0:
        return 1

    # check to make sure that n is a quadratic residue of a mod prime p
    # using the Legendre symbol
    legendre_symbol = n ** ((p - 1) / 2) % p

    if legendre_symbol != 1:
        return 0

    # write p-1 such that p-1 = q2^s with q being odd
    s = 1
    q = (p - 1) / (2 ** s)
    while q % 2 == 0:
        s += 1
        q = (p - 1) / (2 ** s)

    # pick a quadratic non-residue z starting from 1 to p-1
    z = 0
    for i in range(1, p):
        legendre_symbol = i ** ((p - 1) / 2) % p
        if legendre_symbol == (p - 1):
            z = i
            break

    # initialize variables
    r = n ** ((q + 1) / 2) % p
    t = n ** q
    m = s
    c = z ** q % p

    while t % p != 1:
        i = 0
        # find the lowest i in (0, m) such that t^(2^i) == 1 (mod p)
        for j in range(1, m):
            k = 2 ** j
            if t ** k % p == 1:
                i = j
                break

        b = c ** (2 ** (m - i - 1)) % p
        r = (r * b) % p
        t = (t * (b ** 2)) % p
        c = b ** 2 % p
        m = i

    r = r % p

    return r

def sieve(factor_base, sieving_values):
    """
    Find numbers in 'sieving_values' that are smooth. If a number
    can be divided entirely over the factor base, it is
    considered smooth. At the end, the index of each 1 in 'smooth_numbers'
    will coorespond to the index of the number in 'sieving_values' which is
    smooth.

    Args:
        factor_base: Our factor base of primes numbers.
        sieving_values: A list of values to be checked for smoothness.

    Returns:
        smooth_numbers: A list of the same length as 'sieving_values'
                        where the index of each 1 in this list will
                        coorespond to the index of a number in
                        'sieving_values' that is smooth.
    """
    smooth_numbers = sieving_values[:]
    lower_limit = -HALF_SIEVING_LENGTH

    for prime in factor_base:
        print '{0:1} {1:4d}'.format('SIEVE WITH', prime)
        new_array = []
        temp = tonelli_shanks(NUM_TO_FAC, prime)
        s1 = temp - SQUARE_ROOT
        s2 = (prime - temp) - SQUARE_ROOT
        offset1 = math.ceil(float(lower_limit - s1) / prime)
        offset2 = math.ceil(float(lower_limit - s2) / prime)
        x1 = s1 + prime * offset1
        x2 = s2 + prime * offset2

        x_values = range(lower_limit, -lower_limit)
        if x1 in x_values:
            x1_index = x_values.index(x1)
            for i in range(x1_index, SIEVING_LENGTH, prime):
                while smooth_numbers[i] % prime == 0:
                    temp = smooth_numbers[i]
                    smooth_numbers[i] /= prime
                    print '{0:5} {1:8d} {2:8d} {3:8d}'.format('DIVIDE', prime,
                           abs(temp), abs(smooth_numbers[i]))

        if x2 != x1 and x2 in x_values:
            x2_index = x_values.index(x2)
            for j in range(x2_index, SIEVING_LENGTH, prime):
                while smooth_numbers[j] % prime == 0:
                    temp = smooth_numbers[j]
                    smooth_numbers[j] /= prime
                    print '{0:5} {1:8d} {2:8d} {3:8d}'.format('DIVIDE', prime,
                           abs(temp), abs(smooth_numbers[j]))

        print_matrix(new_array)

    return smooth_numbers

def refactor_smooth_values(factor_base, sieving_values, smooth_numbers):
    """
    Refactor the numbers in 'sieving_values' for the numbers that are
    smooth. For each smooth number in 'sieving_values', build a list
    that holds the exponents that of the primes in the 'factor_base'
    that factor the given number. Not all numbers in 'sieving_values'
    are smooth, so we track the row indices for the smooth numbers in
    'sieving_values' so we can easily find them again without recomputing
    what the smooth numbers are.

    Args:
        factor_base: Our factor base of prime numbers.
        sieving_values: A list of values that will be refactored.
        smooth_numbers: A list of the same length as 'sieving_values'
                        where a value of 1 means the number at that index
                        in 'sieving_values' is smooth.

    Returns:
        exponent_matrix: The matrix of the exponents used to factor
                         each of the smooth numbers in 'sieving_values'.
        row_indices: A list containing the indices of the smooth numbers.
    """
    print '\nROOT IS  ', SQUARE_ROOT
    exponent_matrix = []
    smooth_numbers = map(abs, smooth_numbers)
    sieving_values = sieving_values[:]
    row_indices = []

    print 'FB   ', factor_base, '\n'

    for i in range(len(smooth_numbers)):
        if smooth_numbers[i] == 1:
            exponent_matrix.append([])
            exponent_matrix_last_row = len(exponent_matrix) - 1
            print '{0:5} {1:8d} {2:8d} {3:8d}'.format('SMOOTH', i - HALF_SIEVING_LENGTH, i, sieving_values[i])
            for j in range(len(factor_base)):
                exponent_matrix[exponent_matrix_last_row].append(0)
                if factor_base[j] == -1:
                    if sieving_values[i] < 0:
                        sieving_values[i] = abs(sieving_values[i])
                        exponent_matrix[exponent_matrix_last_row][j] += 1
                        continue
                    continue
                while sieving_values[i] % factor_base[j] == 0:
                    print '{0:5} {1:8d} {2:8d} {3:8d}'.format('DIVIDE', j, factor_base[j], sieving_values[i])
                    exponent_matrix[exponent_matrix_last_row][j] += 1
                    sieving_values[i] /= factor_base[j]
            print '{0:5} {1:8d} {2:8d} {3:8d}'.format('SMOOTH', i - HALF_SIEVING_LENGTH, i, sieving_values[i])
            print exponent_matrix[exponent_matrix_last_row], '\n'
            row_indices.append(i - HALF_SIEVING_LENGTH)
        else:
            print '{0:6} {1:8d} {2:8d} {3:8d} {4}'.format('NOT', i - HALF_SIEVING_LENGTH, i, sieving_values[i], '\n')

    # calculations are done, just printing output now
    iterator = iter(row_indices)
    counter = 0
    for i in range(len(exponent_matrix)):
        for col in exponent_matrix[i]:
            if col != 0:
                print '{0:5d} {1:5d} {2}'.format(counter, next(iterator), exponent_matrix[i])
                counter += 1
                break

    return [exponent_matrix, row_indices]

def get_binary_matrix(matrix):
    """
    Mod 2 every value in a matrix.

    Args:
        matrix: A 2D list with numeric values.

    Returns:
        binary_matrix: The 'matrix' with all values mod 2.
    """
    binary_matrix = [map(mod_2, x) for x in matrix]
    return binary_matrix

def mod_2(num):
    """
    Mod a value by 2.

    Args:
        num: Number to be modded by 2.

    Returns:
        Mod 2 of 'num'.
    """
    return num % 2

def swap_rows(matrix, row1, row2):
    """
    Function to take two rows in a matrix and swap their positions.

    Args:
        matrix: The entire matrix containing the two rows.
        row1, row2: The two rows to be swapped.

    Returns:
        matrix: The matrix with row1 and row2 swapped.
    """
    row = matrix[row1]
    matrix[row1] = matrix[row2]
    matrix[row2] = row
    return matrix

def multiplicative_inverse(n, p):
    """
    Compute the multiplicative inverse of a field element n in the finite
    field Z_p.

    Args:
        n: an element of field p
        p: the field

    Returns:
        an element x1 of field Z_p such that n*x1 == 1 (mod p)
    """
    # n only has a multiplicative inverse if n and p are relatively prime.
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

def gcd(a, b):
    """
    Compute the greatest common divisor of two integers using the Euclidean
    algorithm.

    Args:
        a, b: two integers

    Returns:
        the gcd of a and b
    """
    while b != 0:
        temp = b
        b = a % b
        a = temp
    return a

def finite_field_gauss_elimination(matrix, prime):
    """
    Solve a system of congruences modulo a prime by performing
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
    npivot = min(nrow, ncol)
    prow = 0
    pcol = prow
    pivot = 0

    while prow < npivot:
        # swap rows to get a nonzero pivot
        if matrix[prow][pcol] == 0:
            for j in range(prow + 1, nrow):
                if matrix[j][pcol] != 0:
                    matrix = swap_rows(matrix, prow, j)
                    break

        # if we cannot find a nonzero for our pivot in the current column,
        # this column is free and we move to the next
        while matrix[prow][pcol] == 0 and pcol < (ncol - 1):
            pcol += 1
            for j in range(prow + 1, nrow):
                if matrix[j][pcol] != 0:
                    matrix = swap_rows(matrix ,prow, j)
                    break

        if matrix[prow][pcol] != 0:
            for i in range(nrow):
                if i != prow and matrix[i][pcol] != 0:
                    ratio = (matrix[i][pcol] *
                        multiplicative_inverse(matrix[prow][pcol], prime)
                                               % prime)
                    for j in range(ncol):
                        matrix[i][j] = matrix[i][j] - ratio * matrix[prow][j]
                        matrix[i][j] = matrix[i][j] % prime
        else:
            break

        prow += 1

        if pcol == ncol - 1:
            break

        pcol += 1

    return matrix

def print_matrix(matrix):
    """
    Function to print a formatted matrix to the console.

    Args:
    	matrix: The matrix to be formatted and printed.
    """
    for row in matrix:
        print '[{0}]'.format(''.join(['{0:^3}'.format(str(row[i]))
                                     for i in range(len(row))]))
    return

def find_pivot_columns(solution_matrix):
    """
    Find basic variables that correspond to pivot columns in the reduced
    echelon matrix.

    Args:
        solution_matrix: a reduced echelon matrix resulted from the Gauss-Jordan
                         elimination

    Returns:
        a dictionary that contains the basic variables as the keys and their
        corresponding rows as the values. The basic variables are indicated by the
        column indices in the matrix
    """
    pcol = 0
    pivot_variables = {}

    for prow in range(len(solution_matrix)):
        while solution_matrix[prow][pcol] == 0:
            if pcol == len(solution_matrix[0]) - 1:
                break

            pcol += 1

        if pcol == len(solution_matrix[0]) - 1:
            if solution_matrix[prow][pcol] != 0:
                pivot_variables[pcol + 1] = prow
            break

        pivot_variables[pcol + 1] = prow

        pcol += 1

    return pivot_variables

def find_free_variables(pivot_variables, number_of_variables):
    """
    Find free columns in a reduced echelon matrix by counting the columns
    that do not correspond to pivot variables.

    Args:
    	pivot_variables: basic variables
        number_of_variables: the number of columns in the reduced echelon
                             matrix

    Returns:
    	a list of columns indices that correspond to free variables
    """
    variable_list = range(1, number_of_variables + 1)
    free_variables = [x for x in variable_list
                      if x not in pivot_variables]

    return free_variables

def solve_for_nullspace(solution_matrix, pivot_variables, free_variables):
    """
  	Find ble combination of rows in the original matrix of smooth
    relations that yield a perfect square. This function randomly picks one of
    the free variables and set it to 1 while setting the other free variables
    to 0. Then we compute the values of the pivot variables accordingly.

    Args:
        solution_matrix: a reduced echelon matrix
        pivot_variables: basic variables
        free_variables: free variables

    Returns:
        a solution to the nullspace of the solution_matrix whose indices of
        1s represent the row numbers in the original matric of smooth relations
        that would yield a perfect square
    """
    chosen_variables = random.sample(free_variables, NUM_TO_SET)
    print
    print 'VALUES FOR FREE VARIABLES: '
    print ''.join(['X_{} == 1\n'.format(i) for i in chosen_variables])
    print ''.join(['X_{} == 0\n'.format(i) for i in free_variables
                   if i not in chosen_variables])

    solution = [0 for i in range(len(solution_matrix[0]))]

    # set the chosen free variables to 1
    for n in chosen_variables:
        solution[n - 1] = 1

    # solve for pivot variables using the chosen free variables
    for variable, row in pivot_variables.items():
        for b in chosen_variables:
            solution[variable - 1] += solution_matrix[row][b - 1]
            solution[variable - 1] = solution[variable - 1] % 2

    print 'SOLUTIONS FOR PIVOT VARIABLES USING FREE VARIABLES: '
    print ''.join(['X_{} == {}\n'.format(i + 1, solution[i])
                   for i in range(len(solution)) \
                   if (i + 1) not in free_variables])

    return solution

def find_perfect_square_rows(exponent_matrix, solution_matrix, pivot_variables,
                             free_variables):
    """
    Find a combination of rows of smooth relations that result in a perfect
    square when multiplying all the powers of the factor base. This method
    first finds the pivot variables and the free variables from the input
    echelon form. Then we call the function solve_for_nullspace to find one
    possible solution for the nullspace of the matrix. Finally, we use the
    solution to determine the row numbers in the original matrix of relations
    that would be used for the congruence X^2 == Y^2 (mod p).

    Args:
        solution_matrix: An echelon form of a matrix resulted from the Gauss
                         elimination.
        pivot_variables: A list of the pivot variables in the solution matrix.
        free_variables: A list of the free variables in the solution matrix.

    Returns:
        A list of row numbers in the original matrix of smooth relations.
    """
    solution = solve_for_nullspace(solution_matrix, pivot_variables,
                                   free_variables)
    print 'ONE POSSIBLE SOLUTION:', solution

    subset = []
    for j in range(len(solution)):
        if solution[j] == 1:
            subset.append(j)

    print '\nROWS THAT YIELD A PERFECT SQUARE:'
    for i in subset:
        print exponent_matrix[i]

    return subset

def transpose(matrix):
    """
    Transpose an input matrix by making the rows the columns of the new matrix.

    Args:
        matrix: An array of rows.

    Returns:
        transposed_matrix: A transposed matrix.
    """
    transposed_matrix = [list(x) for x in zip(*matrix)]

    return transposed_matrix

def solve_for_smooth_relations(binary_matrix):
    """
    Function to solve the linear system of smooth relations using the
    Gauss-Jordan elimination method.

    Args:
    	binary_matrix: a matrix of smooth relations mod 2

    Returns:
    	a reduced echelon matrix
    """
    print '\nBINARY MATRIX'
    print_matrix(binary_matrix)
    print

    # transpose the matrix
    transposed_matrix = transpose(binary_matrix)
    print 'TRANSPOSED MATRIX'
    print ''.join(['{0:3}'.format(i) for i in range(1, len(binary_matrix) + 1)])
    print
    print_matrix(transposed_matrix)

    # solve for the binary matrix using Gauss elimination
    solution_matrix = finite_field_gauss_elimination(transposed_matrix, 2)

    print '\nSOLUTION:'
    print ''.join(['{0:3}'.format(i) for i in range(1, len(binary_matrix) + 1)])
    print
    print_matrix(solution_matrix)

    return solution_matrix

def find_lhs(row_indices, perfect_square_rows):
    """
    Finds the left-hand side by computing (r + x)^ 2 for each row that forms a
    perfect square and taking its square root. 'r' is the square root of the
    factored number while 'x' is the sieving interval value for each perfect
    square row.

    Args:
    	row_indices: a list of the indices of smooth numbers
        perfect_square_rows: a list of row numbers in the original matrix of
                             smooth relations

    Returns:
    	a whole number for the left-hand side
    """
    lhs_value = 1
    for r in perfect_square_rows:
        lhs_value *= ((SQUARE_ROOT + row_indices[r]) ** 2)

    return math.sqrt(lhs_value)

def find_rhs(row_indices, perfect_square_rows, sieving_values):
    """
    Finds the right-hand side by continuously multiplying previous values by
    the sieving values with the sum of row_indices and half of the sieving
    length as the index.

    Args:
    	row_indices: a list of the indices of smooth numbers
        perfect_square_rows: a list of row numbers in the original matrix of
        					 smooth relations
        sieving_values: a list of values used to sieve

    Returns:
		a whole number for the right-hand side
    """
    rhs = 1
    for i in perfect_square_rows:
        rhs *= sieving_values[row_indices[i] + HALF_SIEVING_LENGTH]

    return math.sqrt(rhs)

def find_factors(exponent_matrix, solution_matrix, pivot_variables,
                 free_variables, row_indices, sieving_values):
    """
    Find the factors for 'NUM_TO_FAC'. Calculate the x and y values
    for the relation x^2 == y^2 (mod N) with helper functions. If x and y
    are not equal and x and -y are not equal, then we have found the factors.
    If these conditions are not met, find a new set of perfect square rows
    and recalculate x and y. The two factors of 'NUM_TO_FAC'
    are gcd(x +/- y, NUM_TO_FAC).

    Args:
        solution_matrix: a reduced echelon matrix
        pivot_variables: basic variables
        free_variables: free variables
        row_indices: The indices of the smooth numbers in 'sieving_values'.
        sieving_values: The original list of sieving values to be used to
                        find the x and y values.

    Returns:
        The two factors of 'NUM_TO_FAC'.
    """
    perfect_square_rows = find_perfect_square_rows(exponent_matrix,
                                                   solution_matrix,
                                                   pivot_variables,
                                                   free_variables)

    x = find_lhs(row_indices, perfect_square_rows)
    y = find_rhs(row_indices, perfect_square_rows, sieving_values)

    # make sure that x !== y (mod N) or x! == -y (mod N), or else we won't
    # obtain a factor
    while x % NUM_TO_FAC == y % NUM_TO_FAC:
        perfect_square_rows = find_perfect_square_rows(exponent_matrix,
                                                       solution_matrix,
                                                       pivot_variables,
                                                       free_variables)
        x = find_lhs(row_indices, perfect_square_rows)
        y = find_rhs(row_indices, perfect_square_rows, sieving_values)

    print
    print 'x =', x
    print 'y =', y
    print
    print 'FACTORS'
    first_factor = gcd(x - y, NUM_TO_FAC)
    print first_factor
    second_factor = NUM_TO_FAC / first_factor
    print second_factor

    return

def main():
    """
    The main method to perform quadratic sieve.
    """
    factor_base = find_factor_base()
    print 'FB   ', factor_base

    sieving_values = compute_sieving_values()
    smooth_numbers = sieve(factor_base, sieving_values)

    factor_base.insert(0,-1)
    # refactor the smooth numbers
    refactored_smooth_values = refactor_smooth_values(factor_base,
                                                      sieving_values,
                                                      smooth_numbers)
    # get the matrix of smooth relations
    exponent_matrix = refactored_smooth_values[0]

    # check if there are smooth relations.
    if len(exponent_matrix) == 0:
        print 'There are no smooth numbers. Try another interval or factor base.'
        sys.exit(0)


    # find one possible solution for the system of relations
    binary_matrix = get_binary_matrix(exponent_matrix)
    solution_matrix = solve_for_smooth_relations(binary_matrix)
    pivot_variables = find_pivot_columns(solution_matrix)
    free_variables = find_free_variables(pivot_variables.keys(),
                                         len(solution_matrix[0]))

    print '\nPIVOT VARIABLES:', ''.join(['X_{}  '.format(i)
                                         for i in pivot_variables.keys()])
    print '\nFREE VARIABLES:', ''.join(['X_{}  '.format(i)
                                        for i in free_variables])

    # check if there are free variables. If not, the only solution we can get
    # from the solution_matrix is all 0s, which means we cannot find any rows
    if len(free_variables) == 0:
        print 'There are no free variables. Cannot find perfect square rows.'
        sys.exit(0)

    # get the x values that give smooth values
    row_indices = refactored_smooth_values[1]
    find_factors(exponent_matrix, solution_matrix, pivot_variables,
                 free_variables, row_indices, sieving_values)

    return

# Call the main method to begin program
if __name__ == '__main__':
    main()
