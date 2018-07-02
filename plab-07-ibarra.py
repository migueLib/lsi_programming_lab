import numpy as np


def get_scores(i, j, seq, loop=1):
    if j-i+1 < loop+2:
        return 0
    else:
        # Case 1,2
        unpaired = get_scores(i+1, j-1, seq, loop=loop) + delta.get((seq[i], seq[j]), 0)

        # Case 3,4
        paired = [get_scores(i, k, seq, loop=loop) +
                  get_scores(k + 1, j, seq, loop=loop)
                  for k in range(i, j) if delta.get((seq[k], seq[j]), 0)]

        if not paired:
            paired = [0]

        paired = max(paired)

    return max(unpaired, paired)


def matched_base_pairs(i, j, seq, matrix, loop=1):
    matched_positions = []
    # Case 1: the sequence length must be greater than 2 to continue
    if j - i + 1 >= loop + 2:
        if matrix[i][j] == matrix[i + 1][j - 1] + delta.get((seq[i], seq[j]), 0):
            # Case 2 where if new match made then append to list
            if delta.get((seq[i], seq[j]), 0):
                matched_positions.append((i, j))
            if matched_base_pairs(i + 1, j - 1, seq, matrix, loop=loop):
                matched_positions.extend(matched_base_pairs(i + 1, j - 1, seq, matrix,
                                                            loop=loop))
        else:
            # Case 3, find k in range (i,j) and get list of paired of comp bases then join them
            for k in range(i, j - 1):
                if matrix[i][j] == matrix[i][k] + matrix[k + 1][j]:
                    first_half = matched_base_pairs(i, k, seq, matrix, loop=loop)
                    second_half = matched_base_pairs(k + 1, j, seq, matrix, loop=loop)
                    if first_half and second_half:
                        matched_positions.extend((first_half, second_half))
                # Once we have found that k value we break from the loop
                break

    return matched_positions


def backtrack(seq, matrix, loop=1):
    positions = list()
    for i in range(len(seq)):
        for j in range(i, len(seq)):
            match = matched_base_pairs(i, j, seq, matrix, loop=loop)
            if len(match) > 0:
                positions.append(match)
    return positions


def pretty_plot(seq, optimal):
    # Create string full of dots
    out_string = ["." for i in range(len(seq))]

    # Place the opening and closing parenthesis
    while optimal:
        opened, closed = optimal.pop()
        out_string[opened] = "("
        out_string[closed] = ")"

    return ''.join(out_string)


def score_nussinov(seq, loop=1):
    matrix = np.zeros((len(seq), len(seq)), dtype='int')

    # Generate scores
    for i in range(len(seq)):
        for j in range(i, len(seq)):
            # OK, here is the main part, here es where we get the scores.
            matrix[i][j] = get_scores(i, j, seq, loop=loop)

    return matrix

def scoremethis(seq):
    # Create matrix
    matrix = np.zeros((len(seq),len(seq)), dtype="int")

    # Iterate over the matrix diagonals starting from (0,1)
    for d in range(len(seq)):
        for i, j in zip(range(len(seq)-1), range(d, len(seq))):
            # Case 1 (down)
            c1 = matrix[i+1, j]

            # Case 2 (left)
            c2 = matrix[i, j-1]

            # Case 3 (diagonal +  delta(match))
            c3 = matrix[i+1, j-1] + delta.get((seq[i], seq[j]), 0)

            # Case 4 (the weird one)
            c4 = [matrix[i, k] + matrix[k+1, j] for k in range(i+1, j)]
            c4 = max(c4) if c4 else 0

            # Add the max to the matrix
            matrix[i, j] = max([c1, c2, c3, c4])

    return matrix

#######
# MAIN
#######
# Set a global variable for delta
delta = {('A', 'U'): 1, ('C', 'G'): 1, ('G', 'C'): 1, ('U', 'A'): 1}

#seq = "AUCGGAGCAUUUUUUGCUCCGACGCAGCCUCAUGCUUUUUU"
seq = "AUCGCAGCAUUUUUUGCUCCGA"
scores = scoremethis(seq)
print(scores)
score_matrix = score_nussinov(seq, loop=1)

# Find number of optimal solutions
solutions = backtrack(seq, score_matrix)

# Display first optimal matching
optimal = max(solutions, key=len)

# Get optimal solutions
n_optimal = sum([1 for s in solutions if len(s) == len(optimal)])

print(seq)
print(pretty_plot(seq, optimal))
