import os
import argparse
import pandas as pd
import numpy as np
from argparse import RawDescriptionHelpFormatter


def get_options():
    """
    Function to pull in arguments
    """

    description = """ Implements the Needle-Wunch algorith for sequence allingment """

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawDescriptionHelpFormatter)
    # Standard Input
    standard = parser.add_argument_group(title='Standard input',
                                         description='Standard input for SECIM tools.')
    standard.add_argument('-sm', "--score", dest="score", action='store',
                          required=True, help="Path for input alignment")
    standard.add_argument('-g', "--gap", dest="gap", action='store', type=int,
                          required=False, help="Gap penalty")
    standard.add_argument('-e', "--gap_extension", dest="gap_extension", action='store', type=int,
                          required=False, help="Gap penalty")
    standard.add_argument('-s1', "--seq1", dest="seq1", action='store',
                          required=False, help="Path for input 1 sequence in fasta format")
    standard.add_argument('-s2', "--seq2", dest="seq2", action='store',
                          required=False, help="Path for input 2 sequence in fasta format")
    standard.add_argument('-o', "--output", dest="output", action='store',
                          required=False, help="Path to output matrix in wide format.")

    args = parser.parse_args()

    # Standardize paths
    args.score = os.path.abspath(args.score)
    args.seq1 = os.path.abspath(args.seq1)
    args.seq2 = os.path.abspath(args.seq2)
    args.output = os.path.abspath(args.output)

    return args


def read_blosum(path):
    """
    Reads int BLOSUM file

    :param path: pathway of the blosum matrix
    :return: DataFrame with the imported BLOSUM matrix
    """
    matrix = list()
    is_head = True
    with open(path, "r") as IN:
        for line in IN:
            line = line.strip()
            line = line.replace("  ", " ")
            line = line.replace(" \t", "\t")
            line = line.replace("\t ", "\t")
            line = line.replace(" ", "\t")
            line = line.split("\t")
            if is_head:
                head = ["row"]+line
                is_head = False
            else:
                matrix.append(line)

    # Creates data frame
    matrix = pd.DataFrame(matrix, columns=head)
    matrix = matrix.set_index(["row"])

    return matrix


def read_seq(path):
    """
    Reads sequence from fasta file

    :param path: path of input fasta file
    :return: string with the sequence contain in the fasta file
    """
    sequence = []
    with open(path, "r") as IN:
        for line in IN:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence)


def get_global_alignment(seq1, seq2, w, substitution):
    """
    Get alignment score for global alignment

    :param seq1: sequence 1
    :param seq2: sequence 2
    :param w: gap penalty
    :param substitution:  substitution matrix
    :return: np.matrix filled with
    """
    # Create empty matrix
    m = np.zeros((len(seq1)+1, len(seq2)+1))

    # Fill gaps i=0 and j=0; this one is for rows
    for i in range(1, len(seq1)+1):
        m[i, 0] = m[i - 1, 0] + w

    # This one is for columns
    for i in range(1, len(seq2)+1):
        m[0, i] = m[0, i-1] + w

    # Fill the rest of the matrix
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            m[i, j] = max([m[i, j-1] + w, m[i-1, j] + w, m[i-1, j-1] + int(substitution[seq1[i-1]][seq2[j-1]])])

    # Creating a data frame with the scores
    # print(m.shape)
    # m = pd.DataFrame(data=m, index=[""]+list(seq1), columns=[""]+list(seq2))

    return m


def get_local_alignment(seq1, seq2, w, substitution):
    """
    Get alignment score for global alignment

    :param seq1: sequence 1
    :param seq2: sequence 2
    :param w: gap penalty
    :param substitution:  substitution matrix
    :return: np.matrix filled with
    """
    # Create empty matrix
    m = np.zeros((len(seq1)+1, len(seq2)+1))

    # Fill gaps i=0 and j=0; this one is for rows
    for i in range(1, len(seq1)+1):
        m[i, 0] = m[i - 1, 0] + w

    # This one is for columns
    for i in range(1, len(seq2)+1):
        m[0, i] = m[0, i-1] + w

    # Fill the rest of the matrix
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            m[i, j] = max([0, m[i, j-1] + w, m[i-1, j] + w, m[i-1, j-1] + int(substitution[seq1[i-1]][seq2[j-1]])])

    # Creating a data frame with the scores
    # print(m.shape)
    # m = pd.DataFrame(data=m, index=[""]+list(seq1), columns=[""]+list(seq2))

    return m


def get_traceback(m, seq1, seq2):
    """
    Trace back the best global alignment

    :param m: scoring matrix
    :param seq1: secuencia 1
    :param seq2: secuencia 2
    :return: list with aligned pairs
    """
    # Initial coordinates
    i, j = len(seq1), len(seq2)
    aligned = list()

    # YA ME QUIERO MORIR, ESTA TAREA ME VA A MATAR
    while i > 0 and j > 0:
        a = m[i-1, j]
        b = m[i-1, j-1]
        c = m[i, j-1]

        # Get matches/mismatches/gaps and re assign coordinates
        if b >= a and b >= c:
            i, j = i-1, j-1
            aligned.append([seq1[i], seq2[j]])
        elif a >= c:
            i, j = i-1, j
            aligned.append([seq1[i], "-"])
        else:
            i, j = i, j-1
            aligned.append(["-", seq2[j]])

    return aligned


def formatted_alignment(aligned, substitution):
    symbol = list()
    for a, b in aligned:

        if a == b:
            symbol.append([a, "|", b])

        elif a == "-" or b == "-":
            symbol.append([a, " ", b])

        elif int(substitution[a][b]) > 0:
            symbol.append([a, ":", b])

        else:
            symbol.append([a, " ", b])

    symbol.reverse()
            
    return symbol


def main(args):
    # Read blosum matrix
    blosum = read_blosum(args.score)

    # Read fasta sequence
    seq1 = read_seq(args.seq1)
    seq2 = read_seq(args.seq2)

    # Get global and local alignment scores
    global_scores = get_global_alignment(seq1, seq2, args.gap, blosum)
    local_scores = get_local_alignment(seq1, seq2, args.gap, blosum)

    # Print alignment matrix
    max_col = np.argmax(local_scores, axis=0)
    print(np.argmax(max_col))

    # print([[i, np.argmax(local_scores, axis=i)] for i in range(len(local_scores))])
    local_scores = pd.DataFrame(local_scores, index=[""]+list(seq1), columns=[""]+list(seq2))
    local_scores.to_csv("local_scores.txt", sep="\t")

    # Traceback
    alignment = get_traceback(global_scores, seq1, seq2)

    # Print aligned sequences
    alignment = formatted_alignment(alignment, blosum)

    # Output alignment
    with open(args.output, "w") as OUT:
        print("Score: {0}".format(global_scores[len(seq1),len(seq2)]), file=OUT)
        print("".join([a for a,_,_ in alignment]), file=OUT)
        print("".join([b for _,b,_ in alignment]), file=OUT)
        print("".join([c for _,_,c in alignment]), file=OUT)

if __name__ == '__main__':
    # Command line options
    arguments = get_options()

    # Running script
    main(arguments)
