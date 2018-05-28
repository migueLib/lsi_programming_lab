import os
import argparse
from math import log
from collections import Counter
from itertools import combinations
from itertools import combinations_with_replacement as cwr
from argparse import RawDescriptionHelpFormatter
import numpy as np
import pandas as pd


def get_options():
    """
    Function to pull in arguments
    """

    description = """ Creates an scoring matrix for a given alignment """
    
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawDescriptionHelpFormatter)
    # Standard Input
    standard = parser.add_argument_group(title='Standard input',
                                         description='Standard input for SECIM tools.')
    standard.add_argument('-i', "--input", dest="input", action='store',
                          required=True, help="Path for input alignment")
    standard.add_argument('-o', "--output", dest="output", action='store',
                          required=True, help="Path to output matrix in wide format.")

    args = parser.parse_args()

    # Standardize paths
    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)

    return args


def get_frequency(sequences):
    """
    Counts the frequency of an aa in a sequence

    :param sequences: string
    :return:
    """
    frequency = Counter()
    for seq in sequences:
        frequency.update({l: seq.count(l) for l in set(seq)})

    return dict(frequency)


def get_probability(letters, n):
    """
    Calculates the probability of a frequency based on the number of elements

    :param letters: dictionary with frequencies
    :param n: total number of elements
    :return: dictionary with probabilities
    """
    return {l: c/n for l, c in letters.items()}


def get_pair_frequency(sequences):
    """

    :param sequences: list of sequences
    :return: dictionary of pairs
    """
    pair_count = dict()
    pairs = combinations(sequences, 2)
    for seq1, seq2 in pairs:
        for l1, l2 in zip(seq1, seq2):
            try:
                pair_count["".join(sorted([l1, l2]))] += 1
            except KeyError:
                pair_count["".join(sorted([l1, l2]))] = 1
    return pair_count


def get_expected_probability(probabilities):
    """
    Calculates the expected probability of two variables
    :param probabilities: dictionary with all the probabilites per element
    :return: expected probabilites per pair of values
    """

    expected = dict()
    for a, b in cwr(probabilities.keys(), 2):
        if a == b:
            expected["".join(sorted([a, b]))] = probabilities[a] * probabilities[b]
        else:
            expected["".join(sorted([a, b]))] = 2 * (probabilities[a] * probabilities[b])

    return expected


def get_log_odds_score(observed, expected):
    """
    Get logs ratios for all pairs of values
    :param observed: Observed values
    :param expected: Expected values
    :return:
    """
    log_ratio = dict()
    for pe in expected.keys():
        try:
            log_ratio[pe] = int(round(2*log(observed[pe]/expected[pe], 2),0))
        except KeyError:
            log_ratio[pe] = int(-99)

    return log_ratio


def read_file(path):
    """
    Function to read input file

    :param path: path to the file
    :return: list of paths
    """
    with open(path, "r") as IN:
        file_seqs = [line.strip() for line in IN]
    return file_seqs


def output_matrix(ratios, path):
    """
    Outputs substitution matrix

    :param ratios: dictionary with log_odds_ratios  for each pair
    :param path: output file path
    :return: Substitution matrix
    """
    # Getting amino acids
    aa = sorted(list(set([list(i)[0] for i in ratios.keys()])))

    # Creating empty data frame
    my_df = pd.DataFrame(data=np.zeros((len(aa),len(aa))), index=aa, columns=aa)

    # Filling df
    for key in ratios.keys():
        a, b = list(key)
        my_df[a][b] = ratios[key]
        my_df[b][a] = ratios[key]

    my_df.to_csv(path, sep="\t")


def main(args):
    # Reading sequences from file
    sequences = read_file(args.input)
    #sequences = ["TSVKTYAKFVTH", "TSVKTYAKFSTH", "TSVKTYAKFVTH", "LSVKKYPKYVVQ", "SSVKKYPKYSVL"]

    # Getting element frequency fa
    fa = get_frequency(sequences)
    print(fa)

    # Getting number of elements
    n = sum(fa.values())
    print(n)

    # Get element probability
    pa = get_probability(fa, n)
    print(pa)

    # Get Pair frequency
    fab = get_pair_frequency(sequences)
    print(fab)

    # Get number of pairs
    n_p = sum(fab.values())
    print(n_p)

    # Get observed pair probability
    pab = get_probability(fab, n_p)
    print(pab)

    # Get expected pair probability
    eab = get_expected_probability(pa)
    print(eab)

    # Get logs ratio
    sab = get_log_odds_score(pab, eab)
    print(sab)

    # Output substitution matrix
    output_matrix(sab, args.output)


if __name__ == '__main__':
    # Command line options
    args = get_options()

    # Running script
    main(args)
