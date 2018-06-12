import os
import re
import argparse
from fastatools import *
from argparse import RawDescriptionHelpFormatter


def get_options():
    """
    Function to pull in arguments
    """

    description = """ Evaluates how the predicted ORF's correlate to actual genes1
     """

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawDescriptionHelpFormatter)
    # Standard Input
    standard = parser.add_argument_group(title='Standard input',
                                         description='Standard input for tool')
    standard.add_argument('-o', "--orfs", dest="orfs", action='store',
                          required=True, help="Path for ORFs FASTA")
    standard.add_argument('-g', "--genes", dest="genes", action='store',
                          required=True, help="Path to Genes FASTA.")
    standard.add_argument('-t', "--threshold", dest="threshold", action='store', type=int,
                          required=True, help="Path to Genes FASTA.")
    arg = parser.parse_args()

    # Standardize paths
    arg.orfs = os.path.abspath(arg.orfs)
    arg.genes = os.path.abspath(arg.genes)

    return arg


def get_sequence_positions(fasta_file):
    """
    Get positions from all headers for a given fasta_file
    :param fasta_file: fasta_file path
    :return: dictionary containing end:start for all headers
    """
    dict_pos = dict()
    with open(fasta_file, "r") as fasta:
        for l in fasta:
            if l.startswith(">"):
                regex = re.search(":c?([0-9]+)-([0-9]+)", l)
                dict_pos.setdefault(int(regex.group(2)), int(regex.group(1)))
    return dict_pos


def filter_by_len(positions, n):
    """
    Filters the sequences with less than lenght n
    :param positions: dictionary with sequences positions
    :param n: threshold for filtering
    :return: positions filtered
    """
    return {i: j for i, j in positions.items() if abs(i-j) >= n}


def get_stats(orfs, genes):
    """
    Get statistics for orfs and genes relation

    :param orfs: dictionary with orfs sequence positions
    :param genes: dictionary with genes sequence positions
    :return: dictionary with
    """
    # Get stats
    n_orfs = len(orfs)
    n_genes = len(genes)

    n_correct = len(set(orfs.items()) & set(genes.items()))
    r_correct = n_correct/n_orfs

    n_ends = len(set(orfs.keys() & set(genes.keys())))
    r_ends = n_ends/n_orfs

    n_mgenes = n_genes - n_ends

    return n_orfs, n_genes, n_correct, r_correct, n_ends, r_ends, n_mgenes


def main(args):
    # Get positions of ORFs
    pos_orfs = get_sequence_positions(args.orfs)

    # Get positions of genes
    pos_genes = get_sequence_positions(args.genes)

    # Filter scores
    pos_orfs = filter_by_len(pos_orfs, args.threshold*3)
    pos_genes = filter_by_len(pos_genes, args.threshold*3)

    # Get stats
    n_orfs, n_genes, n_correct, r_correct, n_ends, r_ends, n_mgenes = get_stats(pos_orfs, pos_genes)

    # Print data
    print("Number of ORFs: {0}".format(n_orfs))
    print("Number of genes: {0}".format(n_genes))
    print("The total number of ORFs correctly predicted: {0}".format(n_correct))
    print("The ratio of ORFs correctly predicted: {0}".format(r_correct))
    print("The total number of ORFs ends correctly predicted: {0}".format(n_ends))
    print("The ratio of ORFs ends correctly predicted: {0}".format(r_ends))
    print("The ratio of missed genes predicted: {0}".format(n_mgenes))


if __name__ == '__main__':
    # Command line options
    parameters = get_options()

    # Running script
    main(parameters)
