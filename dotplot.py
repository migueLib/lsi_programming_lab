#!/usr/bin/env python
################################################################################
# DATE: 2018/05/18
#
# SCRIPT: dotplot.py
#
# VERSION: 1.0
#
# AUTHOR: Miguel A. Ibarra-Arellano (ibarrarellano@gmail.com)
#
# DESCRIPTION:
#
# The output is a set of graphs and spreadsheets of flags
#
################################################################################
# Import native libraries
import os
import argparse
from argparse import RawDescriptionHelpFormatter

# Import add-on libraries
import matplotlib.pyplot as plt
import numpy as np


def getOptions():
    """
    Function to pull in arguments
    """

    description = """ Creates a dotplot either on ASCII art, or matplotlib """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawDescriptionHelpFormatter)
    # Standard Input
    standard = parser.add_argument_group(title='Standard input',
                                         description='Standard input for SECIM tools.')
    standard.add_argument('-sA', "--seqA", dest="seqA", action='store',
                          required=True, help="Input sequence A.")
    standard.add_argument('-sB', "--seqB", dest="seqB", action='store',
                          required=True, help="Input sequence B.")
    standard.add_argument('-w', "--window", dest="window", action='store',
                          required=True, type=int, help="Window size.")
    standard.add_argument('-s', "--stringency", dest="stringency", action='store',
                          required=True, type=int, help="Stringency.")
    standard.add_argument('-t', "--title", dest="title", action='store',
                          required=True, help="dotplot title.")
    standard.add_argument('-o', "--output", dest="output", action='store',
                          required=True, help="Output dataset in wide format.")

    args = parser.parse_args()

    # Standardize paths
    args.output = os.path.abspath(args.output)

    return args


def fasta2str(file):
    """
    Takes a fasta file and converts it to a string
    :param file: file path
    :return: fasta string
    """
    fasta_string = list()
    with open(file, "r") as IN:
        for line in IN:
            if not line.startswith(">"):
                line = line.strip()
                fasta_string.append(line)
        fasta_string = "".join(fasta_string)
    return fasta_string


def dotplot(seqA, seqB, w=1, s=1):
    """
    Generates a dotplot matrix for 2 given strings

    :param seqA: Str A
    :param seqB: Str B
    :param w: int window size
    :param s: int stringency
    :param heading: heading of the file
    :param filename: filename to output the doplot matrix
    :return:
    """
    # Create empty matrix
    dp = np.zeros((len(seqA), len(seqB)), dtype=int)

    # Calcualte window distance
    dis = (w - 1) / 2

    # Lets fill every single point in the matrix
    for i in range(len(seqA)):  # Rows
        for j in range(len(seqB)):  # Columns
            subsA = seqA[int(i - dis):int(i + dis + 1)]
            subsB = seqB[int(j - dis):int(j + dis + 1)]

            # Check every element of the substring for stringency
            n_matches = 0
            for k in range(min(len(subsA), len(subsB))):
                if subsA[k] == subsB[k]:
                    n_matches += 1

            # if the number of matches is bigger than the stringency then put a 1
            # print(n_matches,s)
            if n_matches > s:
                dp[i, j] = 1
    return (dp)


def dotplot2graphics(dp, hdA, hdB, heading="dotplot", filename="my_dotplot.pdf", use_imshow=True):
    """
    Prints dotplot with matplotlib

    :param dp: dotplot matrix
    :param hdA: string A
    :param hdB: string B
    :param heading:  heading
    :param filename: output filename
    :return: file with a dotplot matrix
    """

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    if use_imshow:
        ax.imshow(dp)
    else:
        for i in range(len(seqA)):
            for j in range(len(seqB)):
                if dp[i, j] == 1:
                    ax.scatter(j, i, color="black", marker="+")

        # Invert axis
        plt.gca().invert_yaxis()

    # Mover xticks to top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    # Set xticks to letters
    plt.xticks(np.arange(len(hdB)), list(hdB))
    plt.yticks(np.arange(len(hdA)), list(hdA))

    # Set title
    plt.title(heading, y=1.08)

    # Save figure
    plt.savefig(filename)


def main(args):
    # Convert fasta files to strings
    seq_a = fasta2str(args.seqA)
    seq_b = fasta2str(args.seqB)

    # Create dotplot matrix
    dp = dotplot(seqA=seq_a, seqB=seq_b, w=args.window, s=args.stringency)

    # Create dotplot image
    dotplot2graphics(dp, seq_a, seq_b, heading=args.title, filename=args.output)


if __name__ == '__main__':
    # Command line options
    args = getOptions()

    # Running script
    main(args)

