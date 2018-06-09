import os
import argparse
from argparse import RawDescriptionHelpFormatter


def get_options():
    """
    Function to pull in arguments
    """

    description = """ Gets the complementary dna string of a fasta file """

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawDescriptionHelpFormatter)
    # Standard Input
    standard = parser.add_argument_group(title='Standard input',
                                         description='Standard input for tool')
    standard.add_argument('-i', "--input", dest="input", action='store',
                          required=True, help="Path for input FASTA")
    standard.add_argument('-o', "--output", dest="output", action='store',
                          required=True, help="Path to output FASTA.")

    args = parser.parse_args()

    # Standardize paths
    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)

    return args


def get_complementary_fasta(input, output):
    """
    Get's the complementary string of a FASTA file

    :param input: FASTA file input pathway
    :param output: FASTA file output pathway
    :return: Translated FASTA file
    """
    assert isinstance(input, str)
    assert isinstance(output, str)

    translator = {"A":"T", "T":"A", "C":"G", "G":"C"}

    INPUT = open(input, "r")
    OUTPUT = open(output, "w")

    for line in INPUT:
        line = line.strip()
        if line.startswith(">"):
            print(line+"|cDNA", file=OUTPUT)
        else:
            print("".join([translator[c] for c in line]), file=OUTPUT)


def main(args):
    # Convert fasta file
    get_complementary_fasta(args.input, args.output)


if __name__ == '__main__':
    # Command line options
    args = get_options()

    # Running script
    main(args)

