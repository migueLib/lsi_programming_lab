import os
import argparse
from fastatools import *
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

    arg = parser.parse_args()

    # Standardize paths
    arg.input = os.path.abspath(arg.input)
    arg.output = os.path.abspath(arg.output)

    return arg


def get_cdna(header, sequence):
    """
    Get's the complementary string of a FASTA tuple.

    :param header: header of the FASTA sequence
    :param sequence: sequence of the FASTA sequence
    :return: tuple with the modified header and complementary sequence
    """
    translator = {"A": "T", "T": "A", "C": "G", "G": "C"}
    header = header+"|cDNA"
    sequence = "".join([translator[c] for c in sequence[::-1]])

    return header, sequence


def main(args):
    # Opening input and output files
    output = open(args.output, "w")
    intput = open(args.input, "r")

    # Import data
    fasta = fasta_sequences(intput)
    for hd, seq in fasta:
        hd_c, seq_c = get_cdna(hd, seq)
        write_fasta(output, hd_c, seq_c)

    # Closing input and output files
    output.close()
    intput.close()


if __name__ == '__main__':
    # Command line options
    parameters = get_options()

    # Running script
    main(parameters)
