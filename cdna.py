# Ok, here i don't use fastatools module because the will involve loading whole sequences
# of FASTA sequences to memory and we can do it more memory efficient!

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
    header = header+"|cDNA"
    sequence = "".join([translator[c] for c in sequence])
    fasta_tuple = (header, sequence)

    return fasta_tuple


def get_complementary_fasta(fasta, output):
    """
    Get's the complementary string of a FASTA file

    :param fasta: FASTA file input pathway
    :param output: FASTA file output pathway
    :return: Translated FASTA file
    """
    assert isinstance(fasta, str)
    assert isinstance(output, str)

    translator = {"A": "T", "T": "A", "C": "G", "G": "C"}

    input_f = open(fasta, "r")
    output_f = open(output, "w")

    for line in input_f:
        line = line.strip()
        if line.startswith(">"):
            print(line+"|cDNA", file=output_f)
        else:
            print("".join([translator[c] for c in line]), file=output_f)


def main(args):
    # Convert fasta file
    get_complementary_fasta(args.input, args.output)


if __name__ == '__main__':
    # Command line options
    parameters = get_options()

    # Running script
    main(parameters)
