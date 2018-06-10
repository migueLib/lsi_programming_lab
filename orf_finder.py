import os
import argparse
from fastatools import *
from cdna import *
from argparse import RawDescriptionHelpFormatter
# from answerstolife import solutiontothisexercise


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


def get_orf(hd, seq, f):
    # Setting up stop and start codons
    stop = ["TAA", "TAG", "TGA"]
    start = "ATG"

    # Get complementary DNA
    hd_c, seq_c = get_cdna(hd, seq)

    # Initialize index
    orf1_pos, orf2_pos, orf3_pos = [0, 0], [0, 0], [0, 0]
    orf4_pos, orf5_pos, orf6_pos = [0, 0], [0, 0], [0, 0]


    # Initialize empty ORFs sequences
    orf1_seq, orf2_seq, orf3_seq = "", "", ""
    orf4_seq, orf5_seq, orf6_seq = "", "", ""

    # Iterating over forward and reversed sequences
    for i in range(0, len(seq), 3):
        # Indices for 3 frames
        i, j, k = i, i+1, i+2

        # Forward
        fr1 = seq[i:i+3]
        fr2 = seq[j:j+3]
        fr3 = seq[k:k+3]

        # Reverse
        fr4 = seq_c[i:i+3]
        fr5 = seq_c[j:j+3]
        fr6 = seq_c[k:k+3]

        # FORWARD
        # ORF 1
        if fr1 == start and not orf1_seq:
            orf1_seq = fr1
            orf1_pos[0] = i+1
        elif orf1_seq:
            orf1_seq += fr1

        # ORF 2
        if fr2 == start and not orf2_seq:
            orf2_seq = fr2
            orf2_pos[0] = j+1
        elif orf2_seq:
            orf2_seq += fr2

        # ORF 3
        if fr3 == start and not orf2_seq:
            orf3_seq = fr3
            orf3_pos[0] = k+1
        elif orf3_seq:
            orf3_seq += fr3

        # REVERSE
        # ORF 4
        if fr4 == start and not orf4_seq:
            orf4_seq = fr4
            orf4_pos[0] = i+1
        elif orf4_seq:
            orf4_seq += fr4

        # ORF 2
        if fr5 == start and not orf5_seq:
            orf5_seq = fr5
            orf5_pos[0] = j+1
        elif orf5_seq:
            orf5_seq += fr5

        # ORF 3
        if fr6 == start and not orf6_seq:
            orf6_seq = fr6
            orf6_pos[0] = k+1
        elif orf6_seq:
            orf6_seq += fr6

        # Print founded sequences
        if fr1 in stop and orf1_seq:
            orf1_pos[1] = i+3
            write_fasta(f, hd+":{0}-{1}".format(orf1_pos[0], orf1_pos[1]), orf1_seq)
            orf1_seq = ""
        if fr2 in stop and orf2_seq:
            orf2_pos[1] = j+3
            write_fasta(f, hd+":{0}-{1}".format(orf2_pos[0], orf2_pos[1]), orf2_seq)
            orf2_seq = ""
        if fr3 in stop and orf3_seq:
            orf3_pos[1] = k+3
            write_fasta(f, hd+":{0}-{1}".format(orf3_pos[0], orf3_pos[1]), orf3_seq)
            orf3_seq = ""
        if fr4 in stop and orf4_seq:
            orf4_pos[1] = i+3
            write_fasta(f, hd_c+":c{0}-{1}".format(len(seq_c)-orf4_pos[0], len(seq_c)-orf4_pos[1]), orf4_seq)
            orf4_seq = ""
        if fr5 in stop and orf5_seq:
            orf5_pos[1] = j+3
            write_fasta(f, hd_c+":c{0}-{1}".format(len(seq_c)-orf5_pos[0], len(seq_c)-orf5_pos[1]), orf5_seq)
            orf5_seq = ""
        if fr6 in stop and orf6_seq:
            orf6_pos[1] = k+3
            write_fasta(f, hd_c+":c{0}-{1}".format(len(seq_c)-orf6_pos[0], len(seq_c)-orf6_pos[1]), orf6_seq)
            orf6_seq = ""



def main(args):

    # Load fasta file
    fasta = open(args.input, "r")
    output = open(args.output, "w")

    # Reading fasta sequences
    sequences = fasta_sequences(fasta)

    # Get ORFs
    for hd, seq in sequences:
        get_orf(hd, seq, output)

    # Close files
    fasta.close()
    output.close()

if __name__ == '__main__':
    # Command line options
    parameters = get_options()

    # Running script
    main(parameters)
