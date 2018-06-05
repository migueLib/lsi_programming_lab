import os
import re


def single_fasta_sequence(filename):
    """
    Reads fasta file with only one entry, returns duple of (header,sequence)

    :param filename: File object of the FASTA file
    :return: duple header,sequence for the given FASTA file
    """
    seq = ""
    hd = ""

    for line in filename:
        assert isinstance(line, str)
        line = line.strip()
        if line.startswith(">"):
            hd = line[1:]
        else:
            seq += line

    return hd, seq


def fasta_list(filename):
    """
    Reads FASTA file and returns a list of duples (header, sequence)
    :param filename: File object of the FASTA file.
    :return: list of duples containing (header, sequence) structure
    """
    fastas = list()
    hd = ""
    seq = ""

    for line in filename:
        assert isinstance(line, str)
        line = line.strip()
        if line.startswith(">"):
            if hd and seq:
                duple = (hd, seq)
                fastas.append(duple)
            hd = line[1:]
            seq = ""
        else:
            seq += line

    duple = (hd, seq)
    fastas.append(duple)

    return fastas

# path_fasta = os.path.abspath("handout_05/ecoli-genome-sample.fna")
# with open(path_fasta, "r") as FASTA:
#    header, sequence = single_fasta_sequence(FASTA)


path_fasta = os.path.abspath("handout_05/ecoli-genes-non-standard.ffn")
with open(path_fasta, "r") as FASTA:
    headers = fasta_list(FASTA)
