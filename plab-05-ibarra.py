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
    with open(filename, "r") as FASTA:
        for line in FASTA:
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
    with open(filename, "r") as FASTA:
        fastas = list()
        hd = ""
        seq = ""

        for line in FASTA:
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


def fasta_sequences(filename):
    """
    Reads FASTA file and returns an iterable containing (header,sequence)
    duples for each sequence on the fasta file.

    :param filename: File object of the FASTA file.
    :return: iterable with  (header, sequence) s tructure
    """
    with open(filename, "r") as FASTA:
        hd = ""
        seq = ""
        for line in FASTA:
            assert isinstance(line, str)
            line = line.strip()
            if line.startswith(">"):
                if hd and seq:
                    yield (hd, seq)
                hd = line[1:]
                seq = ""
            else:
                seq += line

        yield (hd, seq)


def slice_string_by_n(s, n):
    """
    Slice string by a given n (number of characters)

    :param s: String to be sliced
    :param n: Number of characters per slice
    :return: iterable with sliced strings.
    """
    while s:
        yield s[:n]
        s = s[n:]


def write_fasta(f, hd, seq):
    """
    Writes a fasta file from a header and sequence
    :param f: file object
    :param hd: header string
    :param seq: sequence string
    :return: output FASTA file
    """
    assert isinstance(hd, str)
    assert isinstance(seq, str)

    # Split to sequence
    sliced_sequences = slice_string_by_n(seq, 70)

    # Print header
    print(">" + hd, file=f)
    for s in sliced_sequences:
        print(s, file=f)

# Easiest way
# path_fasta = os.path.abspath("handout_05/ecoli-genome-sample.fna")
# header, sequence = single_fasta_sequence(path_fasta)

# Lists to get them all on the same place
# path_fasta = os.path.abspath("handout_05/ecoli-genes-non-standard.ffn")
# headers = fasta_list(path_fasta)

# Generators to make it memory efficient
# path_fasta = os.path.abspath("handout_05/ecoli-genes-non-standard.ffn")
# sequences = fasta_sequences(path_fasta)

# part d): test b) and c) on ecoli-proteome.faa
# path_fasta = os.path.abspath("handout_05/ecoli-proteome.faa")
# max_list = max(fasta_list(path_fasta), key=lambda x: len(x[1]))
# min_list = min(fasta_list(path_fasta), key=lambda x: len(x[1]))
# print(len(max_list[1]), len(min_list[1]))

# max_iter = max(fasta_sequences(path_fasta), key=lambda x: len(x[1]))
# min_iter = min(fasta_sequences(path_fasta), key=lambda x: len(x[1]))
# print(len(max_iter[1]), len(min_iter[1]))

# Writing FASTA sequence


with open("handout_05/output.fasta", "w") as f:
    write_fasta(f, "albatross", "WHATFLAVQRISIT")
    write_fasta(f, "lumberjack", "ISLEEPALLNIGHTANDIWQRKALLDAY")
    write_fasta(f, "deadparrot", "NQRWEGIANPLVE")