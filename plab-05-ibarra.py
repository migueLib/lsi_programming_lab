from fastatools import *
import os

# Easiest way
path_fasta = os.path.abspath("handout_05/ecoli-genome-sample.fna")
with open(path_fasta, "r") as pf:
    header, sequence = single_fasta_sequence(pf)
    print(header)

# Lists to get them all on the same place
path_fasta = os.path.abspath("handout_05/ecoli-genes-non-standard.ffn")
with open(path_fasta, "r") as pf:
    headers = fasta_list(pf)
    print(headers)

# Generators to make it memory efficient
path_fasta = os.path.abspath("handout_05/ecoli-genes-non-standard.ffn")
with open(path_fasta, "r") as pf:
    sequences = fasta_sequences(pf)
    print(sequences)

# part d): test b) and c) on ecoli-proteome.faa
path_fasta = os.path.abspath("handout_05/ecoli-proteome.faa")
with open(path_fasta, "r") as pf:
    max_list = max(fasta_list(pf), key=lambda x: len(x[1]))
with open(path_fasta, "r") as pf:
    min_list = min(fasta_list(pf), key=lambda x: len(x[1]))
print(len(max_list[1]), len(min_list[1]))

with open(path_fasta, "r") as pf:
    max_iter = max(fasta_sequences(pf), key=lambda x: len(x[1]))
with open(path_fasta, "r") as pf:
    min_iter = min(fasta_sequences(pf), key=lambda x: len(x[1]))
print(len(max_iter[1]), len(min_iter[1]))

# Writing FASTA sequence
with open("handout_05/output.fasta", "w") as f:
    write_fasta(f, "albatross", "WHATFLAVQRISIT")
    write_fasta(f, "lumberjack", "ISLEEPALLNIGHTANDIWQASDFASDFMALMSDFLAMSDFNKAVAKSDKALDSKFADFLKASRKALLDASLEEPALLNIGHTANDIWQASDFASDFMALMSDFLAMSDFNKAVAKSDKALDSKFASDFLKASRKALLDASLEEPALLNIGHTANDIWQASDFASDFMALMSDFLAMSDFNKAVAKSDKALDSKFASDFLKASKALLDASLEEPALLNIGHTANDIWQASDFASDFMALMSDFLAMSDFNKAVAKSDKALDSKFASDFLKASRKALLDAY")
    write_fasta(f, "deadparrot", "NQRWEGIANPLVE")


# Creating  module: Putting things together
with open("handout_05/ecoli-genome.fna") as f:
    species, genome = single_fasta_sequence(f)
    print("The genome of", species, "contains", len(genome), "nucleotides.")

with open("handout_05/truth.faa", "w") as f:
    write_fasta(f,"theking","ELVISISALIVE")
    write_fasta(f,"liverpoolfour","PAVLISDEAD")