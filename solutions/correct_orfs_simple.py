import re
import sys


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
    dropped = dict()
    keepped = dict()
    for i, j in positions.items():
        if abs(i-j) >= n:
            keepped.setdefault(i, j)
        else:
            dropped.setdefault(i, j)

    return keepped, dropped


def get_best_acc_threshold(orfs, genes):
    max_orf = max(orfs.items(), key=lambda x: abs(x[0]-x[1]))
    max_len = abs(max_orf[0] - max_orf[1])

    results = list()
    for i in range(0, max_len, 3):
        # Filter
        keep, drop = filter_by_len(orfs, i)

        # Accuracy
        results.append([get_balanced_accuracy(set(genes.keys()), set(keep.keys()), set(drop.keys())), i])

    acc, i = max(results)
    keep, drop = filter_by_len(orfs, i)

    return keep, acc, i


def get_stats(orfs, genes, threshold, acc):
    """
    Get statistics for orfs and genes relation

    :param orfs: dictionary with orfs sequence positions
    :param genes: dictionary with genes sequence positions
    :param threshold: threshold
    :param acc: accuracy value
    :return: dictionary with
    """
    # Get length
    n_orfs = len(orfs)
    n_genes = len(genes)

    # Get correct start/end ORFs
    n_correct = len(set(orfs.items()) & set(genes.items()))
    r_correct = n_correct/n_orfs

    # Get correct end ORFs
    n_ends = len(set(orfs.keys() & set(genes.keys())))
    r_ends = n_ends/n_orfs

    # Missed genes
    n_mgenes = n_genes - n_ends

    # Print it!
    print("------------------------------------------------------")
    print("Threshold: {0}".format(threshold))
    print("Accuracy: {0}".format(acc))
    print("Number of ORFs: {0}".format(n_orfs))
    print("Number of genes: {0}".format(n_genes))
    print("ORFs correctly predicted: {0}".format(n_correct))
    print("ORFs correctly predicted (ratio): {0}".format(r_correct))
    print("ORFs ends correctly predicted: {0}".format(n_ends))
    print("ORFs ends correctly predicted (ratio): {0}".format(r_ends))
    print("Missed genes predicted: {0}".format(n_mgenes))


def get_balanced_accuracy(OP, PP, PN):
    """
    Calculates the balanced accuracy for 3 sets
    :param OP: Real positives
    :param PP: Predicted positives
    :param PN: Predicted negatives
    :return: balanced accuracy
    """
    TP = len(OP & PP)
    FP = len(PP - OP)
    TN = len(PN - OP)
    FN = len(OP & PN)

    aac_balanced = 0.5 * ((TP/(TP+FN)) + (TN/(TN+FP)))

    return aac_balanced


# TO run in command line,  ORF path, GENE path, Threshold (just a number, in aa, it)
orf_pathway = sys.argv[1]
gen_pathway = sys.argv[2]
threshold = int(sys.argv[3]) * 3

# IF U WANT 2 RUN THE EX05 CHANGE THIS TO True
find_best_threshold = False


# Get positions of ORFs
orfs = get_sequence_positions(orf_pathway)

# Get positions of genes
genes = get_sequence_positions(gen_pathway)


# HERE IS FOR EX05
if find_best_threshold:
    orfs, acc, threshold = get_best_acc_threshold(orfs, genes)

else:
    # Filter scores
    orfs, orfs_dropped = filter_by_len(orfs, threshold)

    # Get accuracy
    acc = get_balanced_accuracy(set(genes.keys()), set(orfs.keys()), set(orfs_dropped.keys()))

# Get stats
get_stats(orfs, genes, threshold, acc)


