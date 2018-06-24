import os
import numpy as np
import pandas as pd
from itertools import combinations


def read_distance_matrix(path):
    """
    Reads file and converts it to a dictionary with nodes as
    tuples and distances as values
    :param path: path to open the file
    :return: dictionary
    """
    distances = pd.read_csv(path, sep="\s+", index_col=0, header=0, engine="python")
    distances_dict = distances.to_dict()
    #tuples = combinations(distances.columns.tolist(), 2)
    #matrix = {(i, j): int(distances[i][j]) for i, j in tuples}
    #print(distances_dict)

    return distances_dict


def get_n_elements(tup):
    """
    Gets the number of elements of a given tuple
    :param tuple: tuple with elements
    :return: number of elements in a given tuple
    """
    n = 0
    for elem in tup:
        if type(elem) == type(tuple()):
            nl = get_n_elements(elem)
            n += nl
        if type(elem) == type(str()):
            n += 1
    return n


def get_closest(matrix):
    """
    Get the clusters keys with the smallest distances among them
    :param matrix: distance matrix
    :return: cluster with the smallest distance
    """
    return min(matrix.items(), key= lambda x:x[1])[0]


def merge(matrix, closest):
    """
    Re-calculates distance matrix based on a pair of values
    :param matrix: matrix distance
    :param closest: closest clusters in matrix
    :return: update matrix
    """
    # Removing min key from elements to_change (this will avoid problems later on)
    del matrix[closest]

    # Subset matrix into to_changed and unchanged
    to_change = {(i, j): matrix[(i, j)] for i, j in matrix if i in closest or j in closest}
    unchanged = ({(i, j): matrix[(i, j)] for i, j in matrix if i not in closest and j not in closest})

    # Get single letters (except the ones in min)
    singles = set()
    for i,j in matrix.keys():
        if i not in closest:
            singles.update(i)
        if j not in closest:
            singles.update(j)

    # Re calculate distance to merged values in to_changed values
    print(closest)
    print(singles)
    print(matrix.keys())
    print(to_change.keys())
    print(unchanged.keys())

    for i in singles:
        print((i, closest[0]),(i, closest[1]))
        try:
            val_a = to_change[(closest[0], i)]
        except KeyError:
            val_a = to_change[(i, closest[0])]
        try:
            val_b = to_change[(closest[1], i)]
        except KeyError:
            val_b = to_change[(i, closest[1])]

        # Calculates |R|,|S| and |R|+|S|
        n = get_n_elements(closest)
        r = get_n_elements(closest[0])
        s = get_n_elements(closest[1])
        unchanged[(closest, i)] = (r*val_a + s*val_b)/n

    return unchanged


# path_matrix = "handout_06/small-distances.txt"
path_matrix = "handout_06/wiki"
distance_matrix = read_distance_matrix(path_matrix)
smallest = get_closest(distance_matrix)
first = merge(distance_matrix, smallest)
# Second
print(first)
s2 = get_closest(first)
second = merge(first, s2)
print(second)
s3 = get_closest(second)
third = merge(first, s3)
print(third)