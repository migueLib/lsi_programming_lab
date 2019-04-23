import pandas as pd


def read_distance_matrix(path):
    """
    Reads file and converts it to a dictionary with nodes as
    tuples and distances as values
    :param path: path to open the file
    :return: dictionary
    """
    distances = pd.read_csv(path, sep="\s+", index_col=0, header=0, engine="python")
    distances = distances.to_dict()

    # Remove the diagonal from the distance matrix
    distances_dict = dict()
    for i in distances.keys():
        distances_dict[i] = dict()
        for j in distances.keys():
            if i != j:
                distances_dict[i][j] = distances[i][j]

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
    # Transform {key:{key:value}} to  {(key,key):value}
    dis_matrix = dict()
    for i in matrix.keys():
        for j in matrix[i]:
            dis_matrix[(i, j)] = matrix[i][j]

    return min(dis_matrix.items(), key=lambda x: x[1])


def merge(matrix, clusters):
    """
    Re-calculates distance matrix based on a pair of values
    :param matrix: matrix distance
    :param clusters: closest clusters in matrix
    :return: update matrix
    """
    # Removing min key from elements to_change (this will avoid problems later on)
    del matrix[clusters[0]][clusters[1]]
    del matrix[clusters[1]][clusters[0]]

    # Taking unchanged distances (the ones not included in the closest nodes
    new_distances = dict()
    c_to_merge = dict()

    # Splitting matrix into the keys that are going to remain intact\
    #  and the distances that need to be updated
    cleaned = dict()
    for i in matrix.keys():
        cleaned[i] = dict()
        for j in matrix[i].keys():
            if j not in clusters:
                cleaned[i][j] = matrix[i][j]

    # Split cleaned matrix between the ones to merge and the ones that should not be touched
    for i in cleaned:
        if i in clusters:
            c_to_merge[i] = cleaned[i]
        else:
            new_distances[i] = cleaned[i]

    # make the formula
    keys_not_in_cluster = list(new_distances.keys())


    new_distances[clusters] = dict()

    for k in keys_not_in_cluster:
        # Get val_r
        try:
            val_r = cleaned[clusters[0]][k]
        except KeyError:
            val_r = cleaned[k][clusters[0]]

        # Get vel2
        try:
            val_s = cleaned[clusters[1]][k]
        except KeyError:
            val_s = cleaned[k][clusters[1]]

        # Calculates |R|,|S| and |R|+|S|
        n = get_n_elements(clusters)
        r = get_n_elements(clusters[0])
        s = get_n_elements(clusters[1])
        new_distances[k][clusters] = (r*val_r + s*val_s)/n
        new_distances[clusters][k] = (r*val_r + s*val_s)/n

    return new_distances


def hierarchical_clustering(matrix):
    """
    Calculates HC Hierarchical clustering for a matrix distance

    :param matrix: distance matrix
    :return:
    """
    heights=dict()
    n_nodes = len(matrix)
    n_cluster = 1
    while n_cluster < n_nodes:
        clusters = get_closest(matrix)
        matrix = merge(matrix, clusters[0])
        n_cluster = get_n_elements(clusters)
        heights[clusters[0]] = clusters[1]

    return heights


# path_matrix = "handout_06/small-distances.txt"
# Reading distance matrix
path_matrix = "handout_06/wiki"
distance_matrix = read_distance_matrix(path_matrix)

# Get hierarchical clustering for distance matrix
hc = hierarchical_clustering(distance_matrix)
print(hc)