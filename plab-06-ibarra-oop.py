import pandas as pd

class DistanceMatrix:
    def __init__(self, path):
        distances = pd.read_csv(path, sep="\s+", index_col=0, header=0, engine="python")
        distances = distances.to_dict()

        # Remove the diagonal from the distance matrix
        self.matrix = dict()
        for i in distances.keys():
            self.matrix[i] = dict()
            for j in distances.keys():
                if i != j:
                    self.matrix[i][j] = distances[i][j]


    def get_closest(self):
        dis_matrix = dict()
        for i in self.matrix.keys():
            for j in self.matrix[i].keys():
                A = self.get_n_elements(i)
                B = self.get_n_elements(j)
                dis_matrix[(i, j)] = self.matrix[i][j]

        self.min = min(dis_matrix.items(), key=lambda x: x[1])
        self.min_key, self.min_value = self.min[0], self.min[1]

    def get_n_elements(self, tup):
        n = 0
        for elem in tup:
            if type(elem) == type(tuple()):
                nl = self.get_n_elements(elem)
                n += nl
            if type(elem) == type(str()):
                n += 1
        return n

    def merge(self):
        # Calculating
        self.get_closest()

        # Removing min key from elements to_change (this will avoid problems later on)
        del self.matrix[self.min_key[0]][self.min_key[1]]
        del self.matrix[self.min_key[1]][self.min_key[0]]

        # Taking unchanged distances (the ones not included in the closest nodes
        new_distances = dict()
        distances_to_merge = dict()

        # Splitting matrix into the keys that are going to remain intact\
        #  and the distances that need to be updated
        cleaned = dict()
        for i in self.matrix.keys():
            cleaned[i] = dict()
            for j in self.matrix[i].keys():
                if j not in self.min_key:
                    cleaned[i][j] = self.matrix[i][j]

        # Split cleaned matrix between the ones to merge and the ones that should not be touched
        for i in cleaned:
            if i in self.min_key:
                distances_to_merge[i] = cleaned[i]
            else:
                new_distances[i] = cleaned[i]

        # make the formula
        keys_not_in_cluster = list(new_distances.keys())

        new_distances[self.min_key] = dict()

        for k in keys_not_in_cluster:
            # Get val_r
            try:
                val_r = cleaned[self.min_key[0]][k]
            except KeyError:
                val_r = cleaned[k][self.min_key[0]]

            # Get vel2
            try:
                val_s = cleaned[self.min_key[1]][k]
            except KeyError:
                val_s = cleaned[k][self.min_key[1]]

            # Calculates |R|,|S| and |R|+|S|
            r = self.get_n_elements(self.min_key[0])
            s = self.get_n_elements(self.min_key[1])
            new_distances[k][self.min_key] = (r*val_r + s*val_s) / (r+s)
            new_distances[self.min_key][k] = (r*val_r + s*val_s) / (r+s)

        self.matrix = new_distances

    def hierarchical_clustering(self):
        self.heights = dict()
        n_nodes = len(self.matrix)
        n_cluster = 1

        while n_cluster < n_nodes:
            self.get_closest()
            self.heights[self.min_key] = self.min_value
            self.merge()
            n_cluster = self.get_n_elements(self.min)

    def print_dendrogram(self, sep=3):
        def is_pair(T):
            return type(T) == tuple and len(T) == 2

        def max_height(T):
            if is_pair(T):
                h = max(max_height(T[0]), max_height(T[1]))
            else:
                h = len(str(T))
            return h + sep

        # Initialize active levels
        active_levels = {}

        def traverse(T, h, is_first):
            # Get separations for the pairs (start)
            if is_pair(T):
                traverse(T[0], h - sep, 1)
                string = [' '] * (h - sep)
                string.append('|')
            else:
                string = list(str(T))
                string.append(' ')

            # Get heights
            while len(string) < h:
                string.append('-')

            # Get corners
            if is_first >= 0:
                string.append('+')
                if is_first:
                    active_levels[h] = 1
                else:
                    del active_levels[h]

            A = list(active_levels)
            A.sort()
            for L in A:
                if len(string) < L:
                    while len(string) < L:
                        string.append(' ')
                    # Get separations for the pairs of values (end)
                    string.append('|')

            print(''.join(string))

            if is_pair(T):
                traverse(T[1], h - sep, 0)

        traverse(self.min_key, max_height(self.min_key), -1)


path_matrix = "handout_06/wiki"
dm = DistanceMatrix(path_matrix)
dm.hierarchical_clustering()
print(dm.heights)
dm.print_dendrogram()

