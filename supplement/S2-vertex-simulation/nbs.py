from collections import defaultdict

import networkx as nx
import numpy as np
from hyppo.independence import Dcorr
from joblib import Parallel, delayed
from scipy.stats import ttest_ind


def test(X, Y, n_nodes, test, threshold):

    nbs = NBS(test, threshold, X, Y)
    components = nbs.compute()
    nbs_pvalue = {node: pvalue for pvalue, nodes in components for node in list(nodes)}

    pvalues = defaultdict(list)
    for node in range(n_nodes):
        try:
            pvalues[node].append(nbs_pvalue[node])
        except:
            pvalues[node].append(1)

    return pvalues


def ttest_statistic(sample):
    return ttest_ind(*sample, equal_var=False)[0]


def dcorr_statistic(sample):
    sample = [X.reshape(-1, 1) for X in sample]
    return Dcorr().statistic(*sample)


def get_number_of_edges(G, component):
    return G.subgraph(component).number_of_edges()


class NBS:
    def __init__(self, test, threshold, X, Y):
        self.test = test
        self.threshold = threshold
        self.X = X
        self.Y = Y

        # Get number of vertices in graphs
        _, Ix, Jx = self.X.shape
        _, Iy, Jy = self.Y.shape
        assert Ix == Jx == Iy == Jy
        self.n_vertices_ = Ix

    def get_thresholded_test_statistics(self, X, Y):

        # Matrix to store test statistics for each edge
        stat_matrix = np.zeros((self.n_vertices_, self.n_vertices_))

        # Compute a test statistic for each edge
        indices = zip(*np.triu_indices(self.n_vertices_, 1))
        for i, j in indices:
            sample = [X[:, i, j], Y[:, i, j]]
            try:
                stat = np.abs(self.test(sample))
                if np.isnan(stat):
                    raise ValueError
            except ValueError:
                stat = -np.inf
            stat_matrix[i, j] = stat

        # Treshold the test statistics
        tresholded_edges = np.zeros((self.n_vertices_, self.n_vertices_))
        if self.test == ttest_statistic:
            rel_edges = np.where(stat_matrix >= self.threshold)
        elif self.test == dcorr_statistic:
            rel_edges = np.where(
                np.logical_and(0 < stat_matrix, stat_matrix <= self.threshold)
            )
        tresholded_edges[rel_edges] = 1
        tresholded_edges = tresholded_edges + tresholded_edges.T

        return tresholded_edges

    def get_null_distribution(self):

        # Permute the order of the graphs
        graphs = np.vstack([self.X, self.Y])
        n_samples, _, _ = graphs.shape
        permuted_index = np.random.permutation(n_samples)
        permuted_graphs = graphs[permuted_index]
        X, Y = permuted_graphs[: n_samples // 2], permuted_graphs[n_samples // 2 :]

        # Compute test statistic for each edge and threshold
        stat_matrix = self.get_thresholded_test_statistics(X, Y)

        # Get the size of the largest connected component
        G = nx.from_numpy_array(stat_matrix)
        largest_cc = max(nx.connected_components(G), key=len)
        return get_number_of_edges(G, largest_cc)

    def compute(self, n_reps=200):
        """
        Compute two-tailed network-based statistic (NBS). X, Y are stacks of
        graphs with shape (n_samples, n_vertices_, n_vertices_).
        """

        # Perform a two-sample test at each edge independently and threshold
        stat_matrix = self.get_thresholded_test_statistics(self.X, self.Y)

        # Identify connectected components with at least 1 edge
        G = nx.from_numpy_array(stat_matrix)
        components = list(nx.connected_components(G))
        sizes = [get_number_of_edges(G, c) for c in components]

        # Get a null distribution for max component sizes
        null = Parallel(n_jobs=-1)(
            delayed(self.get_null_distribution)() for _ in range(n_reps)
        )
        null = np.array(null)

        # Compute the p-value
        tested_components = []
        for c, size in zip(components, sizes):
            p = np.sum(null >= size) / len(null)
            tested_components.append([p, c])

        return tested_components
