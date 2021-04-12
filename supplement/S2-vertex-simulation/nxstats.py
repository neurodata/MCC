# %%
from collections import defaultdict

import numpy as np
import pandas as pd
import networkx as nx

from hotellings import test as hotelling

# %% Calculate network statistics for a single graph
def summarize(graph):

    # Define the network statistics to compute
    functions = [
        nx.clustering,
        nx.betweenness_centrality,
        nx.closeness_centrality,
        nx.triangles,
    ]

    # Get a dictionary of network statistics for each node in the graph
    G = nx.from_numpy_array(graph)
    stats = defaultdict(list)
    for func in functions:
        for node, scalar in func(G).items():
            stats[node].append(scalar)

    # Organize statistics in a pandas df
    df = pd.DataFrame.from_dict(stats, orient="index")
    df["node"] = df.index
    df.columns = [func.__name__ for func in functions] + ["node"]

    return df


# %% Embed all graphs
def embed(graphs):
    frames = [summarize(graph) for graph in graphs]
    df = pd.concat(frames)
    df.reset_index(level=0, inplace=True, drop=True)
    return df


# %%
def test(X, Y, n_nodes):

    Xhat = embed(X)
    Yhat = embed(Y)

    pvalues = defaultdict(list)
    for node in range(n_nodes):
        Xhat_node = Xhat.query(f"node == {node}")
        Yhat_node = Yhat.query(f"node == {node}")
        try:
            _, pvalue = hotelling(Xhat_node, Yhat_node)
            pvalues[node].append(pvalue)
        except:
            pvalues[node].append(1)

    return pvalues
