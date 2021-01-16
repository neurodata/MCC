# %%
from collections import defaultdict
from itertools import product

import graspologic as gp
import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import confusion_matrix, f1_score
from tqdm import tqdm

from hotellings import Hotelling


# %%
def estimate_embeddings(X, Y, method):
    """
    Parameters
    ----------
    method : str
        Must be {'mase', 'omni'}
    """
    stacked = np.vstack([X, Y])

    if method == "mase":
        embedder = gp.embed.MultipleASE(2)
        embeddings = embedder.fit_transform(stacked)
    elif method == "omni":
        embedder = gp.embed.OmnibusEmbed(2)
        embeddings = embedder.fit_transform(stacked).mean(axis=0)
    else:
        assert ValueError("Invalid embedding method")

    return embeddings


# %%
def rotate(x, angle):
    rads = (angle * np.pi) / 180
    rotated = np.array(
        [
            np.cos(rads) * x[0] - np.sin(rads) * x[1],
            np.sin(rads) * x[0] + np.cos(rads) * x[1],
        ]
    )

    return rotated


def find_angle(x, y, target):
    angle = np.arccos(target / (np.linalg.norm(x) * np.linalg.norm(y))) / np.pi * 180
    return angle


def generate_graphs_1(
    p, effect_size, block_size, num_graphs, initial_angle=60.0, **kwargs
):
    """
    Change magnitude, keep angle same
    Initial angle of 60 means half of p is off diagonal
    """
    assert (1.0 + effect_size) <= (1 / p)
    X1 = np.array([p, p])
    X2 = rotate(X1, initial_angle)
    X3 = X2 * np.sqrt(1 + effect_size)

    X = np.vstack(
        [
            np.repeat(np.array([X1]), block_size[0], axis=0),
            np.repeat(np.array([X2]), block_size[1], axis=0),
        ]
    )
    Y = np.vstack(
        [
            np.repeat(np.array([X1]), block_size[0], axis=0),
            np.repeat(np.array([X3]), block_size[1], axis=0),
        ]
    )

    P1 = X @ X.T
    P2 = Y @ Y.T

    G1 = np.array([gp.simulations.sample_edges(P1) for _ in range(num_graphs)])
    G2 = np.array([gp.simulations.sample_edges(P2) for _ in range(num_graphs)])

    return G1, G2


def generate_graphs_4(
    p, effect_size, block_size, num_graphs, initial_angle=60.0, **kwargs
):
    """
    Initial angle of 60 means half of p is off diagonal
    """
    assert (1.0 + effect_size) <= (1 / p)
    X1 = np.array([p, p])
    X2 = rotate(X1, initial_angle)
    off_diag = X1 @ X2
    X3 = X2 * np.sqrt(1 + effect_size)
    X3 = rotate(X3, find_angle(X1, X3, off_diag) - initial_angle)

    X = np.vstack(
        [
            np.repeat(np.array([X1]), block_size[0], axis=0),
            np.repeat(np.array([X2]), block_size[1], axis=0),
        ]
    )
    Y = np.vstack(
        [
            np.repeat(np.array([X1]), block_size[0], axis=0),
            np.repeat(np.array([X3]), block_size[1], axis=0),
        ]
    )

    P1 = X @ X.T
    P2 = Y @ Y.T

    G1 = np.array([gp.simulations.sample_edges(P1) for _ in range(num_graphs)])
    G2 = np.array([gp.simulations.sample_edges(P2) for _ in range(num_graphs)])

    return G1, G2


# %%
def calculate_statistics(graph):
    G = nx.from_numpy_array(graph)

    stats = defaultdict(list)

    for func in (
        nx.clustering,
        nx.betweenness_centrality,
        nx.closeness_centrality,
        nx.triangles,
    ):
        for node, scalar in func(G).items():
            stats[node].append(scalar)

    df = pd.DataFrame.from_dict(stats, orient="index")
    df["node"] = df.index

    df.columns = [
        func.__name__
        for func in [
            nx.clustering,
            nx.betweenness_centrality,
            nx.closeness_centrality,
            nx.triangles,
        ]
    ] + ["node"]
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]

    return df


def make_df(graphs):
    frames = [calculate_statistics(graph) for graph in graphs]
    df = pd.concat(frames)
    df.reset_index(level=0, inplace=True, drop=True)
    return df


# %%
def compute_pr_at_k(pvals, k, true_labels):
    """
    1 is considered the significant nodes
    """
    idx = np.argsort(pvals)
    sorted_labels = true_labels[idx]

    precision = sorted_labels[:k].mean()
    # recall = sorted_labels[:k].sum() / sorted_labels.sum()

    return precision


def experiment(m, block_1, block_2, p, effect_size, generate_func, reps=100):
    m_per_pop = m // 2
    block_size = np.array([block_1, block_2])

    n = block_1 + block_2
    labels = np.array([0] * block_1 + [1] * block_2)

    def worker():
        pvals = np.zeros((3, n))

        X, Y = generate_func(p, effect_size, block_size, m_per_pop)
        X_stats = make_df(X)
        Y_stats = make_df(Y)

        for node in range(n):
            try:
                X_node = X_stats.query(f"node == {node}")
                Y_node = Y_stats.query(f"node == {node}")
                test = Hotelling().test(X_node, Y_node)
                # test = KSample('Dcorr').test(X_node, Y_node, reps=500)
                if np.isnan(test[1]):
                    pvals[0, node] = 0
                else:
                    pvals[0, node] = test[1]
            except:
                pvals[0, node] = 1

        for j, method in enumerate(["omni", "mase"], start=1):
            embeddings = estimate_embeddings(X, Y, method, 2, sample_space=True)
            Xhat = embeddings[:m_per_pop]
            Yhat = embeddings[m_per_pop:]
            for node in range(n):
                try:
                    test = Hotelling().test(Xhat[:, node, :], Yhat[:, node, :])
                    pvals[j, node] = test[1]
                except:
                    pvals[j, node] = 1

        p_at_k = np.array(
            [compute_pr_at_k(pvals[i], block_2, labels) for i in range(3)]
        )

        pvals_binary = (pvals <= 0.05) * 1

        f1 = np.array([f1_score(labels, pvals_binary[i]) for i in range(3)])

        conf = np.array(
            [confusion_matrix(labels, pvals_binary[i]).ravel() for i in range(3)]
        )
        # tn, fp, fn, tp
        fp = conf[:, 1] / n
        tp = conf[:, 3] / n

        return pvals, p_at_k, f1, fp, tp

    res = Parallel(-1)(delayed(worker)() for _ in range(reps))
    precs = np.array([res[i][1] for i in range(reps)]).mean(axis=0)
    f1s = np.array([res[i][2] for i in range(reps)]).mean(axis=0)
    fps = np.array([res[i][3] for i in range(reps)]).mean(axis=0)
    tps = np.array([res[i][4] for i in range(reps)]).mean(axis=0)

    pvals = np.array([res[i][0] for i in range(reps)])
    pvals = np.array(pvals)
    pvals = (pvals < 0.05).mean(axis=0)
    avg_pval_1 = pvals[:, :block_1].mean(axis=1)
    avg_pval_2 = pvals[:, block_1:].mean(axis=1)

    to_append = [
        m,
        block_1,
        block_2,
        effect_size,
        generate_func.__name__,
        *avg_pval_1,
        *avg_pval_2,
        *precs,
        *f1s,
        *fps,
        *tps,
    ]
    return to_append


# %%
m = 200
# block_1s = np.arange(9, 76, 3) #51
block_2s = np.arange(1, 70)  # 51
block_1s = 70 - block_2s
block_sizes = list(zip(block_1s, block_2s))
effect_size = 1
generate_funcs = [generate_graphs_1, generate_graphs_4]
p = 0.25

args = [
    dict(
        p=p,
        m=m,
        block_1=block_1,
        block_2=block_2,
        effect_size=effect_size,
        generate_func=generate_func,
    )
    for (block_1, block_2), generate_func in product(block_sizes, generate_funcs)
]


# %%
res = []
for arg in tqdm(args):
    res.append(experiment(**arg))


# %%
columns = [
    "m",
    "block_1",
    "block_2",
    "effect_size",
    "model",
    "scalar_same",
    "omni_same",
    "mase_same",
    "scalar_diff",
    "omni_diff",
    "mase_diff",
    "scalar_precision",
    "omni_precision",
    "mase_precision",
    "scalar_f1",
    "omni_f1",
    "mase_f1",
    "scalar_fp",
    "omni_fp",
    "mase_fp",
    "scalar_tp",
    "omni_tp",
    "mase_tp",
]

df = pd.DataFrame(res, columns=columns)
df.to_csv("../../results/vertex_simulation.csv", index=False)
