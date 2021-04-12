# %%
from itertools import product

import graspologic as gp
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

from omni import test as test_omni
from mdmr import test as test_mdmr
from nxstats import test as test_nxstats

# %% Generating functions
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


# %% Experiment function
def experiment(sample_size, block_1, block_2, p, effect_size, generate_func, reps=100):

    m_per_pop = sample_size // 2
    block_size = np.array([block_1, block_2])
    n_nodes = block_1 + block_2
    labels = np.array([0] * block_1 + [1] * block_2)

    def worker(i):
        X, Y = generate_func(p, effect_size, block_size, m_per_pop)
        pvalues_omni = pd.DataFrame.from_dict(test_omni(X, Y, n_nodes))
        pvalues_mdmr = pd.DataFrame.from_dict(test_mdmr(X, Y, n_nodes))
        pvalues_nx = pd.DataFrame.from_dict(test_nxstats(X, Y, n_nodes))
        df = make_df(
            pvalues_omni, pvalues_mdmr, pvalues_nx, labels, block_2, generate_func
        )
        df.to_csv(f"results/tmp/{block_2}-{generate_func.__name__}-{i}.csv")

    Parallel(-1)(delayed(worker)(i) for i in range(reps))


def make_df(pvalues_omni, pvalues_mdmr, pvalues_nx, labels, block_2, generate_func):
    df = pd.concat([pvalues_omni, pvalues_mdmr, pvalues_nx])
    df = df.T
    df.columns = ["omni", "mdmr", "stats"]
    df["node"] = df.index
    df["label"] = labels
    df["n_signal_vertices"] = int(block_2)
    df["generating_func"] = generate_func.__name__
    return df


# %% Simulation parameters
n_nodes = 70
block_2s = np.array([0, 10, 20, 30, 40])
block_1s = 70 - block_2s
block_sizes = list(zip(block_1s, block_2s))
effect_size = 1
sample_size = 200
p = 0.25

generate_funcs = [generate_graphs_1, generate_graphs_4]

args = [
    dict(
        sample_size=sample_size,
        block_1=block_1,
        block_2=block_2,
        p=p,
        effect_size=effect_size,
        generate_func=generate_func,
    )
    for (block_1, block_2), generate_func in product(block_sizes, generate_funcs)
]

# %% Run the experiment
for arg in tqdm(args):
    experiment(**arg)
