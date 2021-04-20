# %%
from itertools import product

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

from generate import generate_graphs_1, generate_graphs_4
from omni import test as test_omni
from mdmr import test as test_mdmr
from nxstats import test as test_nxstats

# %% Experiment function
def experiment(
    sample_size, block_1, block_2, p, effect_size, second_angle, generate_func, reps=100
):

    m_per_pop = sample_size // 2
    block_size = np.array([block_1, block_2])
    n_nodes = block_1 + block_2
    labels = np.array([0] * block_1 + [1] * block_2)

    def worker(i):
        X, Y = generate_func(
            p, effect_size, block_size, m_per_pop, second_angle=second_angle
        )
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
n_nodes = 50
block_2s = np.array([0, 5, 10, 15, 20, 25])
block_1s = 50 - block_2s
block_sizes = list(zip(block_1s, block_2s))
effect_size = 0.5
second_angle = 65.0
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
        second_angle=second_angle,
        generate_func=generate_func,
    )
    for (block_1, block_2), generate_func in product(block_sizes, generate_funcs)
]

# %% Run the experiment
for arg in tqdm(args):
    experiment(**arg)

# %%
