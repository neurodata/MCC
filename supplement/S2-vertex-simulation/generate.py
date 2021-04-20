# %%
from graspologic.simulations import rdpg
import numpy as np


# %% Utility functions
def rotate(x, angle):
    rads = (angle * np.pi) / 180
    rotated = np.array(
        [
            np.cos(rads) * x[0] - np.sin(rads) * x[1],
            np.sin(rads) * x[0] + np.cos(rads) * x[1],
        ]
    )
    return rotated


def stretch(x, effect_size):
    return x * np.sqrt((1 + effect_size))


def find_angle(x, y, target):
    angle = np.arccos(target / (np.linalg.norm(x) * np.linalg.norm(y))) / np.pi * 180
    return angle


def get_latent_positions(X1, X2, X3, block_size):
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
    return X, Y


# %%
def generate_graphs_1(
    p, effect_size, block_size, num_graphs, initial_angle=60.0, second_angle=80.0
):
    """
    Change magnitude, keep angle same
    Initial angle of 60 means half of p is off diagonal
    """
    assert (1.0 + effect_size) <= (1 / p)
    X1 = np.array([p, p])
    X2 = rotate(X1, initial_angle)
    X3 = stretch(X2, effect_size)

    X, Y = get_latent_positions(X1, X2, X3, block_size)

    G1 = np.array([rdpg(X) for _ in range(num_graphs)])
    G2 = np.array([rdpg(Y) for _ in range(num_graphs)])

    return G1, G2


def generate_graphs_4(
    p, effect_size, block_size, num_graphs, initial_angle=60.0, second_angle=120.0
):
    """
    Initial angle of 60 means half of p is off diagonal
    """
    assert (1.0 + effect_size) <= (1 / p)
    X1 = np.array([p, p])
    X2 = rotate(X1, initial_angle)
    X3 = rotate(X1, second_angle)

    X, Y = get_latent_positions(X1, X2, X3, block_size)

    G1 = np.array([rdpg(X) for _ in range(num_graphs)])
    G2 = np.array([rdpg(Y) for _ in range(num_graphs)])

    return G1, G2
