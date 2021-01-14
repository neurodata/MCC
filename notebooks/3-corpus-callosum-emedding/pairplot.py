"""
Implement pairplots with probability ellipses on the lower diagonal.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.stats.distributions import chi2


def _get_parameters(x, y):
    mu = np.array([x.mean(), y.mean()])
    cov = np.cov(x, y)
    return mu, cov


def _get_eigen(cov):
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    return eigvals, eigvecs


def _get_ellipse_parameters(cov, ci):

    eigvals, eigvecs = _get_eigen(cov)

    # Calculate angle of displacement from x-axis
    theta = np.arctan2(*eigvecs[:, 0][::-1])
    theta = np.degrees(theta)

    # Calculate scaling factor based on probability
    dof = len(eigvals)
    alpha = 1 - (1 - ci) / 2
    scale = chi2.ppf(alpha, dof)
    width, height = 2 * np.sqrt(scale * eigvals)

    return width, height, theta


def make_ellipse(x, y, ci, **kwargs):
    mu, cov = _get_parameters(x, y)
    width, height, theta = _get_ellipse_parameters(cov, ci)
    ellipse = Ellipse(mu, width, height, theta, alpha=0.5, **kwargs)
    return ellipse


def draw_probability_ellipse(x, y, ci, **kwargs):
    ellipse = make_ellipse(x, y, ci, **kwargs)
    ax = plt.gca()
    ax.add_artist(ellipse)

    plot_kws = dict(alpha=0.75, s=15, linewidth=0, marker="o")
    ax.scatter(x, y, **plot_kws, **kwargs)


def munge_embedding(embedding, labels):
    col_names = [f"Dimension {i}" for i in range(1, embedding.shape[1] + 1)]
    df = pd.DataFrame(embedding, columns=col_names)
    df["Strain"] = labels
    return df


def _pairgrid(embedding, labels):
    df = munge_embedding(embedding, labels)

    palette = ["#e7298a", "#1b9e77", "#d95f02", "#7570b3"]
    plot_kws = dict(alpha=0.75, s=15, linewidth=0, marker="o")

    with sns.plotting_context("paper", font_scale=1):
        g = sns.PairGrid(df, hue="Strain", palette=palette, height=1.5)
        g = g.map_upper(sns.scatterplot, **plot_kws)
        g = g.map_diag(sns.kdeplot, lw=1, shade=True)
        g.set(xticks=[], yticks=[])

    return g


def ellipse_pairgrid(embedding, labels, ci, **kwargs):
    g = _pairgrid(embedding, labels)
    kwargs["ci"] = ci
    g = g.map_lower(draw_probability_ellipse, **kwargs)
    return g


def kde_pairgrid(embedding, labels):
    g = _pairgrid(embedding, labels)
    g = g.map_lower(sns.kdeplot)
    return g
