# %%
from collections import defaultdict

import numpy as np
from graspologic.embed import OmnibusEmbed

from hotellings import test as hotelling

# %%
def embed(X, Y):

    stacked = np.vstack([X, Y])
    n_per_pop = len(stacked) // 2

    embedder = OmnibusEmbed(2)
    embeddings = embedder.fit_transform(stacked)

    Xhat = embeddings[:n_per_pop]
    Yhat = embeddings[n_per_pop:]

    return Xhat, Yhat


# %%
def test(X, Y, n_nodes):
    Xhat, Yhat = embed(X, Y)

    pvalues = defaultdict(list)
    for node in range(n_nodes):
        try:
            _, pvalue = hotelling(Xhat[:, node, :], Yhat[:, node, :])
            pvalues[node].append(pvalue)
        except:
            pvalues[node].append(1)

    return pvalues


# %%
