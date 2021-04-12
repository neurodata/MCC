# %%
from collections import defaultdict

from hotellings import test as hotelling

# %%
def embed(X, Y, node):
    Xhat = X[:, node, :]
    Yhat = Y[:, node, :]
    return Xhat, Yhat


# %%
def test(X, Y, n_nodes):

    pvalues = defaultdict(list)

    for node in range(n_nodes):
        Xhat, Yhat = embed(X, Y, node)
        try:
            _, pvalue = hotelling(Xhat, Yhat)
            pvalues[node].append(pvalue)
        except:
            pvalues[node].append(1)

    return pvalues
