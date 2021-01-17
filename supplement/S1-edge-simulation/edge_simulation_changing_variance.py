# %%
import warnings
from itertools import product
from functools import partial

import numpy as np
import pandas as pd
from graspologic.utils import symmetrize
from joblib import Parallel, delayed
from scipy.stats import ttest_ind, mannwhitneyu
from hyppo.ksample import KSample

from src import generate_truncnorm_sbms, compute_pr_at_k

warnings.filterwarnings("ignore")


# %%
def compute_statistic(test, pop1, pop2):
    if test.__name__ == "ttest_ind":
        test_statistics, pvals = ttest_ind(pop1, pop2, axis=0)
        np.nan_to_num(test_statistics, copy=False)
        np.nan_to_num(pvals, copy=False)
    else:  # for other tests, do by edge
        n = pop1.shape[-1]
        test_statistics = np.zeros((n, n))
        pvals = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                x_ij = pop1[:, i, j]
                y_ij = pop2[:, i, j]

                if test.__name__ == "multiscale_graphcorr":
                    tmp, pval, _ = test(x_ij, y_ij, is_twosamp=True, reps=1)
                elif test.__name__ == "test":
                    tmp, pval = KSample("Dcorr").test(x_ij, y_ij)
                else:
                    tmp, pval = test(x_ij, y_ij)

                test_statistics[i, j] = tmp
                pvals[i, j] = pval

        test_statistics = symmetrize(test_statistics, method="triu")
        pvals = symmetrize(pvals, method="triu")

    return test_statistics, pvals


def run_experiment(tests, m, block_1, block_2, mean_1, mean_2, var_1, var_2, ks, reps):
    precisions = []
    recalls = []

    for _ in range(reps):
        tmp_precisions = []
        tmp_recalls = []
        pop1, pop2, true_labels = generate_truncnorm_sbms(
            m=m,
            block_1=block_1,
            block_2=block_2,
            mean_1=mean_1,
            mean_2=mean_2,
            var_1=var_1,
            var_2=var_2,
        )

        for test in tests:
            test_statistics, pvalues = compute_statistic(test, pop1, pop2)
            if test.__name__ == "multiscale_graphcorr":
                precision, recall = compute_pr_at_k(
                    k=ks, true_labels=true_labels, test_statistics=test_statistics
                )
            else:
                precision, recall = compute_pr_at_k(
                    k=ks, true_labels=true_labels, pvalues=pvalues
                )
            tmp_precisions.append(precision)
            tmp_recalls.append(recall)

        precisions.append(tmp_precisions)
        recalls.append(tmp_recalls)

    precisions = np.array(precisions).mean(axis=0)
    recalls = np.array(recalls).mean(axis=0)

    to_append = [m, mean_1, mean_2, var_1, var_2, *precisions, *recalls]

    return to_append


# %%
tests = [ttest_ind, mannwhitneyu, KSample("Dcorr").test]

spacing = 100
block_1 = 5
block_2 = 15
mean_1 = 0
mean_2 = 0
var_1 = 0.25
var_2s = np.linspace(var_1, 3 + var_1, spacing + 1)
ms = np.linspace(0, 500, spacing + 1).astype(int)[1:]
ks = range(5, 11)
reps = 100

args = [dict(m=m, var_2=var_2) for (m, var_2) in product(ms, var_2s)]

partial_func = partial(
    run_experiment,
    tests=tests,
    block_1=block_1,
    block_2=block_2,
    mean_1=mean_1,
    mean_2=mean_2,
    var_1=var_1,
    ks=ks,
    reps=reps,
)

res = Parallel(n_jobs=-1, verbose=1)(delayed(partial_func)(**arg) for arg in args)


# %%
new_res = []

for r in res:
    constants = r[:5]
    results = [b for a in r[5:] for b in a]
    new_res.append(constants + results)


# %%
cols = [
    "m",
    "mean1",
    "mean2",
    "var_1",
    "var_2",
    *[f"{test.__name__}_precision_at_{k}" for test in tests for k in ks],
    *[f"{test.__name__}_recall_at_{k}" for test in tests for k in ks],
]
res_df = pd.DataFrame(new_res, columns=cols)
res_df.to_csv("../../results/edge_simulation_changing_variance.csv", index=False)
