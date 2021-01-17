import numpy as np
from graspologic.simulations import er_np, sbm
from scipy.stats import truncnorm


def generate_truncnorm_sbms(
    m, block_1, block_2, mean_1, mean_2, var_1, var_2, a=-1, b=1
):
    """
    Function for generating two populations of undirected, weighted SBMs.
    The weight function is truncated normal such that all values are in [-1, 1].
    Population 1 is sampled with `mean_1` and `variance_1` for all blocks. Population
    2 is sampled with `mean_2` and `variance_2` for community 1, and with `mean_1`
    and `variance_1` for all other communities. This function is used for
    Dos and Don'ts experiments 2 and 4.
    Parameters
    ----------
    m : int
        Number of samples per population
    block_1, block_2 : int
        Number of vertices in community 1 and 2, respectively.
    mean_1, mean_2 : float
        Means of truncated normal for community 1 and 2, respectively.
    var_1, var_2 : float
        Variances of truncated normal for community 1 and 2, respectively.
    Returns
    -------
    pop1, pop2 : 3d-array with shape (m, n, n)
        Sampled undirected, binary graphs.
    labels : 1d-array with shape (n,)
        True community assignments.
    """
    # Parameters for er and sbm functions
    total_n = block_1 + block_2
    n = [block_1, block_2]
    p = [[1, 1], [1, 1]]
    sd_1 = np.sqrt(var_1)
    sd_2 = np.sqrt(var_2)

    # deal with clip values
    a_1 = (a - mean_1) / sd_1
    b_1 = (b - mean_1) / sd_1
    a_2 = (a - mean_2) / sd_2
    b_2 = (b - mean_2) / sd_2

    pop_1 = []
    pop_2 = []
    for _ in range(m):
        # seeds are needed for joblib and scipy random functions
        # numpy random is not affected by joblib
        seeds = np.random.randint(0, 2147483647, size=4)

        wt_func = [[truncnorm.rvs, truncnorm.rvs], [truncnorm.rvs, truncnorm.rvs]]
        wt_args_1 = (dict(a=a_1, b=b_1, loc=mean_1, scale=sd_1, random_state=seeds[0]),)
        wt_args_2 = [
            [
                dict(a=a_2, b=b_2, loc=mean_2, scale=sd_2, random_state=seeds[1]),
                dict(a=a_1, b=b_1, loc=mean_1, scale=sd_1, random_state=seeds[2]),
            ],
            [
                dict(a=a_1, b=b_1, loc=mean_1, scale=sd_1, random_state=seeds[2]),
                dict(a=a_1, b=b_1, loc=mean_1, scale=sd_1, random_state=seeds[3]),
            ],
        ]

        pop_1.append(
            er_np(total_n, 1.0, directed=False, wt=truncnorm.rvs, wtargs=wt_args_1)
        )
        pop_2.append(sbm(n, p, directed=False, wt=wt_func, wtargs=wt_args_2))

    labels = np.array([0] * block_1 + [1] * block_2)

    return np.array(pop_1), np.array(pop_2), labels


def compute_pr_at_k(k, true_labels, test_statistics=None, pvalues=None):
    """
    Computes precision and recall at various k. Can compute based on p-values
    or test-statistics, but not both.
    Parameters
    ----------
    test_statistics : 2d-array like, shape (n, n)
        Test-statistics obtained from some test.
    pvalues : 2d-array like, shape (n, n)
        P-values obtained from some test.
    k : int or array-like
        Values @k to compute precision and recall for. If list, compute P/R for
        each value in list. Otherwise, it computes P/R from range(1, k+1).
    true_labels : 1d-array with shape (n,)
        True community assignments.
    Returns
    -------
    precisions, recalls : array-like
        Computed precisions and recalls.
    """
    if (test_statistics is not None) and (pvalues is not None):
        raise ValueError("You cannot supply both `test_statistics` and `pvalues`.")

    if test_statistics is not None:
        res = test_statistics
        reverse_sorting = True
    else:
        res = pvalues
        reverse_sorting = False

    label_matrix = np.zeros((len(true_labels), len(true_labels)))
    c1 = (true_labels == 0).sum()
    label_matrix[:c1, :c1] = 1

    triu_idx = np.triu_indices_from(res, k=1)
    labels_vec = label_matrix[triu_idx]
    res_vec = res[triu_idx]

    idx = np.argsort(res_vec)
    if reverse_sorting:
        idx = idx[::-1]
    sorted_labels = labels_vec[idx]

    if isinstance(k, int):
        ks = range(1, k + 1)
    else:
        ks = k

    precisions = [sorted_labels[:k].mean() for k in ks]
    recalls = [sorted_labels[:k].sum() / sorted_labels.sum() for k in ks]

    return precisions, recalls
