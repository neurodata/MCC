import numpy as np
from scipy.stats import f


def statistic(x, y):
    r"""
    Calulates the Hotelling :math:`T^2` test statistic.

    Parameters
    ----------
    x,y : ndarray
        Input data matrices. ``x`` and ``y`` must have the same number of
        dimensions. That is, the shapes must be ``(n, p)`` and ``(m, p)`` where
        `n` is the number of samples and `p` and `q` are the number of
        dimensions.

    Returns
    -------
    stat : float
        The computed Hotelling :math:`T^2` statistic.
    """
    # ported from Hotteling packge in R
    nx, p = x.shape
    ny = y.shape[0]

    meanx = np.mean(x, axis=0)
    meany = np.mean(y, axis=0)

    covx = np.cov(x, rowvar=False)
    covy = np.cov(y, rowvar=False)

    covs = ((nx - 1) * covx + (ny - 1) * covy) / (nx + ny - 2)
    m = (nx + ny - p - 1) / (p * (nx + ny - 2))
    if p > 1:
        inv_covs = np.linalg.pinv(covs)
        stat = m * (meanx - meany).T @ inv_covs @ (meanx - meany) * nx * ny / (nx + ny)
    else:
        inv_covs = 1 / p
        stat = m * (meanx - meany) ** 2 * inv_covs * nx * ny / (nx + ny)

    return stat


def test(x, y):
    r"""
    Calculates the Hotelling :math:`T^2` test statistic and p-value.

    Parameters
    ----------
    x,y : ndarray
        Input data matrices. ``x`` and ``y`` must have the same number of
        dimensions. That is, the shapes must be ``(n, p)`` and ``(m, p)`` where
        `n` is the number of samples and `p` and `q` are the number of
        dimensions.

    Returns
    -------
    stat : float
        The computed Hotelling :math:`T^2` statistic.
    pvalue : float
        The computed Hotelling :math:`T^2` p-value.

    Examples
    --------
    >>> import numpy as np
    >>> from hyppo.ksample import Hotelling
    >>> x = np.arange(7)
    >>> y = x
    >>> stat, pvalue = Hotelling().test(x, y)
    >>> '%.3f, %.1f' % (stat, pvalue)
    '0.000, 1.0'
    """
    stat = statistic(x, y)
    nx, p = x.shape
    ny = y.shape[0]
    pvalue = f.sf(stat, p, nx + ny - p - 1)

    return stat, pvalue
