import matplotlib.pyplot as plt
import seaborn as sns


def stripplot(df, data_column):

    kwargs = {
        "alpha": 0.75,
        "edgecolor": None,
        "linewidth": 0,
        "marker": "o",
        "palette": ["#e7298a", "#1b9e77", "#d95f02", "#7570b3"],
    }

    fig, ax = plt.subplots(figsize=(2.5, 2.5), dpi=150)

    g = sns.stripplot(
        x="Strain",
        y=data_column,
        data=df,
        jitter=True,
        orient="v",
        ax=ax,
        **kwargs,
    )

    return g
