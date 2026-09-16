import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def visualise_shift_matrix(
    SHIFT_df,
    significance_threshold=None,
    figsize=(12,8),
    cmap="RdBu_r",
    vmax=None,
    ax=None,
):
    """
    Visualise SHIFT scores as a heatmap.

    Parameters
    ----------
    SHIFT_df : pandas.DataFrame
        DataFrame containing columns:
        ct1, ct2, SHIFT

    significance_threshold : float, optional
        p-value threshold for significance. If provided, non-significant SHIFT scores will be masked in the heatmap.

    figsize : tuple, optional
        Figure size.

    cmap : str, optional
        Matplotlib colormap.

    vmax : float, optional
        Maximum colour scale value. If None, uses the
        largest absolute SHIFT score.

    ax : matplotlib.axes.Axes, optional
        Existing axes to plot into.

    Returns
    -------
    matplotlib.axes.Axes
    """

    glue = SHIFT_df.pivot(
        index="ct1",
        columns="ct2",
        values="SHIFT"
    )

    if vmax is None:
        vmax = abs(glue.values).max()

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if significance_threshold is not None:
        threshold = -np.log10(significance_threshold)
        glue = glue.where(np.abs(glue) >= threshold)

    sns.heatmap(
        glue,
        cmap=cmap,
        vmin=-vmax,
        vmax=vmax,
        square=True,
        cbar_kws={"label": "SHIFT score"},
        ax=ax,
    )

    ax.set_xlabel("Celltype 1")
    ax.set_ylabel("Celltype 2")
    plt.tight_layout()

    return ax
    