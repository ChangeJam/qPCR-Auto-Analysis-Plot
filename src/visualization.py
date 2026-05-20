"""Visualization — bar + strip plots with significance brackets."""

from __future__ import annotations
from pathlib import Path
import logging
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


def plot_qpcr_results(
    df: "pd.DataFrame",
    target_gene: str,
    stats_results: dict[str, dict[str, float]],
    control_group: str,
    figure_config: dict | None = None,
    output_path: str | Path | None = None,
) -> plt.Figure:
    """Generate a bar+strip plot with significance annotations for one target gene.

    Args:
        df: Result DataFrame with 'Group' and 'Value' columns.
        target_gene: Name of the target gene (used in the y-axis label).
        stats_results: Dict from run_all_comparisons().
        control_group: Label of the control group.
        figure_config: Dict of figure settings (width, height, palette, etc.).
        output_path: If set, save the figure to this path instead of showing it.

    Returns:
        The matplotlib Figure object.
    """
    if figure_config is None:
        figure_config = {}

    width = figure_config.get("width", 8)
    height = figure_config.get("height", 6)
    dpi = figure_config.get("dpi", 150)
    palette = figure_config.get("palette", "coolwarm")
    bar_alpha = figure_config.get("bar_alpha", 0.8)
    strip_alpha = figure_config.get("strip_alpha", 0.6)
    strip_size = figure_config.get("strip_size", 6)

    fig, ax = plt.subplots(figsize=(width, height))
    sns.set_style("ticks")

    # Enforce a consistent group order: control first, then treatments
    groups = df["Group"].unique()
    ordered = [control_group] + [g for g in groups if g != control_group]

    sns.barplot(
        x="Group",
        y="Value",
        data=df,
        order=ordered,
        hue="Group",
        palette=palette,
        legend=False,
        width=0.5,
        capsize=0.1,
        err_kws={"linewidth": 1.5},
        errorbar=("ci", 68),
        alpha=bar_alpha,
        ax=ax,
    )

    sns.stripplot(
        x="Group",
        y="Value",
        data=df,
        order=ordered,
        color="black",
        alpha=strip_alpha,
        size=strip_size,
        jitter=True,
        ax=ax,
    )

    # Draw significance brackets
    y_max = df["Value"].max()
    h = y_max * 0.05
    treatment_groups = [g for g in ordered if g != control_group]
    n_treatments = len(treatment_groups)

    for i, group in enumerate(treatment_groups):
        stat = stats_results.get(group)
        if stat is None:
            continue

        x_control = ordered.index(control_group)
        x_treatment = ordered.index(group)
        bracket_y = y_max + (5 + i * 3) * h
        text_y = bracket_y + h

        ax.plot(
            [x_control, x_treatment],
            [bracket_y, bracket_y],
            lw=1.5,
            c="k",
        )
        ax.text(
            (x_control + x_treatment) * 0.5,
            text_y,
            stat["asterisks"],
            ha="center",
            va="bottom",
            color="k",
            fontsize=11,
        )

    ax.set_xlabel("")
    ax.set_ylabel(f"Relative quantity of {target_gene} mRNA (Normalized)")
    ax.set_ylim(0, y_max * (1.5 + n_treatments * 0.15))

    sns.despine()
    fig.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        logger.info("Figure saved to %s", output_path)

    return fig
