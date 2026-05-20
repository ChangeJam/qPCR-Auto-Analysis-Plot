"""Statistical tests for qPCR analysis."""

from __future__ import annotations
import logging
from scipy import stats

logger = logging.getLogger(__name__)


def run_ttest(
    control_values: "list[float] | np.ndarray",
    treatment_values: "list[float] | np.ndarray",
) -> tuple[float, float]:
    """Run an independent two-sample t-test (two-tailed).

    Returns:
        (t_statistic, p_value)
    """
    import numpy as np

    control_values = np.asarray(control_values, dtype=float)
    treatment_values = np.asarray(treatment_values, dtype=float)
    t_stat, p_val = stats.ttest_ind(control_values, treatment_values)
    logger.debug("t-test: stat=%.4f, p=%.6f", t_stat, p_val)
    return float(t_stat), float(p_val)


def pvalue_to_asterisks(p: float) -> str:
    """Convert a p-value to significance asterisks (Mann–Whitney style).

    p <= 0.0001  ->  ****
    p <= 0.001   ->  ***
    p <= 0.01    ->  **
    p <= 0.05    ->  *
    p >  0.05    ->  ns
    """
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    return "ns"


def run_all_comparisons(
    df: "pd.DataFrame",
    control_group: str,
) -> dict[str, dict[str, float]]:
    """Run t-tests for every treatment group against the control.

    Args:
        df: Result DataFrame with 'Group' and 'Value' columns.
        control_group: Label of the control group.

    Returns:
        Dict mapping treatment_group -> {"t_stat": ..., "p_value": ..., "asterisks": ...}
    """
    control_vals = df.loc[df["Group"] == control_group, "Value"]
    if len(control_vals) == 0:
        raise ValueError(f"Control group '{control_group}' not found in DataFrame.")

    treatment_groups = [g for g in df["Group"].unique() if g != control_group]
    results = {}
    for group in treatment_groups:
        treat_vals = df.loc[df["Group"] == group, "Value"]
        t_stat, p_val = run_ttest(control_vals, treat_vals)
        asterisks = pvalue_to_asterisks(p_val)
        results[group] = {
            "t_stat": t_stat,
            "p_value": p_val,
            "asterisks": asterisks,
        }
        logger.info(
            "%s vs %s: t=%.4f, p=%.6f (%s)",
            group, control_group, t_stat, p_val, asterisks,
        )
    return results
