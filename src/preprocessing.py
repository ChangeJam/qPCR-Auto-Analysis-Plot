"""Data preprocessing — grouping, sorting, and sample-name parsing."""

import logging
import pandas as pd

logger = logging.getLogger(__name__)


def group_by_target(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Split the DataFrame into a dict keyed by unique Target (gene) names.

    Each value is sorted by Sample and Cq values are extracted as a Series.
    """
    targets = {}
    for name in df["Target"].unique():
        subset = df[df["Target"] == name].sort_values("Sample")
        targets[name] = subset
    logger.info("Found %d targets: %s", len(targets), list(targets.keys()))
    return targets


def parse_sample_groups(
    sample_series: pd.Series, separator: str = "#"
) -> pd.Series:
    """Extract group labels from sample names by stripping the replicate suffix.

    e.g. "DMSO#1" -> "DMSO", "PFOA001#2" -> "PFOA001"

    Args:
        sample_series: Series of sample names.
        separator: Character that separates group name from replicate number.

    Returns:
        Series of group labels (same index as input).
    """
    groups = sample_series.apply(
        lambda s: s.rsplit(separator, 1)[0] if separator in str(s) else str(s)
    )
    unique_groups = groups.unique()
    logger.info(
        "Parsed %d unique groups from sample names: %s",
        len(unique_groups),
        list(unique_groups),
    )
    return groups


def get_control_mask(
    groups: pd.Series, control_label: str
) -> "np.ndarray":
    """Return a boolean mask for rows belonging to the control group.

    Raises:
        ValueError: If the control label is not found.
    """
    import numpy as np

    mask = (groups == control_label).values
    if not mask.any():
        available = list(groups.unique())
        raise ValueError(
            f"Control group '{control_label}' not found. Available: {available}"
        )
    logger.info(
        "Control group '%s': %d of %d rows.",
        control_label, mask.sum(), len(mask),
    )
    return mask


def align_target_to_ref(
    target_df: "pd.DataFrame",
    ref_df: "pd.DataFrame",
) -> tuple["np.ndarray", "np.ndarray", "pd.Series"]:
    """Align target and reference DataFrames by Sample for positional subtraction.

    Both DataFrames are sorted by Sample.  This function validates that their
    Sample columns match row-for-row and returns the aligned Cq arrays plus
    the sample groups (derived from sample names).

    Returns:
        (target_cq, ref_cq, groups) — numpy arrays for Cq values and a Series
        of group labels, all aligned to the same row order.
    """
    import numpy as np

    target_samples = target_df["Sample"].values
    ref_samples = ref_df["Sample"].values

    if len(target_samples) != len(ref_samples):
        raise ValueError(
            f"Row count mismatch: target has {len(target_samples)} rows "
            f"but reference has {len(ref_samples)} rows."
        )

    if not np.array_equal(target_samples, ref_samples):
        # Find mismatches for a helpful error message
        diffs = np.where(target_samples != ref_samples)[0]
        n_diffs = len(diffs)
        raise ValueError(
            f"Sample alignment failed: {n_diffs} row(s) differ between "
            f"target and reference. First mismatch at row {diffs[0]}: "
            f"'{target_samples[diffs[0]]}' vs '{ref_samples[diffs[0]]}'."
        )

    target_cq = target_df["Cq"].values.astype(float)
    ref_cq = ref_df["Cq"].values.astype(float)

    logger.info("Aligned %d samples between target and reference.", len(target_cq))
    return target_cq, ref_cq, target_df["Sample"]
