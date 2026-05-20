"""qPCR data analysis — ΔCq, fold change, and normalization."""

import logging
import numpy as np

logger = logging.getLogger(__name__)


def calculate_delta_cq(
    target_cq: np.ndarray, ref_cq: np.ndarray
) -> np.ndarray:
    """Compute ΔCq = target_Cq - reference_Cq (element-wise subtraction).

    Args:
        target_cq: Array of Cq values for the target gene.
        ref_cq: Array of Cq values for the reference gene (same length).

    Returns:
        Array of ΔCq values.
    """
    if len(target_cq) != len(ref_cq):
        raise ValueError(
            f"Target and reference arrays must have the same length "
            f"(got {len(target_cq)} vs {len(ref_cq)})."
        )
    return target_cq - ref_cq


def calculate_fold_change(delta_cq: np.ndarray) -> np.ndarray:
    """Compute 2^(-ΔCq) for each element."""
    return np.power(2.0, -delta_cq)


def normalize_to_control(
    fold_changes: np.ndarray, control_mask: np.ndarray
) -> np.ndarray:
    """Normalize fold-change values by dividing by the control group mean.

    Args:
        fold_changes: Array of 2^(-ΔCq) values.
        control_mask: Boolean array marking control-group rows.

    Returns:
        Normalized array (same shape as fold_changes).
    """
    control_mean = fold_changes[control_mask].mean()
    if control_mean == 0:
        raise ValueError("Control group mean is zero — cannot normalize.")
    logger.info("Control group mean 2^(-ΔCq) = %.4f", control_mean)
    return fold_changes / control_mean


def build_result_dataframe(
    normalized: np.ndarray,
    group_labels: list[str],
) -> "pd.DataFrame":
    """Build the final result DataFrame with Group and Value columns."""
    import pandas as pd

    return pd.DataFrame({"Group": group_labels, "Value": normalized})
