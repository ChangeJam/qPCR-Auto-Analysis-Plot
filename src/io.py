"""File I/O operations for qPCR data."""

from pathlib import Path
import logging
import pandas as pd

logger = logging.getLogger(__name__)


class QPCRDataError(Exception):
    """Raised when the input qPCR data fails validation."""


def read_qpcr_csv(
    filepath: str | Path,
    skiprows: int = 19,
    usecols: list[str] | None = None,
) -> pd.DataFrame:
    """Read a Bio-Rad CFX Maestro CSV export and return a cleaned DataFrame.

    Args:
        filepath: Path to the CSV file.
        skiprows: Number of metadata rows to skip before the column header.
        usecols: List of column names to keep. Defaults to Target, Sample, Cq.

    Returns:
        DataFrame with the requested columns.

    Raises:
        FileNotFoundError: If the file does not exist.
        QPCRDataError: If required columns are missing or the data is empty.
    """
    if usecols is None:
        usecols = ["Target", "Sample", "Cq"]

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    try:
        df = pd.read_csv(filepath, skiprows=skiprows)
    except UnicodeDecodeError:
        df = pd.read_csv(filepath, skiprows=skiprows, encoding="latin-1")

    missing = [c for c in usecols if c not in df.columns]
    if missing:
        available = list(df.columns)
        raise QPCRDataError(
            f"Required column(s) {missing} not found. Available columns: {available}"
        )

    df = df[usecols].copy()

    if df.empty:
        raise QPCRDataError("Input file contains no data rows after skipping metadata.")

    # Drop rows where Cq is NaN (wells with no amplification)
    n_before = len(df)
    df = df.dropna(subset=["Cq"])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        logger.info("Dropped %d row(s) with missing Cq values.", n_dropped)

    df["Cq"] = pd.to_numeric(df["Cq"], errors="coerce")
    df = df.dropna(subset=["Cq"])
    if df.empty:
        raise QPCRDataError("No valid Cq values remain after cleaning.")

    logger.info("Loaded %d rows from %s", len(df), filepath.name)
    return df


def write_result_csv(
    df: pd.DataFrame, output_path: str | Path
) -> None:
    """Write a result DataFrame to CSV.

    Args:
        df: DataFrame with 'Group' and 'Value' columns.
        output_path: Destination path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.info("Wrote result to %s", output_path)
