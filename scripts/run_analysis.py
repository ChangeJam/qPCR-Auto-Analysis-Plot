#!/usr/bin/env python3
"""qPCR Auto Analysis — Command-line entry point.

Example:
    python scripts/run_analysis.py \\
        --input data/admin_2025-12-10.csv \\
        --config config.yaml \\
        --ref-gene Rp49 \\
        --control DMSO
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running from the repo root without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import logging
import numpy as np

from src.config import load_config, merge_cli_args
from src.io import read_qpcr_csv, write_result_csv
from src.preprocessing import (
    group_by_target,
    align_target_to_ref,
    parse_sample_groups,
    get_control_mask,
)
def _parse_label_mapping(raw: str | None) -> dict[str, str]:
    """Parse a 'key=value,...' string into a dict."""
    if not raw:
        return {}
    mapping = {}
    for pair in raw.split(","):
        pair = pair.strip()
        if "=" not in pair:
            continue
        k, v = pair.split("=", 1)
        mapping[k.strip()] = v.strip()
    return mapping


from src.analysis import (
    calculate_delta_cq,
    calculate_fold_change,
    normalize_to_control,
    build_result_dataframe,
)
from src.statistics import run_all_comparisons
from src.visualization import plot_qpcr_results


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(levelname)-8s %(message)s",
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Automated qPCR data analysis and plotting."
    )
    p.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the Bio-Rad CFX Maestro CSV export.",
    )
    p.add_argument(
        "-c", "--config",
        default="config.yaml",
        help="Path to YAML configuration file (default: config.yaml).",
    )
    p.add_argument(
        "--ref-gene",
        help="Name of the reference (housekeeping) gene in the Target column.",
    )
    p.add_argument(
        "--control",
        help="Label of the control group (e.g. DMSO).",
    )
    p.add_argument(
        "--separator",
        help="Character separating group name from replicate number in Sample names.",
    )
    p.add_argument(
        "--skiprows",
        type=int,
        help="Number of metadata rows to skip before the column header.",
    )
    p.add_argument(
        "-o", "--output-dir",
        default="output",
        help="Directory for result CSV files (default: output/).",
    )
    p.add_argument(
        "-f", "--figures-dir",
        default="figures",
        help="Directory for output figures (default: figures/).",
    )
    p.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip generating plots.",
    )
    p.add_argument(
        "--group-labels",
        help="Comma-separated raw=display name mappings for group labels. "
             "Example: PFOA001=PFOA 0.01mg/L,PFOA1=PFOA 1mg/L",
    )
    p.add_argument(
        "--list-targets",
        action="store_true",
        help="Only list the Target genes found in the file and exit.",
    )
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug-level logging.",
    )
    # Figure overrides
    p.add_argument("--fig-width", type=float, help="Figure width in inches.")
    p.add_argument("--fig-height", type=float, help="Figure height in inches.")
    p.add_argument("--fig-dpi", type=int, help="Figure DPI.")
    p.add_argument("--fig-palette", help="Seaborn colour palette name.")
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.verbose)
    logger = logging.getLogger("run_analysis")

    # --- Configuration ---
    config = load_config(args.config)
    config = merge_cli_args(config, args)

    ref_gene = config["reference_gene"]
    control_group = config["control_group"]
    separator = config["replicate_separator"]
    skiprows = config["skiprows"]

    # Parse group label overrides: CLI > config
    raw_label_overrides = dict(config.get("group_labels", {}))
    cli_overrides = _parse_label_mapping(args.group_labels)
    raw_label_overrides.update(cli_overrides)

    # Map control group through labels if needed
    display_control = raw_label_overrides.get(control_group, control_group)

    # --- Read data ---
    df = read_qpcr_csv(args.input, skiprows=skiprows)

    # --- List targets mode ---
    if args.list_targets:
        targets = df["Target"].unique()
        print("Targets found in the file:")
        for t in targets:
            marker = " (reference)" if t == ref_gene else ""
            print(f"  - {t}{marker}")
        return

    # --- Validate reference gene ---
    if ref_gene not in df["Target"].values:
        available = list(df["Target"].unique())
        logger.error(
            "Reference gene '%s' not found in Target column. Available: %s",
            ref_gene, available,
        )
        sys.exit(1)

    # --- Group by target ---
    targets = group_by_target(df)

    ref_df = targets[ref_gene]
    target_genes = {k: v for k, v in targets.items() if k != ref_gene}

    if not target_genes:
        logger.error("No target genes found besides the reference gene.")
        sys.exit(1)

    logger.info(
        "Reference gene: %s | Target genes: %s",
        ref_gene, list(target_genes.keys()),
    )

    output_dir = Path(args.output_dir)
    figures_dir = Path(args.figures_dir)

    for gene_name, gene_df in target_genes.items():
        logger.info("--- Processing %s ---", gene_name)

        # Align target and reference DataFrames by Sample (positional)
        try:
            target_cq, ref_cq, aligned_samples = align_target_to_ref(
                gene_df, ref_df
            )
        except ValueError as exc:
            logger.error("Alignment failed for %s: %s", gene_name, exc)
            continue

        # Parse group labels from the aligned sample names
        aligned_groups = parse_sample_groups(aligned_samples, separator)
        group_list = aligned_groups.map(
            lambda g: raw_label_overrides.get(g, g)
        ).tolist()
        control_mask = get_control_mask(aligned_groups, control_group)

        # --- Analysis pipeline ---
        delta_cq = calculate_delta_cq(target_cq, ref_cq)
        fold_change = calculate_fold_change(delta_cq)
        normalized = normalize_to_control(fold_change, control_mask)

        result_df = build_result_dataframe(normalized, group_list)

        # --- Export CSV ---
        csv_path = output_dir / f"{gene_name}_qPCR_result.csv"
        write_result_csv(result_df, csv_path)

        # --- Statistics ---
        stats_results = run_all_comparisons(result_df, display_control)

        # --- Plot ---
        if not args.no_plot:
            fig_path = figures_dir / f"{gene_name}_qPCR_plot.png"
            plot_qpcr_results(
                result_df,
                target_gene=gene_name,
                stats_results=stats_results,
                control_group=display_control,
                figure_config=config.get("figure"),
                output_path=fig_path,
            )
            import matplotlib.pyplot as plt
            plt.close("all")

    logger.info("All done. Output: %s  Figures: %s", output_dir.resolve(), figures_dir.resolve())


if __name__ == "__main__":
    main()
