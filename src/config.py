"""Configuration loading — YAML file + CLI argument merging."""

from pathlib import Path
from typing import Any, Optional
import logging
import yaml

logger = logging.getLogger(__name__)

DEFAULT_CONFIG = {
    "reference_gene": "Rp49",
    "control_group": "DMSO",
    "replicate_separator": "#",
    "skiprows": 19,
    "use_cols": ["Target", "Sample", "Cq"],
    "figure": {
        "width": 8,
        "height": 6,
        "dpi": 150,
        "palette": "coolwarm",
        "bar_alpha": 0.8,
        "strip_alpha": 0.6,
        "strip_size": 6,
    },
}


def load_config(config_path: Optional[str | Path] = None) -> dict[str, Any]:
    """Load configuration from a YAML file, falling back to defaults.

    Args:
        config_path: Path to a YAML config file. If None or the file is missing,
                     the built-in defaults are returned.

    Returns:
        A dict of configuration values.
    """
    if config_path is None:
        logger.info("No config file specified, using built-in defaults.")
        return dict(DEFAULT_CONFIG)

    config_path = Path(config_path)
    if not config_path.exists():
        logger.warning(
            "Config file %s not found, using built-in defaults.", config_path
        )
        return dict(DEFAULT_CONFIG)

    with open(config_path, "r", encoding="utf-8") as fh:
        user_config = yaml.safe_load(fh) or {}

    merged = dict(DEFAULT_CONFIG)
    merged.update(user_config)

    if "figure" in user_config:
        merged["figure"] = dict(DEFAULT_CONFIG["figure"])
        merged["figure"].update(user_config["figure"])

    logger.info("Loaded config from %s", config_path)
    return merged


def merge_cli_args(
    config: dict[str, Any], args: Any
) -> dict[str, Any]:
    """Override config values with CLI arguments (non-None values only)."""
    cli_overrides = {
        "reference_gene": getattr(args, "ref_gene", None),
        "control_group": getattr(args, "control", None),
        "replicate_separator": getattr(args, "separator", None),
        "skiprows": getattr(args, "skiprows", None),
    }
    for key, value in cli_overrides.items():
        if value is not None:
            config[key] = value
            logger.debug("CLI override: %s = %s", key, value)

    for fig_key in ("width", "height", "dpi", "palette"):
        cli_val = getattr(args, f"fig_{fig_key}", None)
        if cli_val is not None:
            config.setdefault("figure", {})[fig_key] = cli_val

    return config
