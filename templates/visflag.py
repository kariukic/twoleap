#!/usr/bin/env python3

# from nenucal import delayflag

# ms_file = "/net/node300/data/users/lofareor/chege/3C196/L192832/zbin3/data/fullband_dd_smooth_data_6hrs.MS"
# config = "/home/codex/chege/projects/3C196/notebooks_v2/flags/default_vis_flagger_settings.toml"
# pdir = (
#     "/net/node300/data/users/lofareor/chege/3C196/L192832/zbin3/data/visflagger_plots/"
# )

# delayflag.apply_vis_filter(ms_file, config, plot_dir=pdir, dry_run=False)


"""Delay flagging application for radio astronomy data."""

import click
import logging
from pathlib import Path
from nenucal import delayflag

# from typing import Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def validate_path(ctx, param, value):
    """Click callback for path validation."""
    if value is None:
        return None
    path = Path(value)
    if not path.exists():
        raise click.BadParameter(f"Path {value} does not exist")
    return path


@click.command()
@click.argument(
    "ms_file",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    required=True,
)
@click.option(
    "--config",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to flagging configuration TOML file",
)
@click.option(
    "--plot-dir",
    type=click.Path(file_okay=False, writable=True),
    required=True,
    help="Directory to save diagnostic plots",
)
@click.option(
    "--dry-run", is_flag=True, help="Run without applying flags (simulation only)"
)
def main(ms_file: str, config: str, plot_dir: str, dry_run: bool):
    """
    Apply delay-based flagging to measurement set data.

    Example:
        python flagger.py /path/to/data.MS \
            --config /path/to/config.toml \
            --plot-dir /output/plots/ \
            --dry-run
    """
    try:
        # Ensure plot directory exists
        plot_path = Path(plot_dir)
        plot_path.mkdir(parents=True, exist_ok=True)

        logger.info(f"Flagging: {ms_file}")
        logger.info(f"Using config: {config}")
        logger.info(f"Saving plots to: {plot_dir}")

        delayflag.apply_vis_filter(
            ms_file=ms_file,
            config_file=config,
            plot_dir=str(plot_path),
            dry_run=dry_run,
        )

        logger.info("pspipe visflagger completed successfully")
    except Exception as e:
        logger.error(f"visflagger failed: {str(e)}")
        raise click.ClickException(str(e))


if __name__ == "__main__":
    main()
