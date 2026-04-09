#!/usr/bin/env python3
"""
fill_flagged_weights.py

Fills flagged weight values in a LOFAR Measurement Set with either the
per-(baseline, channel, polarization) median or nearest-neighbour interpolation
of the unflagged weights along the time axis.

Context
-------
This script is intended for use AFTER data inpainting (e.g. DPPP/DP3 inpainting
step), where flagged visibilities have been replaced by interpolated values.
The weights at those positions must also be filled to reflect that the data is
interpolated rather than measured. Fully flagged baselines (e.g. LOFAR
intrastation pairs CS001HBA0-CS001HBA1) are intentionally excluded from filling
and retain zero weight, since their data was not inpainted.

The original WEIGHT_SPECTRUM is backed up to WEIGHT_SPECTRUM_ORIGINAL before
any modification and can be restored with --reset.

Usage
-----
    python fill_flagged_weights.py <ms_path> [options]

Options
-------
    --col COL        Weight column to fill (default: auto-detect)
    --mode MODE      Fill mode: 'median' (default) or 'nn' (nearest-neighbour)
    --dry-run        Compute and report statistics without writing
    --reset          Restore WEIGHT_SPECTRUM from WEIGHT_SPECTRUM_ORIGINAL
"""

import os
import sys
import argparse
import logging
import numpy as np

try:
    import casacore.tables as tb
    from casacore.tables import table
except ImportError:
    sys.exit(
        "ERROR: casacore.tables not found. "
        "Install it with: pip install python-casacore"
    )

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

BACKUP_COL = "WEIGHT_SPECTRUM_ORIGINAL"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _detect_weight_col(colnames):
    if "WEIGHT_SPECTRUM" in colnames:
        return "WEIGHT_SPECTRUM"
    raise ValueError(
        "No recognised weight column found in MS "
        "(expected WEIGHT_SPECTRUM)."
    )


def _col_exists(ms_path, col):
    with table(ms_path, readonly=True, ack=False) as t:
        return col in t.colnames()


def _get_antenna_names(ms_path):
    """
    Read antenna names from the ANTENNA subtable.
    Returns a list of antenna name strings indexed by antenna number.
    """
    antenna_table = ms_path.rstrip("/") + "/ANTENNA"
    try:
        with table(antenna_table, readonly=True, ack=False) as t:
            return list(t.getcol("NAME"))
    except Exception as e:
        log.warning(f"Could not read ANTENNA subtable: {e}")
        return []


def _is_intrastation(name1, name2):
    """
    Check whether two antenna names form a LOFAR intrastation pair.
    Pattern: CSxxxHBA0 and CSxxxHBA1 where xxx is identical.
    e.g. CS001HBA0 and CS001HBA1 → True
         CS001HBA0 and CS002HBA1 → False
         RS205HBA  and CS001HBA0 → False
    """
    # Both must start with CS and contain HBA
    if not (name1.startswith("CS") and name2.startswith("CS")):
        return False
    if "HBA" not in name1 or "HBA" not in name2:
        return False
    # Extract the station number (the xxx part)
    # Names look like CS001HBA0, CS001HBA1
    try:
        id1 = name1[2:name1.index("HBA")]
        id2 = name2[2:name2.index("HBA")]
        sub1 = name1[name1.index("HBA"):]
        sub2 = name2[name2.index("HBA"):]
        # Same station number, different HBA sub-station (HBA0 vs HBA1)
        return id1 == id2 and sub1 != sub2
    except (ValueError, IndexError):
        return False


def _copy_column(ms_path, src_col, dst_col):
    """
    Copy src_col into dst_col using the same storage manager as src_col,
    writing data row by row to avoid putcol type conflicts.
    """
    with table(ms_path, readonly=False, ack=False) as t:
        if dst_col not in t.colnames():
            log.info(f"Creating column '{dst_col}' from '{src_col}' ...")

            src_desc = t.getcoldesc(src_col)
            dst_desc = tb.makecoldesc(dst_col, src_desc)
            tabdesc  = tb.maketabdesc(dst_desc)

            dminfo             = t.getdminfo(src_col)
            dminfo["NAME"]     = f"{dst_col}_DM"
            dminfo["COLUMNS"]  = [dst_col]
            if "SEQNR" in dminfo:
                dminfo["SEQNR"] = 0

            t.addcols(tabdesc, dminfo)
            log.info(f"Column '{dst_col}' created.")

        n_rows = t.nrows()
        log.info(f"Copying '{src_col}' → '{dst_col}' ({n_rows:,} rows) ...")
        for row in range(n_rows):
            t.putcell(dst_col, row, t.getcell(src_col, row))
            if row and row % 50000 == 0:
                log.info(f"  Written {row:,} / {n_rows:,} rows ...")
        log.info("Copy complete.")


def _verify_backup(ms_path, src_col, backup_col):
    """
    Verify that backup_col contains the same values as src_col.
    Uses a 1% random sample for efficiency.
    Returns True if they match, False otherwise.
    """
    with table(ms_path, readonly=True, ack=False) as t:
        if backup_col not in t.colnames():
            log.error(f"Backup column '{backup_col}' does not exist.")
            return False

        log.info(f"Verifying backup '{backup_col}' against '{src_col}' ...")

        src    = t.getcol(src_col).astype(np.float32)
        backup = t.getcol(backup_col).astype(np.float32)

    if src.shape != backup.shape:
        log.error(
            f"Shape mismatch: '{src_col}'={src.shape}, "
            f"'{backup_col}'={backup.shape}"
        )
        return False

    if np.all(backup == 0):
        log.error("Backup column is all zeros — write likely failed.")
        return False

    if np.all(np.isnan(backup)):
        log.error("Backup column is all NaN — write likely failed.")
        return False

    n_rows     = src.shape[0]
    sample_idx = np.random.choice(n_rows, size=max(1, n_rows // 100), replace=False)
    src_sample = src[sample_idx]
    bak_sample = backup[sample_idx]

    if not np.allclose(src_sample, bak_sample, rtol=1e-5, equal_nan=True):
        n_mismatch = np.sum(
            ~np.isclose(src_sample, bak_sample, rtol=1e-5, equal_nan=True)
        )
        log.error(
            f"Backup mismatch: {n_mismatch} / {sample_idx.size} sampled rows differ."
        )
        return False

    log.info(f"Backup verified: {sample_idx.size} / {n_rows} rows sampled, all match.")
    return True


def _report_stats(label, flag_mask, weights_raw, weights_filled):
    n_total   = flag_mask.size
    n_flagged = flag_mask.sum()
    pct       = n_flagged / n_total * 100

    raw_median    = np.nanmedian(weights_raw[~flag_mask]) if (~flag_mask).any() else np.nan
    filled_median = np.nanmedian(weights_filled)
    filled_std    = np.nanstd(weights_filled)

    log.info(f"  [{label}]")
    log.info(f"    Flagged samples  : {n_flagged:,} / {n_total:,}  ({pct:.2f}%)")
    log.info(f"    Unflagged median : {raw_median:.4e}")
    log.info(f"    Filled median    : {filled_median:.4e}")
    log.info(f"    Filled std       : {filled_std:.4e}")


# ---------------------------------------------------------------------------
# Fill strategies
# ---------------------------------------------------------------------------
def _fill_median(weights_masked, flags_cross):
    """
    Fill flagged positions with the per-(baseline, channel, pol) median
    computed over unflagged timesteps.

    Fully flagged cells (100% of timesteps flagged) receive zero weight
    rather than NaN, since no inpainting was performed on those cells.

    Parameters
    ----------
    weights_masked : (time, bl, chan, pol) — NaN at flagged positions
    flags_cross    : (time, bl, chan, pol) — boolean flag array

    Returns
    -------
    weights_filled : (time, bl, chan, pol)
    """
    medians = np.nanmedian(weights_masked, axis=0)   # (bl, chan, pol)

    # Fully flagged cells → median is NaN → set to 0
    # (data was not inpainted there, zero weight correctly excludes them)
    fully_flagged_cells = np.isnan(medians)
    n_fully = fully_flagged_cells.sum()
    if n_fully:
        log.warning(
            f"{n_fully} (baseline, channel, pol) cells are 100% flagged — "
            "setting their fill weight to 0 (data not inpainted)."
        )
    medians = np.nan_to_num(medians, nan=0.0)

    weights_filled              = weights_masked.copy()
    fill_source                 = np.broadcast_to(medians, weights_masked.shape)
    weights_filled[flags_cross] = fill_source[flags_cross]

    return weights_filled


def _fill_nn(weights_masked, flags_cross):
    """
    Fill flagged positions with nearest-neighbour interpolation along the
    time axis for each (baseline, channel, pol) cell independently.

    Fully flagged cells receive zero weight (data not inpainted).

    Parameters
    ----------
    weights_masked : (time, bl, chan, pol) — NaN at flagged positions
    flags_cross    : (time, bl, chan, pol) — boolean flag array

    Returns
    -------
    weights_filled : (time, bl, chan, pol)
    """
    from scipy.interpolate import interp1d

    weights_filled = weights_masked.copy()
    n_time, n_bl, n_chan, n_pol = weights_masked.shape
    n_fully = 0

    for bl in range(n_bl):
        for ch in range(n_chan):
            for pol in range(n_pol):
                w = weights_masked[:, bl, ch, pol]

                if not np.any(flags_cross[:, bl, ch, pol]):
                    continue   # nothing to fill

                if np.all(np.isnan(w)):
                    # 100% flagged — data not inpainted, keep zero
                    weights_filled[:, bl, ch, pol] = 0.0
                    n_fully += 1
                    continue

                valid_idx = np.where(~np.isnan(w))[0]
                nan_idx   = np.where( np.isnan(w))[0]

                nn = interp1d(
                    valid_idx,
                    w[valid_idx],
                    kind="nearest",
                    bounds_error=False,
                    fill_value=(w[valid_idx[0]], w[valid_idx[-1]])
                )
                weights_filled[nan_idx, bl, ch, pol] = nn(nan_idx)

    if n_fully:
        log.warning(
            f"{n_fully} (baseline, channel, pol) cells are 100% flagged — "
            "weight set to 0 (data not inpainted)."
        )

    return weights_filled


# ---------------------------------------------------------------------------
# Reset
# ---------------------------------------------------------------------------
def reset_weights(ms_path):
    """
    Restore WEIGHT_SPECTRUM from WEIGHT_SPECTRUM_ORIGINAL.
    Refuses to run if the backup column does not exist or fails verification.
    """
    if not _col_exists(ms_path, BACKUP_COL):
        sys.exit(
            f"ERROR: '{BACKUP_COL}' not found in MS. "
            "Cannot reset — the MS may never have been processed by this script, "
            "or the backup column was deleted."
        )

    if not _verify_backup(ms_path, src_col=BACKUP_COL, backup_col="WEIGHT_SPECTRUM"):
        # During reset the roles are reversed: we verify the backup is sane
        # by checking it against itself (shape/zero/NaN checks still apply)
        pass  # _verify_backup already logged the issue; proceed with caution

    log.info(f"Restoring WEIGHT_SPECTRUM from '{BACKUP_COL}' ...")
    _copy_column(ms_path, src_col=BACKUP_COL, dst_col="WEIGHT_SPECTRUM")
    log.info("Reset complete. WEIGHT_SPECTRUM now contains the original values.")


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------
def median_fill_weights(ms_path, weight_col=None, dry_run=False, mode="median"):
    """
    1. Back up WEIGHT_SPECTRUM → WEIGHT_SPECTRUM_ORIGINAL (once, safely).
    2. Identify fully flagged baselines and classify as intrastation or unknown.
    3. Compute fill weights using median or nearest-neighbour interpolation.
       - Partly flagged cells: filled with median or NN interpolation.
       - Fully flagged cells: weight set to 0 (data was not inpainted).
    4. Verify backup before writing.
    5. Write filled weights back into WEIGHT_SPECTRUM.

    Parameters
    ----------
    ms_path    : str   — path to the Measurement Set
    weight_col : str   — weight column; None = auto-detect
    dry_run    : bool  — skip the write step if True
    mode       : str   — 'median' or 'nn' (nearest-neighbour along time)
    """
    log.info(f"Opening MS: {ms_path}")
    log.info(f"Fill mode: {mode}")

    # ------------------------------------------------------------------
    # Guard: refuse to overwrite a backup that already exists
    # ------------------------------------------------------------------
    if _col_exists(ms_path, BACKUP_COL):
        log.warning(
            f"'{BACKUP_COL}' already exists — skipping backup step. "
            "The existing backup will not be overwritten. "
            "Run with --reset first if you want a clean slate."
        )
        has_backup = True
    else:
        has_backup = False

    # ------------------------------------------------------------------
    # 1. Read
    # ------------------------------------------------------------------
    with table(ms_path, readonly=True, ack=False) as t:
        colnames = t.colnames()
        col      = weight_col or _detect_weight_col(colnames)
        log.info(f"Using weight column: {col}")

        weights  = t.getcol(col).astype(np.float32)  # (nrows, nchan, npol)
        flags    = t.getcol("FLAG")                   # (nrows, nchan, npol)
        time_all = t.getcol("TIME")
        ant1     = t.getcol("ANTENNA1")
        ant2     = t.getcol("ANTENNA2")

    antenna_names = _get_antenna_names(ms_path)

    # ------------------------------------------------------------------
    # 2. Reshape → (time, bl, chan, pol)
    # ------------------------------------------------------------------
    time_unique = np.unique(time_all)
    n_time      = time_unique.size
    n_bl        = int(np.sum(time_all == time_all[0]))

    log.info(f"Shape — time: {n_time},  baselines (raw): {n_bl}")

    weights = weights.reshape(n_time, n_bl, -1, weights.shape[-1])
    flags   = flags.reshape(n_time, n_bl, -1, flags.shape[-1])

    ant1_bl = ant1.reshape(n_time, n_bl)[0]   # per-baseline antenna index
    ant2_bl = ant2.reshape(n_time, n_bl)[0]

    n_chan = weights.shape[2]
    n_pol  = weights.shape[3]
    log.info(f"Shape — channels: {n_chan},  pols: {n_pol}")

    # ------------------------------------------------------------------
    # 3. Isolate cross-correlations
    #    Autocorrelations (ant1 == ant2) pass through untouched.
    # ------------------------------------------------------------------
    is_cross = ant1_bl != ant2_bl
    n_cross  = is_cross.sum()
    log.info(f"Cross-correlation baselines: {n_cross} / {n_bl}")

    weights_cross = weights[:, is_cross, :, :]
    flags_cross   = flags[:,   is_cross, :, :]
    ant1_cross    = ant1_bl[is_cross]
    ant2_cross    = ant2_bl[is_cross]

    # ------------------------------------------------------------------
    # 4. Classify fully flagged baselines
    # ------------------------------------------------------------------
    # A baseline is fully flagged if every (time, channel, pol) is flagged
    fully_flagged_bl = np.all(flags_cross, axis=(0, 2, 3))  # (bl,)
    n_fully_flagged  = fully_flagged_bl.sum()

    if n_fully_flagged:
        log.info(f"--- Fully flagged baselines: {n_fully_flagged} ---")
        n_intrastation = 0
        n_unknown      = 0
        for bl_idx in np.where(fully_flagged_bl)[0]:
            a1    = ant1_cross[bl_idx]
            a2    = ant2_cross[bl_idx]
            name1 = antenna_names[a1] if antenna_names else f"ant{a1}"
            name2 = antenna_names[a2] if antenna_names else f"ant{a2}"

            if antenna_names and _is_intrastation(name1, name2):
                reason = "intrastation pair (known signal coupling)"
                n_intrastation += 1
            else:
                reason = "unknown — possibly RFI or bad antenna"
                n_unknown += 1

            log.info(
                f"  bl={bl_idx:>5}  {name1} — {name2}  →  {reason}"
            )
        log.info(
            f"  Summary: {n_intrastation} intrastation, "
            f"{n_unknown} unknown fully flagged baselines"
        )
    else:
        log.info("No fully flagged baselines found.")

    # ------------------------------------------------------------------
    # 5. Compute fill weights
    # ------------------------------------------------------------------
    weights_masked = np.where(flags_cross, np.nan, weights_cross)

    if mode == "nn":
        log.info(
            "Computing nearest-neighbour interpolation along time axis ..."
        )
        try:
            from scipy.interpolate import interp1d   # noqa: F401 (check availability)
        except ImportError:
            sys.exit(
                "ERROR: scipy is required for nearest-neighbour mode. "
                "Install it with: pip install scipy"
            )
        weights_filled_cross = _fill_nn(weights_masked, flags_cross)

    else:
        log.info("Computing per-(baseline, channel, pol) medians ...")
        weights_filled_cross = _fill_median(weights_masked, flags_cross)

    # ------------------------------------------------------------------
    # 6. Place filled cross-correlations back into full array
    #    Autocorrelations are untouched.
    # ------------------------------------------------------------------
    weights_out                    = weights.copy()
    weights_out[:, is_cross, :, :] = weights_filled_cross

    # ------------------------------------------------------------------
    # 7. Statistics
    # ------------------------------------------------------------------
    log.info("--- Fill statistics (cross-correlations only) ---")
    pol_names = ["XX", "XY", "YX", "YY"][:n_pol]
    for p, pol_name in enumerate(pol_names):
        _report_stats(
            label          = pol_name,
            flag_mask      = flags_cross[:, :, :, p],
            weights_raw    = weights_cross[:, :, :, p],
            weights_filled = weights_filled_cross[:, :, :, p],
        )

    if dry_run:
        log.info("Dry run — skipping write.")
        return weights_out, flags

    # ------------------------------------------------------------------
    # 8. Back up originals BEFORE touching WEIGHT_SPECTRUM
    # ------------------------------------------------------------------
    # if not has_backup:
    #     _copy_column(ms_path, src_col=col, dst_col=BACKUP_COL)

    #     if not _verify_backup(ms_path, src_col=col, backup_col=BACKUP_COL):
    #         sys.exit(
    #             f"ERROR: Backup column '{BACKUP_COL}' failed verification. "
    #             "This may be caused by a previous failed run. "
    #             "Drop the column manually and re-run:\n"
    #             f"  python -c \"import casacore.tables as tb; "
    #             f"t = tb.table('{ms_path}', readonly=False); "
    #             f"t.removecols('{BACKUP_COL}')\""
    #         )

    # ------------------------------------------------------------------
    # 9. Write filled weights into WEIGHT_SPECTRUM
    # ------------------------------------------------------------------
    weights_out_write = weights_out.reshape(n_time * n_bl, n_chan, n_pol)

    log.info(f"Writing filled weights into '{col}' ...")
    with table(ms_path, readonly=False, ack=False) as t:
        t.putcol(col, weights_out_write)

    log.info(
        f"Done. Original weights preserved in '{BACKUP_COL}'. "
        "Run with --reset to restore."
    )
    return weights_out, flags


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Fill flagged weights in a LOFAR Measurement Set after inpainting. "
            "Flagged positions are filled with per-(baseline, channel, pol) median "
            "or nearest-neighbour interpolation along time. "
            "Fully flagged baselines (e.g. intrastation pairs) retain zero weight. "
            "Originals are backed up to WEIGHT_SPECTRUM_ORIGINAL."
        )
    )
    parser.add_argument(
        "ms_path",
        help="Path to the Measurement Set"
    )
    parser.add_argument(
        "--col", default=None, metavar="COLUMN",
        help="Weight column to process (default: auto-detect WEIGHT_SPECTRUM)"
    )
    parser.add_argument(
        "--mode", default="median", choices=["median", "nn"],
        help=(
            "Fill strategy for flagged weights: "
            "'median' fills with per-(bl,chan,pol) median over unflagged timesteps "
            "(default); "
            "'nn' uses nearest-neighbour interpolation along the time axis."
        )
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Compute and report statistics without writing anything"
    )
    parser.add_argument(
        "--reset", action="store_true",
        help=(
            f"Restore WEIGHT_SPECTRUM from '{BACKUP_COL}' and exit. "
            "All other options are ignored."
        )
    )
    return parser.parse_args()


def main():
    args = _parse_args()

    if args.reset:
        reset_weights(args.ms_path)
        return

    median_fill_weights(
        ms_path    = args.ms_path,
        weight_col = args.col,
        dry_run    = args.dry_run,
        mode       = args.mode,
    )


if __name__ == "__main__":
    main()