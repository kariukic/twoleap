#!/usr/bin/env python3
"""
Script to split MS into frequency subbands
"""
import os
import numpy as np
import casacore.tables as pt

import logging
from argparse import ArgumentParser

logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.DEBUG)

parser = ArgumentParser(description="Split a measurement set in frequency")
parser.add_argument(
    "-m",
    "--msin",
    help="Path to input MS",
    type=str,
    required=True,
)


parser.add_argument(
    "-n",
    "--nchans_per_msout",
    type=int,
    help="Number of channels making each MS",
    required=False,
)

parser.add_argument(
    "-d",
    "--datacolumn",
    help="The datacolumn to be used",
    required=True,
    type=str,
)


def makeSubbands(
    ms,
    nchans_per_msout,
    datacolumn,
):

    spw_table = pt.table(ms + "/SPECTRAL_WINDOW")
    num_channels = spw_table.getcol(
        "NUM_CHAN"
    )  # get the total number of channel sin the ms file

    assert num_channels.shape == (1,)
    num_channels = num_channels[0]

    assert num_channels > nchans_per_msout
    assert num_channels % nchans_per_msout == 0

    startchans = np.arange(0, num_channels, nchans_per_msout)

    for s, startchan in enumerate(startchans):

        # msout = ms.replace('.MS', f'_C{s:03}.MS')
        msout = f"{ms.replace('.MS', '')}_C{s:03}.MS"

        comm = f"DP3 steps=[] msin={ms} msin.datacolumn={datacolumn} msin.startchan={startchan} msin.nchan={nchans_per_msout} msout.overwrite=true msout={msout}"

        logging.info(f"{s}: {comm}")

        os.system(comm)

    return


if __name__ == "__main__":
    args = parser.parse_args()
    assert os.path.isdir(args.msin)
    makeSubbands(args.msin, args.nchans_per_msout, args.datacolumn)
