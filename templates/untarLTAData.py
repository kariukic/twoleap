#! /usr/bin/env python3
import os
import tarfile
import logging
from argparse import ArgumentParser
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.DEBUG)

parser = ArgumentParser(description="unpack LTA Measurement set tarballs")

parser.add_argument(
    "-i",
    "--ms_tarballs",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="ms_tarballs",
    default="",
)

parser.add_argument(
    "-l",
    "--label",
    help="",
    required=False,
    type=str,
    dest="label",
)


def extract(tfile, msname):
    tar = tarfile.open(tfile)
    tar.extractall(path=msname)
    tar.close()

    return


def all_are_files(ms_tarballs):
    return all([os.path.isfile(tarfile) for tarfile in ms_tarballs])


def parse_tarballs(ms_list):
    ms_lists = [ms_list]
    return [ms for ms_list in ms_lists for ms in ms_list]


def main(args):
    ms_tarballs = parse_tarballs(args.ms_tarballs)
    assert all_are_files(ms_tarballs)

    for tarball in ms_tarballs:
        msname = tarball.replace('.MS', f'_{args.label}.MS')
        extract(tarball, msname)
        logging.info(f"Extracted {tarball} ----> {os.getcwd()}/{msname}")


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
