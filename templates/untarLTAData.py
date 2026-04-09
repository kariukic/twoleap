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
    help="Measurement sets",
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
    return all([os.path.isfile(tfile) for tfile in ms_tarballs])


def parse_tarballs(ms_list):
    ms_lists = [ms_list]
    return [ms for ms_list in ms_lists for ms in ms_list]


def getActualMS(msfile):
    assert os.path.isdir(msfile)
    contents = os.listdir(msfile)
    if len(contents) == 1:
        logging.info(f"inner dir named as: {contents[0]}")
        return os.path.join(msfile, contents[0])
    elif len(contents) > 1:
        assert "FIELD" in contents, "Typical MS table, FIELD, not found in MS directory"
        assert (
            "OBSERVATION" in contents
        ), "Typical MS table, OBSERVATION, not found in MS directory"
        assert (
            "SPECTRAL_WINDOW" in contents
        ), "Typical MS table, SPECTRAL_WINDOW, not found in MS directory"
        assert (
            "ANTENNA" in contents
        ), "Typical MS table, ANTENNA, not found in MS directory"

        return msfile
    else:
        logging.error(f"Empty MS directory: {msfile}")
        raise ValueError(f"Empty MS directory: {msfile}")


def main(args):
    ms_tarballs = parse_tarballs(args.ms_tarballs)
    assert all_are_files(ms_tarballs)

    for tarball in ms_tarballs:
        msname = tarball.replace(".MS", f"_{args.label}.MS")
        extract(tarball, msname)

        actualMS = getActualMS(msname)
        if actualMS != msname:
            logging.info(f"Renaming {actualMS} to {msname}")
            # os.rename(actualMS, msname)
            import shutil

            inner_dir = actualMS
            outer_dir = msname

            # Check if inner directory exists
            if not os.path.exists(inner_dir):
                raise FileNotFoundError(f"Inner directory not found: {inner_dir}")

            # Move all contents from inner_dir to outer_dir
            for item in os.listdir(inner_dir):
                src = os.path.join(inner_dir, item)
                dst = os.path.join(outer_dir, item)

                # Handle case where destination already exists (optional: overwrite or skip)
                if os.path.exists(dst):
                    print(f"Warning: Destination already exists, skipping: {dst}")
                    continue  # or use `shutil.move(src, dst)` to overwrite

                shutil.move(src, dst)  # Moves files/dirs while preserving metadata

            # Remove the now-empty inner directory
            try:
                os.rmdir(inner_dir)  # Will only work if the directory is empty
                print(f"Successfully removed empty directory: {inner_dir}")
            except OSError as e:
                print(f"Error: Could not remove {inner_dir} (may not be empty): {e}")

        # os.system(f"mv {msname}/{tarball}/* {msname}")
        # os.system(f"rm -r {msname}/{tarball}")
        logging.info(f"Extracted {tarball} ----> {os.getcwd()}/{msname}")


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
