#! /usr/bin/env python3

# get the list from dawn
# extract the subbands
# go to the given databands and obtain files with those subbands
# if untarred already distribute them to the required nodes --> modify migrate data script to accet a list of mses as input
# if not untarred, untare thm first, then move them to the required nodes --> import the untaarring function from the untarring file


import os
import re
import logging
import subprocess
import tarfile
from glob import glob
from argparse import ArgumentParser

logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.DEBUG)

parser = ArgumentParser(
    description="Select subbnads belonging to a given redshift bin, unpack and move them"
)

parser.add_argument(
    "-z",
    "--zbin",
    type=int,
    help="EoR redshift bin",
    dest="zbin",
    default="",
)

parser.add_argument(
    "-p",
    "--ms_pattern",
    type=str,
    help="Measurement sets or file with list of MS",
    dest="ms_pattern",
    default="",
)

parser.add_argument(
    "-n",
    "--nodes",
    nargs="+",
    help="Measurement sets",
    dest="nodes",
    default="",
)

parser.add_argument(
    "-x",
    "--to_nodes",
    nargs="+",
    help="Measurement sets",
    dest="to_nodes",
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

parser.add_argument("-y", "--directory_B", help="to this directory", required=False)

parser.add_argument(
    "-f",
    "--nfiles_per_worker",
    type=int,
    help="number of files o distribute to each worker node. the master node gets the remainder",
    required=False,
)

parser.add_argument(
    "-m",
    "--mode",
    type=str,
    help="Either move, copy or symlink. Default is move",
    choices=["move", "copy", "symlink"],
    dest="mode",
    default="mode",
)

parser.add_argument(
    "-t",
    "--dry_run",
    help="test first before actually moving the data",
    dest="dry_run",
    action="store_true",
)


def extract(tfile, msname):
    tar = tarfile.open(tfile)
    tar.extractall(path=msname)
    tar.close()


def allAreExtracted(mses):
    return all([os.path.isdir(m) for m in mses])


def allAreFiles(mses):
    return all([os.path.isfile(m) for m in mses])


def getMSlist(fullpath):
    return [f"{fullpath}/{m}" for m in os.listdir(fullpath)]


def readTxt2List(txt_file):
    with open(txt_file, "r") as txt:
        return txt.readlines()


def getMsSubband(mspath):
    return re.findall(r"SB\d+", mspath)[0][2:]


def missingSubbands(subbands):
    min_sb = int(sorted(subbands)[0])
    max_sb = int(sorted(subbands)[-1])

    missing = [m for m in subbands if int(m) not in range(min_sb, max_sb + 1)]

    return missing


def parseNodes(nodes):
    if ".." in nodes[0]:
        min_node = int(nodes[0].split("..")[0])
        max_node = int(nodes[0].split("..")[1])
        return list(range(min_node, max_node + 1))
    elif "," in nodes[0]:
        return [int(n) for n in nodes[0].split(",")]
    else:
        return nodes


def chunks(lst: list, n: int):
    """
    Yield successive n-sized chunks from a list.

    Parameters:
    - lst (list): The list from which to yield chunks.
    - n (int): The size of each chunk.

    Returns:
    - generator: A generator that yields successive n-sized chunks of the input list.

    Example:
    >>> lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    >>> for chunk in chunks(lst, 3):
    ...     print(chunk)
    [1, 2, 3]
    [4, 5, 6]
    [7, 8, 9]
    [10]
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def distributeMSets(
    files,
    to_nodes,
    directory_B,
    nfiles_per_worker,
    dry_run=True,
    mode="move",
    label=None,
):
    files_per_node = list(chunks(files, nfiles_per_worker))

    logging.info(to_nodes)
    for n, nodeB in enumerate(to_nodes):
        # nodeB = int(nodeB) #.strip("[").strip("]").strip(","))
        datadirB = os.path.join(directory_B)
        files_chunk = files_per_node[n]

        mkdir_cmd = f"ssh node{nodeB} 'mkdir -p {datadirB}'"
        logging.info(f"running {mkdir_cmd}")
        if not dry_run:
            subprocess.run(mkdir_cmd, shell=True)
            logging.info(f"/net/node{nodeB}/{datadirB}")

            assert os.path.isdir(f"/net/node{nodeB}/{datadirB}")

        for fyl in files_chunk:

            # if os.path.isfile(fyl):
                # try:
                #     assert label, "label must be given"
                #     msname = fyl.replace(".MS", f"_{label}.MS")
                #     if not dry_run:
                #         extract(fyl, msname)
                #         # os.system(f"mv {msname}/{fyl}/* {msname}")
                #     logging.info(f"Moving {msname}/{fyl}/* ----> {msname}")
                #     logging.info(f"Extracted {fyl} ----> {msname}")

                # except Exception as e:
                #     logging.error(f"Could not extract {fyl}: {e}")

                # fyl = msname

            symlink_name = (
                os.path.basename(fyl).replace(".MS", f"_{label}.MS")
                if label
                else os.path.basename(fyl)
            )

            if mode == "copy":
                cmd = f"ssh node{nodeB} 'scp -r {fyl} {datadirB}/{symlink_name}'"

            elif mode == "symlink":
                cmd = f"ssh node{nodeB} 'ln -s {fyl} {datadirB}/{symlink_name}'"

            elif mode == "move":
                cmd = f"ssh node{nodeB} 'mv -f {fyl} {datadirB}/{symlink_name}'"

            logging.info(f"running {mode}: {cmd}")
            if not dry_run:
                if not os.path.exists(f"/net/node{nodeB}/{datadirB}/{symlink_name}"):
                    subprocess.run(cmd, shell=True)
                else:
                    logging.info(
                        f"File /net/node{nodeB}/{datadirB}/{symlink_name} exists already. Skipping."
                    )


def main(args):
    subbands_per_redshift_bin = {"1": [], "2": list(range(98, 165)), "3": list(range(33, 101)), "4": list(range(14, 278))}

    redshift_bin_subbands = subbands_per_redshift_bin[str(args.zbin)]
    logging.info(
        f"{len(redshift_bin_subbands)} Redshift bin sbands: {redshift_bin_subbands}"
    )

    nodes = parseNodes(args.nodes)
    logging.debug(f"Nodes: {nodes}")

    mses = []
    for node in nodes:
        glob_pattern = f"/net/node{node}/{args.ms_pattern}"
        logging.debug(f"glob: {glob_pattern}")
        mses += glob(glob_pattern)

    logging.info(f"All MSes found: {len(mses)}")

    #Make sure the mses found to belong the the redshift bin com out sorted
    mses = sorted(mses, key=lambda x: int(getMsSubband(x)))
    red_mses = [ms for ms in mses if int(getMsSubband(ms)) in redshift_bin_subbands]
    assert len(red_mses) > 0, "No msfiles found for that redshift bin."

    logging.debug(f" All redshift bin MSes: {red_mses}")
    logging.info(f" {len(red_mses)} redshift bin MSes found")

    if allAreExtracted(red_mses):
        logging.info(f"All {len(red_mses)} MSes are EXTRACTED")

    elif allAreFiles(red_mses):
        logging.info(f"All {len(red_mses)} MSes are UNEXTRACTED")

    missing_subbands = missingSubbands([getMsSubband(m) for m in red_mses])
    if missing_subbands:
        logging.info(f"MISSING subbands: {missing_subbands}")
    else:
        logging.info("0 missing subbands")

    if args.directory_B:
        assert (
            args.to_nodes
        ), "Provided directory_B but no to_nodes. to_nodes is needed to redistributed data"
        assert (
            args.nfiles_per_worker
        ), "nfiles_per_worker is needed to redistributed data"

        to_nodes = parseNodes(args.to_nodes)
        distributeMSets(
            red_mses,
            to_nodes,
            args.directory_B,
            args.nfiles_per_worker,
            dry_run=args.dry_run,
            mode=args.mode,
            label=args.label,
        )


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)