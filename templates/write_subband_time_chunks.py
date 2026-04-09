#!/usr/bin/env python3
"""
Script to split MS into frequency subbands
"""
import os
import re
from collections import defaultdict
from glob import glob
from pathlib import Path


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
    "--nmses_per_node",
    type=int,
    help="Number of subbands to distribute per node",
    required=False,
)

parser.add_argument(
    "-d",
    "--datapath",
    help="The datacolumn to be used",
    required=True,
    type=str,
)

parser.add_argument(
    "-f",
    "--from_nodes",
    help="List of nodes where the timechunks data is",
    nargs="+",
    required=True,
)

parser.add_argument(
    "-t",
    "--to_nodes",
    help="List of nodes where the subbands data is to be distributed",
    nargs="+",
    required=True,
)


def extract_T_and_C_numbers(strings):
    t_numbers = set()
    c_numbers = set()

    # Define the pattern to match 'T' followed by 3 digits and 'C' followed by 3 digits
    t_pattern = re.compile(r"T(\d{3})")
    c_pattern = re.compile(r"C(\d{3}).MS")

    # Loop through each string in the list
    for s in strings:
        # Find the 'T' number
        t_match = t_pattern.search(s)
        if t_match:
            t_numbers.add(t_match.group(1))

        # Find the 'C' number
        c_match = c_pattern.search(s)
        if c_match:
            c_numbers.add(c_match.group(1))

    return sorted(list(t_numbers)), sorted(list(c_numbers))


def extract_C_number(file_path):
    """Extract the C number from the file path."""
    c_pattern = re.compile(r"C(\d{3}).MS")
    match = c_pattern.search(file_path)
    if match:
        return match.group(1)
    return None


def group_files_by_C(files):
    """Group files by their C number."""
    grouped_files = defaultdict(list)

    for file in files:
        c_number = extract_C_number(file)
        if c_number:
            grouped_files[c_number].append(file)

    return grouped_files


def write_txt_files_to_nodes(grouped_files, datapath, nodes, max_txt_files_per_node):
    """Write the grouped files into JSON files and distribute them across nodes.

    Args:
        grouped_files (dict): Files grouped by C numbers.
        nodes (list): List of node paths where JSON files should be saved.
        max_txt_files_per_node (int): Maximum number of JSON files per node.

    Returns:
        dict: Mapping of nodes to JSON files written on them.
    """
    node_distributions = defaultdict(list)  # Keep track of files assigned to each node
    node_index = 0  # Start with the first node

    # Iterate over each C number and its corresponding files
    for c_number, files in grouped_files.items():
        # Get the base filename from one of the files
        p = Path(files[0])
        # base_filename = p.stem.replace("_T000", "")
        # new_stem = p.stem.replace("_T000", "_subband")
        fname = p.name
        txt_filename = fname.replace("T000_C", "SB").replace(
            ".MS", ".txt"
        )  # f"{new_stem}.txt"

        # Determine the current node's path for saving the JSON file
        node_path = nodes[node_index]
        node_dir = node_path + datapath  # os.path.join(node_path, datapath) doesnt work

        # os.makedirs(node_dir, exist_ok=True)

        txt_filepath = node_dir + "/" + txt_filename

        with open(txt_filepath, "w") as txt_file:
            txt_file.write("\n".join(files) + "\n")

        # Assign this JSON file to a node
        node_distributions[nodes[node_index]].append(txt_filepath)

        # Move to the next node if needed
        if len(node_distributions[nodes[node_index]]) >= max_txt_files_per_node:
            node_index = node_index + 1

    return node_distributions


def main(msin, datapath, nmses_per_node, from_nodes, to_nodes):
    all_msfiles = []
    for node in from_nodes:
        all_msfiles += glob(f"/net/node{node}/{datapath}/{msin}_T???_C???.MS")

    all_msfiles = sorted(all_msfiles)
    logging.info(f"Found {len(all_msfiles)} Total files")

    t_numbers, c_numbers = extract_T_and_C_numbers(all_msfiles)

    logging.info(f"Found {len(t_numbers)} Timechunks IDS: {t_numbers}")
    logging.info(f"Found {len(c_numbers)} Subband IDS: {c_numbers}")

    grouped_files = group_files_by_C(all_msfiles)

    logging.info(f"Made {len(grouped_files)} timechunk groups")

    assert len(to_nodes) * nmses_per_node >= len(
        grouped_files
    ), "Not enough nodes to distribute all subbands!"

    node_distributions = write_txt_files_to_nodes(
        grouped_files, datapath, to_nodes, nmses_per_node
    )

    # Print the distribution of JSON files across nodes
    for node, txt_files in node_distributions.items():
        logging.info(f"Node {node} has the following txt files:")
        for json_file in txt_files:
            logging.info(f"  - {json_file}")


if __name__ == "__main__":
    args = parser.parse_args()
    assert os.path.isdir(args.datapath)
    to_nodes = [f"/net/node{n}/" for n in args.to_nodes]
    main(args.msin, args.datapath, args.nmses_per_node, args.from_nodes, to_nodes)
