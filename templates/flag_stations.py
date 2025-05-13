#! /usr/bin/env python3

import os
import logging
import logging.config

import casacore.tables as tab
from argparse import ArgumentParser

logging.config.fileConfig(os.path.join(os.path.dirname(__file__), 'logging.config'))
logging = logging.getLogger('root')


parser = ArgumentParser("get clip value from XY,YX rms of data")
parser.add_argument(
    "-i",
    "--ms_list",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="ms_list",
    default="",
)

parser.add_argument(
    "-s", "--stations", help="Stations to flag default=None", dest="stations", default="None"
)
parser.add_argument(
    "-l", "--flag_longbaselines", help="flag stations outside 10km of core", action="store_true",
)
parser.add_argument(
    "-n", "--flag_international_stations", help="flag international stations", action="store_true"
)
parser.add_argument(
    "-r", "--flag_intrastations", help="flag CS CS stations", action="store_true"
)


def parse_ms_list(ms_list_or_txt):
    ms_list = [i for i in ms_list_or_txt if tab.tableexists(i)]
    if len(ms_list) == 0:
        ms_lists = []
        for ifile in ms_list_or_txt:
            if os.path.isfile(ifile):
                myf = open(ifile)
                ms_lists.append([i.strip() for i in myf if tab.tableexists(i.strip())])
    else:
        ms_lists = [ms_list]

    return [ms for ms_list in ms_lists for ms in ms_list]


def main(args):
    logging.info("Checking measurement set(s)")
    ms_list = parse_ms_list(args.ms_list)

    logging.info(f"Total measurement sets found: {len(ms_list)}")

    for i in ms_list:
        # TODO: cant find a way to dynamically change the station  to flag within the r'' line
        # if args.stations:
        #     logging.info(f"Flagging stations")
        #         tab.taql(
        #             r'UPDATE $i SET FLAG=True WHERE mscal.baseline("RS208HBA")'
        #         )

        if args.flag_international_stations:
            logging.info("Flagging all non-Dutch stations")
            tab.taql(r'UPDATE $i SET FLAG=True WHERE mscal.baseline("![CR]*&&")')

        if args.flag_intrastations:
            logging.info("Flagging intrastrations")
            tab.taql(
                r'UPDATE $i SET FLAG=True WHERE mscal.baseline("/(.*)HBA0&\1HBA1/")'
            )

        if args.flag_longbaselines:  #TODO: The stations flagged here should change based on the frequency of the MS
            logging.info(f"Flagging long baselines")
            tab.taql(r'UPDATE $i SET FLAG=True WHERE mscal.baseline("RS406HBA")')


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
