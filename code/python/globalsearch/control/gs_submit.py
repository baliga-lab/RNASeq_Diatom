#!/usr/bin/env python3

import argparse
import json
import os, sys, subprocess
from globalsearch.rnaseq.find_files import rnaseq_data_folder_list


DESCRIPTION = """gs_submit.py - submit job file"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('jobfile', help="SLURM job file")
    parser.add_argument('configfile', help="configuration file")
    args = parser.parse_args()

    with open(args.configfile) as infile:
        config = json.load(infile)
    data_folders = rnaseq_data_folder_list(config)
    print("Submitting %d data folders..." % len(data_folders))
    for data_folder in data_folders:
        command_str = 'sbatch %s %s' % (args.jobfile, data_folder)
        print("Submitting data folder '%s'" % data_folder)
        print(command_str)
        compl_proc = subprocess.run(['sbatch', args.jobfile, data_folder], check=True)


