#!/usr/bin/env python3

"""
Post-processing step for STAR Salmon analysis.

This includes

  1. extracting TPM and read information from salmon quant
  2. running MultiQC
"""
import argparse
import json
from rpy2.robjects.packages import importr
import subprocess
import os

DESCRIPTION = """post_star_salmon.py - Post-run step for STAR Salmon"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('configfile', help='configuration file')
    args = parser.parse_args()
    with open(args.configfile) as infile:
        config = json.load(infile)
    output_dir = config['output_dir']
    postrun_outdir = os.path.join(output_dir, "Post_Run_Results")
    genome_dir = config['genome_dir']
    org1, org2 = os.path.basename(genome_dir).split('_')
    print("postrun_outdir: '%s', org1: %s, org2: %s" % (postrun_outdir, org1, org2))
    global_search = importr("GlobalSearch")
    global_search.extract_salmon_quants(org1, org2, output_dir, postrun_outdir)
