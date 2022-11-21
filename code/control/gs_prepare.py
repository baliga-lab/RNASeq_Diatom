#!/usr/bin/env python3 

"""
Script to check parameters and prepare for submission
"""
import argparse
import json
import os, sys, subprocess


DESCRIPTION = """gs_prepare.py - prepare data for workflow submission"""

def check_params(config):
    """ensures integrity of the parameters"""
    print("checking integrity of parameters...", end="")
    if not os.path.exists(config["input_dir"]):
        sys.exit("Input directory '%s' does not exist" % config["input_dir"])

    if len(config["includes"]) > 0:
        # check the existence of the included directories
        for incl in config["includes"]:
            inp_dir = os.path.join(config["input_dir"], incl)
            if not os.path.exists(inp_dir):
                sys.exit("Input directory '%s' does not exist" % inp_dir)
    print("done")


def create_dirs(config):
    if not os.path.exists(config["output_dir"]):
        print("creating output directory '%s'" % config["output_dir"])
        os.makedirs(config["output_dir"])


def check_salmon():
    try:
        print("checking for salmon... ", end="")
        compl_proc = subprocess.run(["salmon", "-v"], check=True, capture_output=True)
        progname, version = compl_proc.stdout.decode('utf-8').strip().split()
        print("(found version '%s') ..." % version, end="")
        if version != "0.13.1":
            sys.exit("Unsupported version %s. Currently, only salmon 0.13.1 is supported" % version)
        print("done")
    except FileNotFoundError:
        sys.exit("Can not find salmon (not installed or not in PATH)")


def check_star():
    try:
        print("checking for STAR... ", end="")
        compl_proc = subprocess.run(["STAR", "--version"], check=True, capture_output=True)
        version = compl_proc.stdout.decode('utf-8').strip()
        print("(found version '%s') ..." %  version, end="")
        print("done")
    except FileNotFoundError:
        sys.exit("Can not find STAR (not installed or not in PATH)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)    
    parser.add_argument('configfile', help="configuration file")
    args = parser.parse_args()

    check_star()
    check_salmon()
    
    with open(args.configfile) as infile:
        config = json.load(infile)
    check_params(config)
    create_dirs(config)

