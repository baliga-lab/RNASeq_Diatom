#!/usr/bin/env python3 

"""
Script to check parameters and prepare for submission
"""
import argparse
import json
import os, sys, glob, subprocess


DESCRIPTION = """gs_prepare.py - prepare data for workflow submission"""

def check_params(config):
    """ensures integrity of the parameters"""
    print("checking integrity of parameters...", end="")
    if not os.path.exists(config["input_dir"]):
        sys.exit("Input directory '%s' does not exist" % config["input_dir"])

    try:
        if len(config['genome_fasta'].strip()) > 0 and not os.path.exists(config["genome_fasta"]):
            sys.exit("ERROR: Specified Genome FASTA '%s' does not exist" % config["genome_fasta"])
    except:
        print("WARNING: genome_gff not specified (htseq will be skipped)")

    if not os.path.exists(config["genome_dir"]):
        sys.exit("Genome directory '%s' does not exist" % config["genome_dir"])
        genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))
        if len(genome_fasta) == 0:
            sys.exit("Genome directory '%s' does not contain any FASTA files" % config["genome_dir"])
    try:
        if not os.path.exists(config["genome_gff"]):
            print("WARNING: Genome GFF '%s' does not exist (htseq will be skipped)" % config["genome_gff"])
    except:
        print("WARNING: genome_gff not specified (htseq will be skipped)")

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
    if not os.path.exists(config["log_dir"]):
        print("creating log directory '%s'" % config["log_dir"])
        os.makedirs(config["log_dir"])


def __check_command(command, num_info_components=1, check_version=None, version_switch='--version',
                    multiline=False):
    """Generic command checker, can check for different version info formats and
    restrict version numbers
    """
    try:
        print("checking for %s... " % command, end="")
        compl_proc = subprocess.run([command, version_switch], check=True, capture_output=True)
        if multiline:
            info_string = compl_proc.stdout.decode('utf-8').split('\n')[0].strip()
        else:
            info_string = compl_proc.stdout.decode('utf-8').strip()
        if num_info_components == 2:
            progname, version = info_string.split()
        else:
            version = info_string

        print("(found version '%s') ..." % version, end="")
        if check_version is not None and check_version != version:
            sys.exit("Unsupported version %s. Currently, only %s %s is supported" % (command, check_version))
        print("done")
    except FileNotFoundError:
        sys.exit("Can not find %s (not installed or not in PATH)" % command)


def check_salmon():
    __check_command("salmon", num_info_components=2, check_version="0.13.1")


def check_star():
    __check_command("STAR")


def check_htseq():
    __check_command("htseq-count")


def check_samtools():
    __check_command("samtools", multiline=True)


def check_kallisto():
    __check_command("kallisto", version_switch="version")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('configfile', help="configuration file")
    args = parser.parse_args()

    check_star()
    check_salmon()
    check_htseq()
    check_samtools()

    with open(args.configfile) as infile:
        config = json.load(infile)
    check_params(config)
    create_dirs(config)

