#!/usr/bin/env python3

"""
Generate a SLURM job file to call the STAR indexer
It accepts a pipeline configuration file and
fills out the template with the appropriate parameters
"""

import jinja2
import os
import argparse
import json


TEMPLATE = """#!/bin/bash

#SBATCH -J star_salmon_{{genome}}
#SBATCH -o {{log_dir}}/"%j".out
#SBATCH -e {{log_dir}}/"%j".out

{{sbatch_options}}

echo "TASK ID: $SLURM_JOB_ID"

{{sbatch_extras}}

python3 -m globalsearch.rnaseq.index_star {{star_index_cmd_options}} {{genome_fasta_option}} {{genome_dir}}
"""

DESCRIPTION = """make_star_salmon_job.py - Create STAR Salmon job file for Slurm"""

def make_sbatch_options(config):
    result = ""
    for option in config['sbatch_options']['star_salmon']['options']:
        result += "#SBATCH %s\n" % option
    return result

def make_sbatch_extras(config):
    result = ""
    for extra in config['sbatch_options']['star_salmon']['extras']:
        result += "%s\n" % extra
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('configfile', help="configuration file")
    args = parser.parse_args()
    with open(args.configfile) as infile:
        config = json.load(infile)

    templ = jinja2.Template(TEMPLATE)
    genome = os.path.basename(os.path.normpath(config['genome_dir']))
    config['genome'] = genome
    config['sbatch_extras'] = make_sbatch_extras(config)
    config['sbatch_options'] = make_sbatch_options(config)

    # see if optional genome_fasta exists
    try:
        config['genome_fasta_option'] = ''
        genome_fasta = config['genome_fasta']
        if os.path.exists(genome_fasta):
            config['genome_fasta_option'] = '--genome_fasta %s' % genome_fasta
    except:
        pass

    # see if optional star_index_options exists
    try:
        config['star_index_cmd_options'] = ''
        star_index_options = config['star_index_options']
        options = []
        try:
            options.append("--runThreadN %d" % star_index_options['runThreadN'])
        except KeyError:
            pass
        try:
            options.append("--genomeChrBinNbits %d" % star_index_options['genomeChrBinNbits'])
        except KeyError:
            pass
        try:
            options.append("--genomeSAindexNbases %d" % star_index_options['genomeSAindexNbases'])
        except KeyError:
            pass

        try:
            options.append("--sjdbGTFfeatureExon %s" % star_index_options['sjdbGTFfeatureExon'])
        except KeyError:
            pass
        try:
            options.append("--sjdbGTFtagExonParentTranscript %s" % star_index_options['sjdbGTFtagExonParentTranscript'])
        except KeyError:
            pass
        try:
            options.append("--sjdbGTFtagExonParentGene %s" % star_index_options['sjdbGTFtagExonParentGene'])
        except KeyError:
            pass


        config['star_index_cmd_options'] = ' '.join(options)
    except KeyError:
        pass

    print(templ.render(config))
