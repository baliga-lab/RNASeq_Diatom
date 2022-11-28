#!/usr/bin/env python3

import jinja2
import os
import argparse
import json

TEMPLATE = """#!/bin/bash

#SBATCH -J star_salmon_{{genome}}
#SBATCH -o /proj/omics4tb2/wwu/slurm_logs/"%j".out
#SBATCH -e /proj/omics4tb2/wwu/slurm_logs/"%j".out

{{conda_activate}}
{{sbatch_options}}

data_folder=$1
star_prefix="star_{{star_options.outFilterMismatchNmax}}_{{star_options.outFilterMismatchNoverLmax}}_{{star_options.outFilterScoreMinOverLread}}_{{star_options.outFilterMatchNmin}}{{dedup_prefix}}"
salmon_prefix="salmon_{{star_options.outFilterMismatchNmax}}_{{star_options.outFilterMismatchNoverLmax}}_{{star_options.outFilterScoreMinOverLread}}_{{star_options.outFilterMatchNmin}}{{dedup_prefix}}"


#$GS_HOME/code/rnaseq/run_star_salmon.py {{genome_dir}} {{input_dir}} $data_folder {{output_dir}} --outFilterMismatchNmax {{star_options.outFilterMismatchNmax}} --outFilterMismatchNoverLmax {{star_options.outFilterMismatchNoverLmax}} --outFilterScoreMinOverLread {{star_options.outFilterScoreMinOverLread}} --outFilterMatchNmin {{star_options.outFilterMatchNmin}} {{dedup_option}} --starPrefix $star_prefix --salmonPrefix $salmon_prefix
$GS_HOME/code/rnaseq/run_star_salmon_old.py {{genome_dir}} {{input_dir}} $data_folder {{output_dir}}
"""

DESCRIPTION = """make_star_salmon_job.py - Create STAR Salmon job file for Slurm"""

def make_sbatch_options(config):
    result = ""
    for option in config['sbatch_options']['star_salmon']:
        result += "#SBATCH %s\n" % option
    return result

def make_conda_activate(config):
    try:
		result = config['sbatch_options']['conda_env']
		if result != '':
			result = 'conda activate %s' % result
	except:
		result = ''
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
    config['sbatch_options'] = make_sbatch_options(config)
    config['conda_activate'] = make_conda_activate(config)

    config['dedup_prefix'] = '_dedup' if config['deduplicate_bam_files'] else ''
    config['dedup_option'] = '--dedup' if config['deduplicate_bam_files'] else ''
    print(templ.render(config))
