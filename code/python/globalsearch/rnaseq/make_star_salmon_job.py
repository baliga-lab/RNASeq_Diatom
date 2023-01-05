#!/usr/bin/env python3

import jinja2
import os
import argparse
import json

TEMPLATE = """#!/bin/bash

#SBATCH -J star_salmon_{{genome}}
#SBATCH -o {{log_dir}}/"%j".out
#SBATCH -e {{log_dir}}/"%j".out
{{sbatch_options}}

data_folder=$1
star_prefix="star_{{star_options.outFilterMismatchNmax}}_{{star_options.outFilterMismatchNoverLmax}}_{{star_options.outFilterScoreMinOverLread}}_{{star_options.outFilterMatchNmin}}{{dedup_prefix}}"
salmon_prefix="salmon_{{star_options.outFilterMismatchNmax}}_{{star_options.outFilterMismatchNoverLmax}}_{{star_options.outFilterScoreMinOverLread}}_{{star_options.outFilterMatchNmin}}{{dedup_prefix}}"

{{sbatch_extras}}

$GS_HOME/code/python/globalsearch/rnaseq/run_star_salmon.py {{twopass_mode}} {{fastq_patterns}} --outSAMattributes {{star_options.outSAMattributes}} --outFilterMismatchNmax {{star_options.outFilterMismatchNmax}} --outFilterMismatchNoverLmax {{star_options.outFilterMismatchNoverLmax}} --outFilterScoreMinOverLread {{star_options.outFilterScoreMinOverLread}} --outFilterMatchNmin {{star_options.outFilterMatchNmin}} {{dedup_option}} --starPrefix $star_prefix --salmonPrefix $salmon_prefix {{genome_gff_option}} {{genome_fasta_option}} {{genome_dir}} {{input_dir}} $data_folder {{output_dir}}
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

    config['dedup_prefix'] = '_dedup' if config['deduplicate_bam_files'] else ''
    config['dedup_option'] = '--dedup' if config['deduplicate_bam_files'] else ''
    config['dedup_option'] = '--dedup' if config['deduplicate_bam_files'] else ''
    config['twopass_mode'] = '--twopassMode' if config['star_options']['twopassMode'] else ''
    config['fastq_patterns'] = '--fastq_patterns "%s"' % ','.join(config['fastq_patterns']) if len(config['fastq_patterns']) > 0 else ''

    # see if optional genome_gff exists
    try:
        config['genome_gff_option'] = ''
        genome_gff = config['genome_gff']
        if os.path.exists(genome_gff):
            config['genome_gff_option'] = '--genome_gff %s' % genome_gff
    except:
        pass

    # see if optional genome_fasta exists
    try:
        config['genome_fasta_option'] = ''
        genome_fasta = config['genome_fasta']
        if os.path.exists(genome_fasta):
            config['genome_fasta_option'] = '--genome_fasta %s' % genome_fasta
    except:
        pass
    print(templ.render(config))
