#!/usr/bin/env python3

import jinja2
import os
import argparse
import json

from globalsearch.rnaseq.find_files import rnaseq_data_folder_list

TEMPLATE = """#!/bin/bash

#SBATCH -J kallisto_{{genome}}
#SBATCH -o {{log_dir}}/"%j".out
#SBATCH -e {{log_dir}}/"%j".out
#SBATCH --array={{array_range}}

{{sbatch_options_comments}}

echo "ARRAY TASK ID: $SLURM_ARRAY_TASK_ID"
data_folders=({{data_folders}})
data_folder=${data_folders[$SLURM_ARRAY_TASK_ID]}

{{sbatch_extras}}

python3 -m globalsearch.rnaseq.run_kallisto {{genome_dir}} {{input_dir}} $data_folder {{genome_fasta}} {{output_dir}}
"""

DESCRIPTION = """make_kallisto_job.py - Create Kallisto job file for Slurm"""


def make_sbatch_options(config):
    result = ""
    for option in config['sbatch_options']['kallisto']['options']:
        result += "#SBATCH %s\n" % option
    return result

def make_sbatch_extras(config):
    result = ""
    for extra in config['sbatch_options']['kallisto']['extras']:
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
    config['sbatch_options_comments'] = make_sbatch_options(config)

    data_folders = rnaseq_data_folder_list(config)
    config["data_folders"] = ' '.join(data_folders)

    # Array specification
    try:
        array_max_tasks = config['sbatch_options']['array_max_tasks']
    except:
        array_max_tasks = 0
    if array_max_tasks > 0:
        array_max_task_spec = "%%%d" % array_max_tasks
    else:
        array_max_task_spec = ""

    config["array_range"] = "0-%d%s" % (len(data_folders) - 1, array_max_task_spec)
    print(templ.render(config))
