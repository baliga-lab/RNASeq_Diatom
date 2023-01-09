#!/usr/bin/env python3

from datetime import datetime, timedelta
from textwrap import dedent

from airflow import DAG
from airflow.models.param import Param

from airflow.operators.python import BranchPythonOperator
from airflow.operators.bash import BashOperator

from airflow.decorators import task

import os
import json
from pprint import pprint

import logging

from globalsearch.rnaseq.run_star_salmon import create_genome_index as idx_star_genome
from globalsearch.rnaseq.run_kallisto import kallisto_index as idx_kallisto_genome

# get the airflow.task logger
task_logger = logging.getLogger('airflow.task')


# DEFAULT VALUES IS the example files
CONFIG_DIR = "/Users/weiju/Projects/ISB/Global_Search/examples"
CONFIG_FILE = os.path.join(CONFIG_DIR, "redsea-001.json")
DEBUGGING = True

with DAG(
        'GlobalSearch',
        default_args={
            'depends_on_past': False,
            'email': ['wwu@systemsbiology.org'],
            'email_on_failure': False,
            'email_on_retry': False,
            'retries': 1,
            'retry_delay': timedelta(minutes=5),
        },
        description='GlobalSearch Pipeline DAG',
        schedule_interval=None,
        start_date=datetime(2022, 12, 1),
        catchup=False,
        tags=['globalsearch', 'star', 'salmon'],
        params={
            'config_file': Param(type='string', default=CONFIG_FILE),
            'debugging': Param(type='boolean', default=DEBUGGING)
        }
) as dag:


    def _choose_algorithm(**kwargs):
        """
        decide the RNA Seq method based on the rnaseq_algorithm parameter in config file
        """
        run_id = kwargs['run_id']
        config_file = kwargs['params']['config_file']
        with open(config_file) as infile:
            config = json.load(infile)
        if config['rnaseq_algorithm'] == "star_salmon":
            return "idx_genome_star"
        else:
            return "idx_genome_kallisto"

    gs_prepare = BashOperator(
        task_id='prepare',
        bash_command='{% if params.debugging %}echo "gs_prepare {{params.config_file}}"{% else %}gs_prepare {{params.config_file}}{% endif %}'
    )

    decide_algo = BranchPythonOperator(
        task_id="decide_algorithm",
        provide_context=True,
        python_callable=_choose_algorithm
    )

    """
    STAR Based Genome Indexing
    """
    @task(task_id="idx_genome_star")
    def idx_genome_star_task(ds=None, **kwargs):
        run_id = kwargs['run_id']
        config_file = kwargs['params']['config_file']
        with open(config_file) as infile:
            config = json.load(infile)
        genome_dir = config['genome_dir']
        genome_fasta = config['genome_fasta']
        task_logger.info("start indexing the genome for STAR")
        if not kwargs['params']['debugging']:
            idx_star_genome(genome_dir, genome_fasta)
        task_logger.info("finished indexing the genome for STAR")
        return "Indexed STAR in dir: '%s' on file '%s' finished" % (genome_dir, genome_fasta)

    idx_genome_star = idx_genome_star_task()

    """
    Kallisto Based Genome Indexing
    """
    @task(task_id="idx_genome_kallisto")
    def idx_genome_kallisto_task(ds=None, **kwargs):
        config_file = kwargs['params']['config_file']
        with open(config_file) as infile:
            config = json.load(infile)
        genome_dir = config['genome_dir']
        genome_fasta = config['genome_fasta']
        organism = os.path.basename(genome_dir)
        # make the index in the genome directory for now
        index_path = os.path.join(genome_dir, '%s_kallistoindex' % organism)
        task_logger.info("start indexing the genome for STAR")
        if not kwargs['params']['debugging']:
            idx_kallisto_genome(index_path, genome_fasta)
        return "Indexed Kallisto"

    idx_genome_kallisto = idx_genome_kallisto_task()


    star_salmon = BashOperator(
        task_id='starsalmon',
        bash_command='echo "Submitting STAR/salmon jobs...done"',
    )

    spladder = BashOperator(
        task_id='spladder',
        bash_command='echo "Submitting SplAdder jobs...done"',
    )

    orthofinder = BashOperator(
        task_id='orthofinder',
        bash_command='echo "Submitting Orthofinder jobs...done"',
    )

    process_star_results = BashOperator(
        task_id='proc_star_results',
        bash_command='echo "Processing STAR results...done"',
    )

    process_spladder_results = BashOperator(
        task_id='proc_spladder_results',
        bash_command='echo "Processing SplAdder results...done"',
    )

    process_orthofinder_results = BashOperator(
        task_id='proc_orthofinder_results',
        bash_command='echo "Processing Orthofinder results...done"',
    )


    kallisto = BashOperator(
        task_id='kallisto',
        bash_command='echo "Submitting Kallisto jobs...done"',
    )
    process_kallisto_results = BashOperator(
        task_id='proc_kallisto_results',
        bash_command='echo "Processing Kallisto results...done"',
    )

    dag.doc_md = __doc__
    dag.doc_md = """
    This is the Airflow implementation of the GlobalSearch pipeline. This is a parameterized pipeline, so
    please trigger this manual with your input data paths as your parameters.
    """

    gs_prepare >> decide_algo >> [idx_genome_star, idx_genome_kallisto]
    idx_genome_star >> [star_salmon, orthofinder]
    star_salmon >> process_star_results
    star_salmon >> spladder >> process_spladder_results
    orthofinder >> process_orthofinder_results
    idx_genome_kallisto >> kallisto >> process_kallisto_results

