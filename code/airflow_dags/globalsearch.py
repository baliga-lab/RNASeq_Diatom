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

# get the airflow.task logger
task_logger = logging.getLogger('airflow.task')


LOG_DIR = "/Users/weiju/tmp"
OUT_DIR = "/Users/weiju/tmp"

ERROR_LOG = os.path.join(LOG_DIR, "globalsearch_errors.log")
STATUS_LOG = os.path.join(LOG_DIR, "globalseach_status.log")


# DEFAULT VALUES IS the example files
CONFIG_DIR = "/Users/weiju/Projects/ISB/Global_Search/examples"
CONFIG_FILE = os.path.join(CONFIG_DIR, "redsea-001.json")


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

    # we can obtain the run id through {{run_id}}
    gs_prepare = BashOperator(
        task_id='prepare',
        bash_command='echo "Preparing parameters for {{params.config_file}}" >> %s 2>&1' % (STATUS_LOG),
    )
    decide_algo = BranchPythonOperator(
        task_id="decide_algorithm",
        provide_context=True,
        python_callable=_choose_algorithm
    )
    index_star = BashOperator(
        task_id='idx_genome_star',
        bash_command='echo "Index the genome...done" >> %s 2>&1' % (STATUS_LOG),
    )
    star_salmon = BashOperator(
        task_id='starsalmon',
        bash_command='echo "Submitting STAR/salmon jobs...done" >> %s 2>&1' % (STATUS_LOG),
    )

    spladder = BashOperator(
        task_id='spladder',
        bash_command='echo "Submitting SplAdder jobs...done" >> %s 2>&1' % (STATUS_LOG),
    )

    orthofinder = BashOperator(
        task_id='orthofinder',
        bash_command='echo "Submitting Orthofinder jobs...done" >> %s 2>&1' % (STATUS_LOG),
    )

    process_star_results = BashOperator(
        task_id='proc_star_results',
        bash_command='echo "Processing STAR results...done" >> %s 2>&1' % (STATUS_LOG),
    )

    process_spladder_results = BashOperator(
        task_id='proc_spladder_results',
        bash_command='echo "Processing SplAdder results...done" >> %s 2>&1' % (STATUS_LOG),
    )

    process_orthofinder_results = BashOperator(
        task_id='proc_orthofinder_results',
        bash_command='echo "Processing Orthofinder results...done" >> %s 2>&1' % (STATUS_LOG),
    )


    index_kallisto = BashOperator(
        task_id='idx_genome_kallisto',
        bash_command='echo "Index the genome...done" >> %s 2>&1' % (STATUS_LOG),
    )
    kallisto = BashOperator(
        task_id='kallisto',
        bash_command='echo "Submitting Kallisto jobs...done" >> %s 2>&1' % (STATUS_LOG),
    )
    process_kallisto_results = BashOperator(
        task_id='proc_kallisto_results',
        bash_command='echo "Processing Kallisto results...done" >> %s 2>&1' % (STATUS_LOG),
    )

    dag.doc_md = __doc__
    dag.doc_md = """
    This is the Airflow implementation of the GlobalSearch pipeline. This is a parameterized pipeline, so
    please trigger this manual with your input data paths as your parameters.
    """

    gs_prepare >> decide_algo >> [index_star, index_kallisto]
    index_star >> [star_salmon, orthofinder]
    star_salmon >> process_star_results
    star_salmon >> spladder >> process_spladder_results
    orthofinder >> process_orthofinder_results
    index_kallisto >> kallisto >> process_kallisto_results

