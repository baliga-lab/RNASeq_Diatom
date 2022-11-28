#!/bin/bash

if [ "$#" -ne 1 ] ; then
	echo "Usage: gs_submit.sh <config-file>"
	exit 1
fi

if [ -z "$GS_HOME" ]; then
	echo "Variable GS_HOME is not set. Please set to Global_Search repository path."
	exit 1
fi

tmpfile=$(mktemp /tmp/slurm_job.XXXXXX)
$GS_HOME/code/control/gs_prepare.py $1 && $GS_HOME/code/rnaseq/make_star_salmon_job.py $1 > $tmpfile && $GS_HOME/code/control/gs_submit.py $tmpfile $1

rm -f $tmpfile
