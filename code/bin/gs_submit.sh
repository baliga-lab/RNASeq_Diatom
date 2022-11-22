#!/bin/bash

if [ "$#" -ne 1 ] ; then
	echo "Usage: gs_submit.sh <config-file>"
	exit 1
fi

if [ -z "$GS_HOME" ]; then
	echo "Variable GS_HOME is not set. Please set to Global_Search repository path."
	exit 1
fi

$GS_HOME/code/control/gs_prepare.py $1
