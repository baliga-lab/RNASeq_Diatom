#!/bin/bash

if [ "$#" -ne 1 ] ; then
	echo "Usage: gs_submit.sh <config-file>"
	exit 1
fi

./gs_prepare.py $1
