#!/bin/bash

#
# Deconvolute and build tree 
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N DECONVOLVEJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

SNVBPNMFT=$BINDIR/snvbpnmft.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

	# get full path of script
	CMD=`readlink -f $0`

    qsub -sync yes -V -v ARGS="$@" $CMD 
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	# run the algorithm
	$SNVBPNMFT $ARGS
fi
