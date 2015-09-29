#!/bin/bash

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N LOUVAINJOB

if [ -z "$JOB_ID" ]
then
    BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    source $BINDIR/bash_init.sh
fi

LOUVAIN=bin/louvain_cluster.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

    if [ $# -ne 3 ]
    then
		echo "Usage: <otype> <input> <output>"
		exit 1
	fi
	echo "Submitting run"

	#trap 'rollback_rm_file $3; exit $?' INT TERM EXIT
	qsub -sync yes -V -v OTYPE=$1,INPUT=$2,OUTPUT=$3 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	if [ -s $INPUT ]
	then
		$LOUVAIN --otype $OTYPE $INPUT $OUTPUT
	else
		echo "Input had no data" > $OUTPUT
	fi

fi
