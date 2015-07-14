#!/bin/bash

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N LOUVAINJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
    BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    source $BINDIR/bash_init.sh
fi

LOUVAIN=bin/louvain_cluster.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then

    if [ $# -ne 3 ]
    then
		echo "Usage: <otype> <input> <output>"
		exit 1
	fi
	echo "Submitting run"

	trap 'rollback_rm_file $3; exit $?' INT TERM EXIT
	qsub -W block=true -v OTYPE=$1,INPUT=$2,OUTPUT=$3 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	if [ -s $INPUT ]
	then
		$LOUVAIN --otype $OTYPE $INPUT $OUTPUT
	else
		echo "Input had no data" > $OUTPUT
	fi

fi
