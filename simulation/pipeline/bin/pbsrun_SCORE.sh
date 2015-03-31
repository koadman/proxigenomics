#!/bin/bash

#
# Map sequences in simulation
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N SCOREJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

F1SCORE=bin/f1score.py
VMEASURE=bin/vmeasure.py
BCUBED=bin/bcubed.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 2 ]
	then
		echo "Usage: <truth table> <clustering>"
		exit 1
	fi
	echo "Submitting run"

    TARGETS=( ${CLUSTERING}.f1, ${CLUSTERING}.vm )
	trap 'rollback_rm_files "${TARGETS[@]}"; exit $?' INT TERM EXIT
	qsub -W block=true -v TRUTH=$1,CLUSTERING=$2 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	$F1SCORE $TRUTH $CLUSTERING ${CLUSTERING}.f1
	$VMEASURE $TRUTH $CLUSTERING ${CLUSTERING}.vm
	$BCUBED -o ${CLUSTERING}.bc $TRUTH $CLUSTERING

fi
