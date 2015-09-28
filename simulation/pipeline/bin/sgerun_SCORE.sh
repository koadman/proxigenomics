#!/bin/bash

#
# Map sequences in simulation
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N SCOREJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

F1SCORE=bin/f1score.py
VMEASURE=bin/vmeasure.py
BCUBED=bin/bcubed.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 2 ]
	then
		echo "Usage: <truth table> <clustering>"
		exit 1
	fi
	echo "Submitting run"

    TARGETS=( ${CLUSTERING}.f1, ${CLUSTERING}.vm, ${CLUSTERING}.bc )
	#trap 'rollback_rm_files "${TARGETS[@]}"; exit $?' INT TERM EXIT
	qsub -sync yes -V -v TRUTH=$1,CLUSTERING=$2 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	# skipping F1 due to computability problems -- runtimes for complicated assemblies
	#$F1SCORE -s 1500 $TRUTH $CLUSTERING ${CLUSTERING}.f1
	touch ${CLUSTERING}.f1
	$VMEASURE $TRUTH $CLUSTERING ${CLUSTERING}.vm
	$BCUBED -o ${CLUSTERING}.bc $TRUTH $CLUSTERING

fi
