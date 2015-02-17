#!/bin/bash

#
# Map sequences in simulation
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N SCOREJOB


JOINTRUTH=bin/mclJoinWithTruth.py
F1SCORE=bin/f1score.py
VMEASURE=bin/vmeasure.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 2 ]
	then
		echo "Usage: <truth table> <clustering>"
		exit 1
	fi
	echo "Submitting run"
	qsub -W block=true -v TRUTH=$1,CLUSTERING=$2 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	$F1SCORE $TRUTH $CLUSTERING ${CLUSTERING}.f1
	$VMEASURE $TRUTH $CLUSTERING ${CLUSTERING}.vm
fi
