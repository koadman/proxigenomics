#!/bin/bash

#
# Map sequences in simulation
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -N MAPJOB

BWAEXE=$HOME/bin/bwa-0.7.6a/bwa

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ]
	then
		echo "Usage: <subject> <query> <output>"
		exit 1
	fi
	echo "Submitting run"
	qsub -W block=true -v SUBJECT=$1,QUERY=$2,OUTPUT=$3 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	# map
	$BWAEXE mem -t 2 $SUBJECT $QUERY > $OUTPUT
fi
