#!/bin/bash

#
# create BWA index
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -N INDEXJOB

BWAEXE=$HOME/bin/bwa-0.7.6a/bwa

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 1 ]
	then
		echo "Usage: <fasta>"
		exit 1
	fi
	echo "Submitting run"
	qsub -W block=true -v FASTA=$1 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	# create indexes
	$BWAEXE index $FASTA
fi
