#!/bin/bash

#
# Create homology search indexes
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N INDEXJOB

BWA=$HOME/bin/bwa-0.7.6a/bwa
LASTDB=$HOME/bin/lastdb

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
	$BWA index $FASTA
	$LASTDB $FASTA $FASTA
fi
