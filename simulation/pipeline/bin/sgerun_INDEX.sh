#!/bin/bash

#
# Create homology search indexes
#

#$ -o logs
#$ -e logs
#$ -cwd
#$ -N INDEXJOB

BWA=$PWD/external/a5_miseq_linux_20140604/bin/bwa

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 1 ]
	then
		echo "Usage: <fasta>"
		exit 1
	fi
	echo "Submitting run"
	CMD=`readlink -f $0`
	qsub -sync yes -b n -v FASTA=$1 $CMD
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	
	# create indexes
	$BWA index $FASTA
fi
