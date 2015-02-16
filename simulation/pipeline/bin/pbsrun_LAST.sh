#!/bin/bash

#
# Align query to subject
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N LASTJOB

LASTAL=$HOME/bin/lastal
MAFCONV=$HOME/bin/maf-convert

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

	# align
	$LASTAL -f 0 $SUBJECT $QUERY | $MAFCONV psl > $OUTPUT
fi
