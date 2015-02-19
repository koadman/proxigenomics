#!/bin/bash

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

#
# Map 2 query sequences in simulation
#
# this is tedious as passing >1 query sequences to qsub is very awkward. I gave up with
# how complicated it was becoming.
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N MAPJOB

BWAEXE=$HOME/bin/bwa-0.7.6a/bwa

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 4 ]
	then
		echo "Usage: <subject> <query1> <query2> <output>"
		exit 1
	fi
	echo "Submitting run"
	trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -W block=true -v SUBJECT=$1,QUERY1=$2,QUERY2=$3,OUTPUT=$4 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	# map
	$BWAEXE mem -t 2 $SUBJECT $QUERY1 $QUERY2 > $OUTPUT
fi
