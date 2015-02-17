#!/bin/bash

#
# Map sequences in simulation
#





#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N MAPJOB

BWAEXE=$HOME/bin/bwa-0.7.6a/bwa

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ] && [ $# -ne 4 ]
	then
		echo "Usage: <subject> <query> <output>"
		echo "or"
		echo "Usage: <subject> <query1> <query2> <output>"
		exit 1
	fi

	if [ $# -eq 3 ]
	then
    	echo "Submitting single query run"
	    qsub -W block=true -v SUBJECT=$1,QUERY1=$2,QUERY2="",OUTPUT=$3 $0
    elif [ $# -eq 4 ]
    then
    	echo "Submitting double query run"
        qsub -W block=true -v SUBJECT=$1,QUERY1=$2,QUERY2=$3,OUTPUT=$4 $0
    fi

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	# Map reads. The second query file can be empty
	$BWAEXE mem -t 2 $SUBJECT $QUERY1 $QUERY2 > $OUTPUT
fi
