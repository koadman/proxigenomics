#!/bin/bash

#
# Create mcl input file
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N MKMETISJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

CMD=bin/edgeToMetis.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

	if [ $# -ne 5 ]
	then
		echo "Usage: [min length] [edges csv] [nodes csv] [output] [nodemap]"
		exit 1
	fi

	echo "Submitting run"
	ONROLLBACK=($4 $5)
	#trap 'rollback_rm_files ${ONROLLBACK[@]}; exit $?' INT TERM EXIT
	qsub -sync yes -V -v MINLEN=$1,EDGES=$2,NODES=$3,OUTPUT=$4,NODEMAP=$5 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

    $CMD --fmt metis -m $MINLEN $EDGES $NODES $OUTPUT $NODEMAP

fi
