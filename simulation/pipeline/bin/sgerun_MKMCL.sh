#!/bin/bash

#
# Create mcl input file
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N MKMCLJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

MKMCL=bin/makeMCLinput.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

	if [ $# -ne 4 ]
	then
		echo "Usage: [min length] [edges csv] [nodes csv] [output]"
		exit 1
	fi

	echo "Submitting run"
	#trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -sync yes -V -v MINLEN=$1,EDGES=$2,NODES=$3,OUTPUT=$4 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

    $MKMCL $MINLEN $EDGES $NODES $OUTPUT

fi
