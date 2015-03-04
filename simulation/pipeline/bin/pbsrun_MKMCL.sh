#!/bin/bash

#
# Create mcl input file
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N MKMCLJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

MKMCL=bin/makeMCLinput.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then

	if [ $# -ne 4 ]
	then
		echo "Usage: [min length] [edges csv] [nodes csv] [output]"
		exit 1
	fi

	echo "Submitting run"
	trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -W block=true -v MINLEN=$1,EDGES=$2,NODES=$3,OUTPUT=$4 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

    $MKMCL $MINLEN $EDGES $NODES $OUTPUT

fi
