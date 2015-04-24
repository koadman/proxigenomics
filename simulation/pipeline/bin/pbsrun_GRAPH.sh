#!/bin/bash

#
# Create graph files
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N GRAPHJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

BAMTOEDGES=bin/bamToEdges.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	
	if [ $# -ne 3 ]
	then
		echo "Usage: [hic2ctg.bam] [edge out] [node out]"
		exit 1
	fi

	echo "Submitting run"
	TARGET=( $2 $3 )
	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	qsub -W block=true -v HIC2CTG=$1,EDGES=$2,NODES=$3 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	$BAMTOEDGES $HIC2CTG $EDGES $NODES
	
fi
