#!/bin/bash

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

#
# Create a minia assembly
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N MINIAJOB

MINIA=$HOME/bin/minia-1.6088/minia

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 6 ]
	then
		echo "Usage: [kmer size] [min abun] [genome size] [R1] [R2] [out base]"
		exit 1
	fi
	echo "Submitting run"
	trap 'rollback_rm_dir $6; exit $?' INT TERM EXIT
	qsub -W block=true -v KMER=$1,ABUN=$2,GSIZE=$3,R1=$4,R2=$5,OUT_BASE=$6 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	echo -e "$R1\n$R2" > files.minia
	
	$MINIA files.minia $KMER $ABUN $GSIZE $OUT_BASE.minia
fi
