#!/bin/bash

#
# Align query to subject
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N LASTJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

BT2_HOME=$HOME/software/bowtie2-2.2.2/
LAP_CALC=$HOME/git/lap_release_1.1/aligner/calc_prob.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 4 ]
	then
		echo "Usage: <REFSEQ> <R1> <R2> <output>"
		exit 1
	fi

	echo "Submitting run"
	trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -W block=true -v ASM=$1,R1=$2,R2=$3,OUTPUT=$4 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	# calc probs
	$LAP_CALC -a $ASM -1 $R1 -2 $R2 -q -I 0 -o fr -X 600 -m 450 -t 100 -p 2 > $OUTPUT
fi
