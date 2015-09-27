#!/bin/bash

#
# Align query to subject
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N LAPJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

export BT2_HOME=$ITHREE_GIT/software/bowtie2-2.2.2/
LAP_CALC=$HOME/git/lap_release_1.1/aligner/calc_prob.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 4 ]
	then
		echo "Usage: <REFSEQ> <R1> <R2> <output>"
		exit 1
	fi

	echo "Submitting run"
	#trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -sync yes -V -pe smp 4 -v ASM=$1,R1=$2,R2=$3,OUTPUT=$4 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	# calc probs
	$LAP_CALC -a $ASM -1 $R1 -2 $R2 -q -I 0 -o fr -X 650 -m 450 -t 100 -p $NSLOTS > $OUTPUT
fi
