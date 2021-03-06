#!/bin/bash

#
# Align query to subject
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N LAPSUMJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

export BT2_HOME=$ITHREE_GIT/software/bowtie2-2.2.6
LAP_SUM=$HOME/git/lap_release_1.1/aligner/sum_prob.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 4 ]
	then
		echo "Usage: <ref_in> <ctg_in>"
		exit 1
	fi

	echo "Submitting run"
	#trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -sync yes -V -v REFPROB=$1,CTGPROB=$2 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	# calc probs
	$LAP_SUM -i $REFPROB > ${REFPROB}.refsum
	$LAP_SUM -i $CTGPROB > ${CTGPROB}.ctgsum
fi
