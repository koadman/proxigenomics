#!/bin/bash

#
# Create a simulated WGS sequencing run
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -e logs/
#PBS -o logs/
#PBS -N METAARTJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

MA_PATH=$HOME/git/proxigenomics/simulation/hic_simulator/src/metaART.py
ART_PATH=$HOME/bin/ART/art_illumina

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 8 ]
	then
		echo "Usage: [seed] [insert len] [insert sd] [read len] [x-fold] [comm table] [ref sequences] [out base]"
		exit 1
	fi
	echo "Submitting run"

    TARGETS=( "${8}1.fq" "${8}2.fq" )
	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	qsub -W block=true -v SEED=$1,INSERT_LEN=$2,INSERT_SD=$3,READ_LEN=$4,X_FOLD=$5,COMM_TABLE=$6,REF_SEQ=$7,OUT_BASE=$8 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	$MA_PATH --log logs/$PBS_JOBID.ma.log --art-path $ART_PATH -S $SEED -m $INSERT_LEN -s $INSERT_SD -l $READ_LEN -M $X_FOLD -t $COMM_TABLE $REF_SEQ $OUT_BASE
fi
