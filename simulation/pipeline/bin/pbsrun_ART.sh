#!/bin/bash

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

#
# Create a simulated WGS sequencing run
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N ARTJOB

ARTEXE=$HOME/bin/ART/art_illumina

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 7 ]
	then
		echo "Usage: [seed] [insert len] [insert sd] [read len] [x-fold] [ref sequences] [out base]"
		exit 1
	fi
	echo "Submitting run"

	echo "qsub -W block=true -v SEED=$1,INSERT_LEN=$2,INSERT_SD=$3,READ_LEN=$4,X_FOLD=$5,REF_SEQ=$6,OUT_BASE=$7 $0"

    TARGETS=( "${7}1.fq" "${7}2.fq" "${7}1.aln" "${7}2.aln" )
	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	qsub -W block=true -v SEED=$1,INSERT_LEN=$2,INSERT_SD=$3,READ_LEN=$4,X_FOLD=$5,REF_SEQ=$6,OUT_BASE=$7 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	$ARTEXE -p -rs $SEED -m $INSERT_LEN -s $INSERT_SD -l $READ_LEN -f $X_FOLD -i $REF_SEQ -o $OUT_BASE
fi
