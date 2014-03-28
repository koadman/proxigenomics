#!/bin/bash

#
# Create a simulated WGS sequencing run
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
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
	qsub -W block=true -v SEED=$1,INSERT_LEN=$2,INSERT_SD=$3,READ_LEN=$4,X_FOLD=$5,REF_SEQ=$6,OUT_BASE=$7 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	$ARTEXE -p -rs $SEED -m $INSERT_LEN -s $INSERT_SD -l $READ_LEN -f $X_FOLD -i $REF_SEQ -o $OUT_BASE
fi
