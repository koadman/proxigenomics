#!/bin/bash

#
# Create simulated Hi-C sequencing data
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N HICJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

HICEXE=$HOME/git/proxigenomics/simulation/hic_simulator/src/simForward.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 7 ]
	then
		echo "Usage: [seed] [n_frag] [read len] [inter prob] [table] [ref seq] [output fasta]"
		exit 1
	fi

	#trap 'rollback_rm_file $7; exit $?' INT TERM EXIT
	qsub -sync yes -V -v SEED=$1,N_FRAG=$2,READ_LEN=$3,INTER_PROB=$4,TABLE=$5,REF_SEQ=$6,OUTPUT=$7 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

    $HICEXE -r $SEED -n $N_FRAG -l $READ_LEN -p $INTER_PROB -t $TABLE -s $REF_SEQ -o ${OUTPUT}
fi
