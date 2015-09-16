#!/bin/bash

#
# Create phylogenetic reconstruction of multi-fasta
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N TREEJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

MKTREE=bin/mktree.sh

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then

	if [ $# -ne 2 ]
	then
		echo "Usage: [input multifasta] [output dir]"
		exit 1
	fi

	echo "Submitting run"
	trap 'rollback_rm_file $7/genomes.xmfa.tre; exit $?' INT TERM EXIT
	qsub -W block=true -v INPUT_FASTA=$1,OUTPUT_DIR=$2 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

    if [ ! -e $INPUT ]
    then
       sleep 5
    fi

    $MKTREE $INPUT_FASTA $OUTPUT_DIR

fi
