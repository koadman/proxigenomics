#!/bin/bash

#
# Create simulated Hi-C sequencing data
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -N HICJOB

HICEXE=$HOME/git/proxigenomics/simulation/hic_simulator/src/simForward.py
PYTHON=python
#PYTHON=$HOME/python/bin/python

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 7 ]
	then
		echo "Usage: [seed] [n_frag] [read len] [inter prob] [table] [ref seq] [out base]"
		exit 1
	fi
	echo "Submitting run"
	qsub -W block=true -v SEED=$1,N_FRAG=$2,READ_LEN=$3,INTER_PROB=$4,TABLE=$5,REF_SEQ=$6,OUT_BASE=$7 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	$PYTHON $HICEXE -r $SEED -n $N_FRAG -l $READ_LEN -p $INTER_PROB -t $TABLE -s $REF_SEQ -o ${OUT_BASE}.fasta
fi
