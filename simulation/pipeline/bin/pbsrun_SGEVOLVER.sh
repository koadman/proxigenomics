#!/bin/bash

#
# Create a simulated community
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N SGEVOLVERJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $HICPIPE/bin/bash_init.sh
fi

PATH=$PATH:$HICPIPE/bin:$HICPIPE/bin/sgevolver
SGBIN=$HICPIPE/bin/sgevolver/simujobrun.pl
SGPARMS=$HICPIPE/bin/sgevolver/prepParams.py


if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 6 ]
	then
		echo "Usage: [seed] [scale] [length] [tree] [input seq] [output seqs]"
		exit 1
	fi
	echo "Submitting run"

	trap 'rollback_rm_file $6; exit $?' INT TERM EXIT
	qsub -W block=true -v SEED=$1,SCALE="$2",LENGTH=$3,TREE=$4,INPUT_SEQ=$5,OUTPUT_SEQ=$6 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
    echo "Running"
    cd $PBS_O_WORKDIR

    OUTDIR=`dirname $OUTPUT_SEQ`

    # copy base files
    cp $TREE $OUTDIR
    cp $INPUT_SEQ $OUTDIR
    cd $OUTDIR

    # create runtime parameter file
    $SGPARMS --tree `basename $TREE` --seq `basename $INPUT_SEQ` --seq-len $LENGTH --tree-scale $SCALE -o .

    $SGBIN $INPUT_SEQ $SEED
    cp evolved_seqs.fas `basename $OUTPUT_SEQ`

fi
