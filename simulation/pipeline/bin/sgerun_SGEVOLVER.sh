#!/bin/bash

#
# Create a simulated community
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N SGEVOLVERJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

echo "BINDIR=$BINDIR"

PATH=$PATH:$BINDIR:$BINDIR/sgevolver
SGBIN=$PWD/bin/sgevolver/simujobrun.pl
SGPARMS=$PWD/bin/sgevolver/prepParams.py


if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 7 ]
	then
		echo "Usage: [seed] [tree scale] [sg scale] [length] [tree] [input seq] [output seqs]"
		exit 1
	fi

	if [ ! -d logs ]
	then
	    echo "Log folder is missing, creating one..."
	    mkdir logs
	fi

	echo "Submitting run"

#	trap 'rollback_rm_file $7; exit $?' INT TERM EXIT
	# get the absolute path for output
	OUT_ABS=`readlink -f $7`
	qsub -sync yes -V -v SEED=$1,TR_SCALE="$2",SG_SCALE="$3",LENGTH=$4,TREE=$5,INPUT_SEQ=$6,OUTPUT_SEQ=$OUT_ABS $0
#	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
    echo "Running"

    OUTDIR=`dirname $OUTPUT_SEQ`

    # copy base files
    cp -n $TREE $OUTDIR
    cp -n $INPUT_SEQ $OUTDIR

    cd $OUTDIR

    # create runtime parameter file
    $SGPARMS --tree `basename $TREE` --seq `basename $INPUT_SEQ` --seq-len $LENGTH --tree-scale $TR_SCALE --sg-scale $SG_SCALE

    $SGBIN $INPUT_SEQ $SEED

    cp evolved_seqs.fas `basename $OUTPUT_SEQ`

fi
