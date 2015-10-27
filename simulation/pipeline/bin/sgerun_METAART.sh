#!/bin/bash
unset module
#
# Create a simulated WGS sequencing run
#

#$ -o logs
#$ -e logs
#$ -cwd
#$ -N METAART_JOB

export PATH=/work/bin:$PATH

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

MA_PATH=$PROXIHOME/simulation/hic_simulator/src/metaART.py
if [ ! -f $MA_PATH ]
then
	echo "$MA_PATH did not exist"
	exit 1
fi

ART_PATH=$PROXIHOME/simulation/pipeline/external/ART/art_illumina
if [ ! -f $ART_PATH ]
then
	echo "$ART_PATH did not exist"
	exit 1
fi

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 9 ]
	then
		echo "Usage: [seed] [insert len] [insert sd] [read len] [x-fold] [comm table] [ref sequences] [base_name] [out_dir]"
		exit 1
	fi
	echo "Submitting run"

	# check for output file existance
	if [ -f $9/${8}1.fq ] && [ -f $9/${8}2.fq ]
	then
		echo "Output files already exist"
		exit 0
	fi

	TARGETS=( "${8}1.fq" "${8}2.fq" )
	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT

	CMD=`readlink -f $0`
	qsub -sync yes -V -b n -v SEED=$1,INSERT_LEN=$2,INSERT_SD=$3,READ_LEN=$4,X_FOLD=$5,COMM_TABLE=$6,REF_SEQ=$7,BASE_NAME=$8,OUT_DIR=$9 $CMD

	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running $MA_PATH --log logs/$PBS_JOBID.ma.log --art-path $ART_PATH -S $SEED -m $INSERT_LEN -s $INSERT_SD -l $READ_LEN -M $X_FOLD -t $COMM_TABLE $REF_SEQ     $OUT_DIR/wgs"

    $MA_PATH --log logs/$PBS_JOBID.ma.log --art-path $ART_PATH -S $SEED -m $INSERT_LEN -s $INSERT_SD -l $READ_LEN -M $X_FOLD -t $COMM_TABLE $REF_SEQ $OUT_DIR/wgs
fi
