#!/bin/bash

#
# Filter bam file
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N FILTERJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

FILTERBAM=bin/filter_bam.py
SAMTOOLS=external/samtools

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

	if [ $# -ne 4 ]
	then
		echo "Usage: [mincov] [minqual] [in bam] [out bam]"
		exit 1
	fi

	echo "Submitting run"
	#trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	CMD=`readlink -f $0`
	qsub -sync yes -V -v MINCOV=$1,MINQUAL=$2,BAMFILE=$3,OUTPUT=$4 $CMD
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

    BASE=${OUTPUT%.bam}

	$FILTERBAM --mincov $MINCOV --minqual $MINQUAL $BAMFILE ${BASE}.tmp
    $SAMTOOLS sort ${BASE}.tmp $BASE
    $SAMTOOLS index $OUTPUT
    $SAMTOOLS idxstats $OUTPUT > ${BASE}.idxstats
    $SAMTOOLS flagstat $OUTPUT > ${BASE}.flagstat
    rm ${BASE}.tmp

fi
