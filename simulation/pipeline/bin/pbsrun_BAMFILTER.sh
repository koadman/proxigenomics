#!/bin/bash

#
# Filter bam file
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N FILTERJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

FILTERBAM=bin/filter_bam.py
SAMTOOLS=$HOME/bin/samtools

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then

	if [ $# -ne 4 ]
	then
		echo "Usage: [mincov] [minqual] [in bam] [out bam]"
		exit 1
	fi

	echo "Submitting run"
	trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -W block=true -v MINCOV=$1,MINQUAL=$2,BAMFILE=$3,OUTPUT=$4 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

    BASE=${OUTPUT%.bam}

	$FILTERBAM --mincov $MINCOV --minqual $MINQUAL $BAMFILE ${BASE}.tmp
    $SAMTOOLS sort ${BASE}.tmp $BASE
    $SAMTOOLS index $OUTPUT
    $SAMTOOLS idxstats $OUTPUT > ${BASE}.idxstats
    $SAMTOOLS flagstat $OUTPUT > ${BASE}.flagstat
    rm ${BASE}.tmp

fi
