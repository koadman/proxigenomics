#!/bin/bash

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

#
# Map sequences in simulation
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N MAPJOB

BWAEXE=$HOME/bin/bwa-0.7.6a/bwa
SAMTOOLS=$HOME/bin/samtools

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ] && [ $# -ne 4 ]
	then
		echo "Usage: <subject> <query> <output>"
		echo "or"
		echo "Usage: <subject> <query1> <query2> <output>"
		exit 1
	fi

	if [ $# -eq 3 ]
	then
    	echo "Submitting single query run"
    	BASE=${3%.sam}
    	TARGET=( ${BASE}.sam ${BASE}.bam ${BASE}.bai ${BASE}.idxstats ${BASE}.flagstats )
    	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	    qsub -W block=true -v SUBJECT=$1,QUERY1=$2,QUERY2="",OUTPUT=$3 $0
    	trap - INT TERM EXIT
    elif [ $# -eq 4 ]
    then
    	echo "Submitting double query run"
    	BASE=${4%.sam}
    	TARGET=( ${BASE}.sam ${BASE}.bam ${BASE}.bai ${BASE}.idxstats ${BASE}.flagstats )
        trap 'rollback_rm_files ${TARGET[@}}; exit $?' INT TERM EXIT
        qsub -W block=true -v SUBJECT=$1,QUERY1=$2,QUERY2=$3,OUTPUT=$4 $0
    	trap - INT TERM EXIT
    fi
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	# Map reads. The second query file can be empty
	$BWAEXE mem -t 2 $SUBJECT $QUERY1 $QUERY2 > $OUTPUT

    BASE=${sf%.sam}
    $SAMTOOLS view -bS $OUTPUT | $SAMTOOLS sort - $BASE
    $SAMTOOLS index ${BASE}.bam
    $SAMTOOLS idxstats ${BASE}.bam > ${BASE}.idxstats
    $SAMTOOLS flagstat ${BASE}.bam > ${BASE}.flagstat

fi
