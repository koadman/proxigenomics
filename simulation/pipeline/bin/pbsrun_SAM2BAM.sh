#!/bin/bash

#
# Samtools on simulation
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -N SAMTOOLSJOB

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne  1 ]
	then
		echo "Usage: [sam file]"
		exit 1
	fi
	
	echo "Submitting run"
	qsub -W block=true -v SAMFILE=$1 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	# Make sorted indexed bam and then some stat files
	bf=${SAMFILE%.sam}
	samtools view -bS $SAMFILE | samtools sort - $bf
	samtools index ${bs}.bam
	samtools idxstats ${bf}.bam > ${bf}.idxstats
	samtools flagstat ${bf}.bam > ${bf}.flagstat
	
fi
