#!/bin/bash

#
# Create graph files
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N GRAPHJOB

SAMTOEDGES=bin/samToEdges.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	
	if [ $# -ne 4 ]
	then
		echo "Usage: [hic2ctg.sam] [wgs2ctg.bam] [edge out] [node out]"
		exit 1
	fi

	echo "Submitting run"
	qsub -W block=true -v HIC2CTG=$1,WGS2CTG=$2,EDGES=$3,NODES=$4 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	$SAMTOEDGES ${HIC2CTG} ${WGS2CTG%.bam}.idxstats $EDGES $NODES
	
	# Need weighted edges or more columns in above
	# Need node.csv
fi
