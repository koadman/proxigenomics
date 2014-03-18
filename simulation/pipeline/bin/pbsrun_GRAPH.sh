#!/bin/bash

#
# Create graph files
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -N GRAPHJOB

PYTHON=python
#PYTHON=$HOME/python/bin/python
SAMTOEDGES="$PYTHON /panfs/panspermia/120274/work/hi-c/sim/samToEdges.py"

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	
	if [ $# -ne 2 ]
	then
		echo "Usage: [hic2scf] [wgs2scf]"
		exit 1
	fi

	echo "Submitting run"
	qsub -W block=true -v HIC2SCF=$1,WGS2SCF=$2 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	$SAMTOEDGES ${HIC2SCF}.sam ${WGS2SCF}.idxstats edges.csv nodes.csv
	
	# Need weighted edges or more columns in above
	# Need node.csv
fi
