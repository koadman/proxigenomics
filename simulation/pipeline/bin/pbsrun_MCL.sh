#!/bin/bash

#
# Map sequences in simulation
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -N MCLJOB

MCLEXE=$HOME/bin/mcl

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ]
	then
		echo "Usage: <inflation> <input> <output>"
		exit 1
	fi
	echo "Submitting run"
	qsub -W block=true -v INFLATION=$1,INPUT=$2,OUTPUT=$3 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	if [ -s $INPUT ]
	then
		$MCLEXE $INPUT -te 2 --abc -I $INFLATION -o $OUTPUT
	else
		echo "Input had no data" > $OUTPUT
	fi

fi
