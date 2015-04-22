#!/bin/bash

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N MCLJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

MCLEXE=$HOME/bin/mcl

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ]
	then
		echo "Usage: <inflation> <input> <output>"
		exit 1
	fi
	echo "Submitting run"

	trap 'rollback_rm_file $3; exit $?' INT TERM EXIT
	qsub -W block=true -v INFLATION=$1,INPUT=$2,OUTPUT=$3 $0
	trap - INT TERM EXIT
	echo "Finished"

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
