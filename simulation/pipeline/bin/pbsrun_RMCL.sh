#!/bin/bash

#
# Map sequences in simulation
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N RMCLJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

CMD=$HOME/bin/rmcl

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	while getopts ":b:i:" opt
	do
		case $opt in
			b)
				CMDOPT="$CMDOPT -b $OPTARG"
				;;
			i)
				CMDOPT="$CMDOPT -i $OPTARG"
				;;
		esac
	done

	shift $(( OPTIND-1 ))

	if [ $# -ne 2 ]
	then
		echo "Usage: [-b balance] [-i inflation] <input> <output>"
		exit 1
	fi
	
	export CMDOPT
	
	echo "Submitting run"
	trap 'rollback_rm_file $2; exit $?' INT TERM EXIT
	qsub -V -W block=true -v INPUT=$1,OUTPUT=$2 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	if [ -s $INPUT ]
	then
		$CMD $CMDOPT -o $OUTPUT $INPUT
	else
		echo "Input had no data" > $OUTPUT
	fi

fi
